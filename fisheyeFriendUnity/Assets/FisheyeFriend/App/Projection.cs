using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading.Tasks;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;

namespace FisheyeFriend.App
{
    /// <summary>
    /// The defish projection and the geometries made from it.
    ///
    /// The projection file is read on a worker thread at start-up. It is normally a model (defish-hemi-model.json), which
    /// the user can adjust; a geometry is then made for each frame shape and adjustment (about half a second, on a
    /// worker; the last few are kept). A file of point pairs (Fisheye Friend's original defish-hemi-596.json) still
    /// works, as fixed splines without adjustment. With no file at all, the built-in model is used.
    ///
    /// The state the interface sees changes only in <see cref="Tick"/> (called from Update), never between the layout
    /// and repaint passes of one GUI event, so the controls drawn for it stay the same within a frame.
    /// </summary>
    internal sealed class Projection
    {
        private sealed class Source
        {
            public DefishModel Model;
            public DefishGeometry Splines;
        }

        private readonly Task<Source> _load;
        private Source _source;
        private bool _failed;

        // geometries: being made, and those the interface may see (updated in Tick)
        private readonly Dictionary<string, Task<DefishGeometry>> _building = new Dictionary<string, Task<DefishGeometry>>();
        private readonly List<string> _order = new List<string>();
        private readonly Dictionary<string, DefishGeometry> _visible = new Dictionary<string, DefishGeometry>();
        private const int Keep = 4;

        public string Name { get; }
        public string Path { get; }

        public Projection(AppConfig cfg)
        {
            Name = cfg.projectionName;
            Path = cfg.ProjectionPath;
            string path = Path;
            _load = Task.Run(() =>
            {
                if (!File.Exists(path)) return new Source { Model = DefishModel.Hemi };
                string text = File.ReadAllText(path);
                var root = Json.Parse(text);
                if (Json.Get(root, "pairs") != null) return new Source { Splines = new DefishGeometry(HemiMapping.Parse(text)) };
                return new Source { Model = DefishModel.Parse(text) };
            });
        }

        /// <summary>Takes in the workers' progress; call once per frame, outside OnGUI.</summary>
        public void Tick()
        {
            if (_source == null && !_failed)
            {
                if (_load.Status == TaskStatus.RanToCompletion) _source = _load.Result;
                else if (_load.IsFaulted || _load.IsCanceled) _failed = true;
            }
            foreach (var kv in _building)
                if (kv.Value.Status == TaskStatus.RanToCompletion && !_visible.ContainsKey(kv.Key)) _visible[kv.Key] = kv.Value.Result;
        }

        public bool Ready => _source != null;
        public bool Failed => _failed;
        /// <summary>True when the projection is a model the user can adjust (not fixed point pairs).</summary>
        public bool Adjustable => _source != null && _source.Model != null;
        public DefishModel Model => _source != null ? _source.Model : null;

        /// <summary>Why the projection could not be prepared, or null.</summary>
        public string Error
        {
            get
            {
                if (!Failed) return null;
                var e = _load.Exception != null ? _load.Exception.GetBaseException() : null;
                return e != null ? e.Message : "cancelled";
            }
        }

        private static string KeyFor(double aspect, DefishAdjustment adj, int nodes)
        {
            return aspect.ToString("0.00000", CultureInfo.InvariantCulture) + "|" + nodes + "|" + (adj ?? new DefishAdjustment()).Key;
        }

        /// <summary>
        /// The geometry for frames of height / width <paramref name="aspect"/> (as the lens saw them) with an adjustment,
        /// made on a worker if it isn't ready. <paramref name="preview"/> gives a coarse one in a few milliseconds.
        /// Call on the main thread.
        /// </summary>
        public async Task<DefishGeometry> GeometryAsync(double aspect, DefishAdjustment adjustment, bool preview = false)
        {
            var src = await _load;
            if (src.Splines != null) return src.Splines;
            int nodes = preview ? 129 : 1025;
            var adj = (adjustment ?? new DefishAdjustment()).Clone();
            adj.Validate();
            if (preview)
            {
                var m = src.Model;
                return await Task.Run(() => DefishGeometry.FromModel(m, adj, aspect, nodes));    // a few ms: not kept
            }
            string key = KeyFor(aspect, adj, nodes);
            Task<DefishGeometry> task;
            if (!_building.TryGetValue(key, out task))
            {
                var model = src.Model;
                task = Task.Run(() => DefishGeometry.FromModel(model, adj, aspect, nodes));
                _building[key] = task;
                _order.Add(key);
                Evict();
            }
            else
            {
                _order.Remove(key);
                _order.Add(key);
            }
            return await task;
        }

        /// <summary>Forgets the least recently asked-for geometries that are finished, beyond the last few.</summary>
        private void Evict()
        {
            for (int i = 0; i < _order.Count && _order.Count > Keep; )
            {
                Task<DefishGeometry> t;
                if (_building.TryGetValue(_order[i], out t) && t.IsCompleted)
                {
                    _building.Remove(_order[i]);
                    _visible.Remove(_order[i]);
                    _order.RemoveAt(i);
                }
                else i++;
            }
        }

        /// <summary>
        /// The full geometry if it is ready (as of the last <see cref="Tick"/>), else null; asks for it to be made.
        /// Safe to call from OnGUI.
        /// </summary>
        public DefishGeometry TryGeometry(double aspect, DefishAdjustment adjustment)
        {
            if (_source == null) return null;
            if (_source.Splines != null) return _source.Splines;
            var adj = (adjustment ?? new DefishAdjustment()).Clone();
            adj.Validate();
            string key = KeyFor(aspect, adj, 1025);
            DefishGeometry g;
            if (_visible.TryGetValue(key, out g)) return g;
            if (!_building.ContainsKey(key)) Request(aspect, adj);
            return null;
        }

        /// <summary>The full geometry if it is ready (as of the last <see cref="Tick"/>), else null; asks for nothing.</summary>
        public DefishGeometry PeekGeometry(double aspect, DefishAdjustment adjustment)
        {
            if (_source == null) return null;
            if (_source.Splines != null) return _source.Splines;
            var adj = (adjustment ?? new DefishAdjustment()).Clone();
            adj.Validate();
            DefishGeometry g;
            return _visible.TryGetValue(KeyFor(aspect, adj, 1025), out g) ? g : null;
        }

        /// <summary>
        /// A coarse geometry made now, on the calling thread (a few milliseconds), for previews; the fixed splines when
        /// the projection is point pairs; null until the projection is loaded.
        /// </summary>
        public DefishGeometry CoarseGeometry(double aspect, DefishAdjustment adjustment)
        {
            var src = _source;
            if (src == null) return null;
            if (src.Splines != null) return src.Splines;
            return DefishGeometry.FromModel(src.Model, adjustment ?? new DefishAdjustment(), aspect, 129);
        }

        private async void Request(double aspect, DefishAdjustment adj)
        {
            try
            {
                await GeometryAsync(aspect, adj);
            }
            catch (Exception)
            {
                // reported when the geometry is used
            }
        }

        public string Describe()
        {
            if (Failed) return Name + " - could not load: " + Error;
            if (!Ready) return Name + " - preparing…";
            if (_source.Splines != null)
            {
                var m = _source.Splines.Mapping;
                return Name + " - " + m.Pairs.Count + " point pairs (fixed: not adjustable)";
            }
            return _source.Model.Name + " model - " + _source.Model.DiagonalFieldOfView(2.0 / 3.0).ToString("0", CultureInfo.InvariantCulture) +
                   "° diagonal fisheye, Panini-style projection";
        }
    }
}
