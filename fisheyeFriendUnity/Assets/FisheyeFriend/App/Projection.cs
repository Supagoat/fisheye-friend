using System;
using System.IO;
using System.Threading.Tasks;
using FisheyeFriend.Core.Defish;

namespace FisheyeFriend.App
{
    /// <summary>
    /// The defish projection, loaded from the point-pair file and prepared on a worker thread at start-up (fitting the
    /// splines and evaluating them on the fine grid takes a second or so), so the first photo does not wait for it.
    /// The state the interface sees changes only in <see cref="Tick"/> (called from Update), never between the layout
    /// and repaint passes of one GUI event, so the controls drawn for it stay the same within a frame.
    /// </summary>
    internal sealed class Projection
    {
        private readonly Task<DefishGeometry> _task;
        private DefishGeometry _current;
        private bool _failed;

        public string Name { get; }
        public string Path { get; }

        public Projection(AppConfig cfg)
        {
            Name = cfg.projectionName;
            Path = cfg.ProjectionPath;
            string path = Path;
            _task = Task.Run(() =>
            {
                if (!File.Exists(path)) throw new FileNotFoundException("The projection file was not found: " + path);
                var mapping = HemiMapping.Parse(File.ReadAllText(path));
                return new DefishGeometry(mapping);
            });
        }

        /// <summary>Takes in the worker's progress; call once per frame, outside OnGUI.</summary>
        public void Tick()
        {
            if (_current != null || _failed) return;
            if (_task.Status == TaskStatus.RanToCompletion) _current = _task.Result;
            else if (_task.IsFaulted || _task.IsCanceled) _failed = true;
        }

        public bool Ready => _current != null;
        public bool Failed => _failed;

        /// <summary>Why the projection could not be prepared, or null.</summary>
        public string Error
        {
            get
            {
                if (!Failed) return null;
                var e = _task.Exception != null ? _task.Exception.GetBaseException() : null;
                return e != null ? e.Message : "cancelled";
            }
        }

        /// <summary>The geometry once prepared (throws if it could not be).</summary>
        public Task<DefishGeometry> GetAsync()
        {
            return _task;
        }

        /// <summary>The geometry if it is ready, else null.</summary>
        public DefishGeometry Current => _current;

        public string Describe()
        {
            if (Failed) return Name + " - could not load: " + Error;
            if (!Ready) return Name + " - preparing…";
            var g = _current;
            return Name + " - " + g.Mapping.Pairs.Count + " point pairs, measured on " + g.Mapping.SourceWidth + " x " + g.Mapping.SourceHeight;
        }
    }
}
