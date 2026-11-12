using System;
using System.Globalization;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;
using UnityEngine;

namespace FisheyeFriend.App.UI
{
    /// <summary>
    /// The defish controls shared by interactive and batch modes: the projection, the output size, what happens at the
    /// frame's edges, the resampling filter and the advanced options.
    /// </summary>
    internal sealed class SettingsPanel
    {
        private struct FilterInfo
        {
            public WarpKernel Kernel;
            public double Radius;
            public string Name, Description;

            public FilterInfo(WarpKernel k, double r, string name, string description)
            {
                Kernel = k;
                Radius = r;
                Name = name;
                Description = description;
            }
        }

        private static readonly FilterInfo[] Filters =
        {
            new FilterInfo(WarpKernel.TensorLanczos, 4, "Lanczos 4",
                "Separable 4-lobe Lanczos: measured best against ground truth, 1.9-4.8 dB above the old renderer's radial filter. Recommended."),
            new FilterInfo(WarpKernel.TensorLanczos, 3, "Lanczos 3",
                "Separable 3-lobe Lanczos: 0.02-0.05 dB behind Lanczos 4, about 25% faster."),
            new FilterInfo(WarpKernel.Bicubic, 2, "Bicubic",
                "Catmull-Rom bicubic: fastest, 0.5-0.8 dB softer."),
        };

        private static readonly string[] LightNames = { "Linear", "Gamma", "Sigmoidal" };
        private static readonly LightSpace[] LightValues = { LightSpace.Linear, LightSpace.Gamma, LightSpace.Sigmoidal };

        private readonly AppConfig _cfg;
        private readonly Action<Action> _defer;
        private readonly Action _onChanged;
        private bool _advancedOpen;

        // The planned size is cached: the crop search maps a few thousand points.
        private int _sizeW = -1, _sizeH = -1, _sizeO = -1;
        private DefishOutputSize _sizeFor;
        private DefishEdges _edgesFor;
        private DefishGeometry _sizeGeometry;
        private int _outW, _outH;

        /// <param name="defer">Runs structural changes outside OnGUI (see FisheyeFriendApp.Defer).</param>
        /// <param name="onChanged">Called after a deferred change to the settings.</param>
        public SettingsPanel(AppConfig cfg, Action<Action> defer, Action onChanged)
        {
            _cfg = cfg;
            _defer = defer;
            _onChanged = onChanged;
        }

        public static string FilterName(DefishSettings s)
        {
            foreach (var f in Filters)
                if (f.Kernel == s.Kernel && Math.Abs(f.Radius - s.Radius) < 1e-9) return f.Name;
            return s.Kernel + " " + s.Radius.ToString("0.#", CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Draws the controls (GUILayout). Returns true if any setting changed. <paramref name="image"/> is the open
        /// (upright) image in interactive mode, or null; <paramref name="projection"/> is the projection being prepared.
        /// </summary>
        public bool Draw(DefishSettings s, UiStyles ui, PixelImage image, Projection projection)
        {
            bool changed = false;
            var g = projection.Current;

            // ---------------------------------------------------------------- projection
            GUILayout.Label("PROJECTION", ui.Section);
            GUILayout.Label(projection.Describe(), projection.Failed ? ui.Error : ui.MutedSmall);
            if (image != null && g != null && Defisher.ShapeDiffers(g, image))
                GUILayout.Label("This photo is not 3:2 like the frame the projection was measured on: it is defished stretched to " +
                                "that frame, so the result may be slightly off.", ui.Error);

            // ---------------------------------------------------------------- output size
            GUILayout.Label("OUTPUT SIZE", ui.Section);
            GUILayout.BeginHorizontal();
            changed |= Choice(ui, s.OutputSize == DefishOutputSize.SameAsPhoto, "Same as photo",
                "The result has the photo's pixel size. The projection squeezes the centre of the frame by about 1.19x, so some of its finest detail is averaged away.",
                () => s.OutputSize = DefishOutputSize.SameAsPhoto);
            changed |= Choice(ui, s.OutputSize == DefishOutputSize.FullDetail, "Full detail",
                "About 1.19x larger each way, so the centre of the frame keeps every pixel of the photo (about 1.42x the pixels).",
                () => s.OutputSize = DefishOutputSize.FullDetail);
            GUILayout.EndHorizontal();
            if (image != null && g != null)
            {
                PlannedSize(g, image, s);
                GUILayout.Label(image.Width + " x " + image.Height + "  »  " + _outW + " x " + _outH + " px", ui.Label);
                GUILayout.Label(((double)_outW * _outH / 1e6).ToString("0.0", CultureInfo.InvariantCulture) + " MP  ·  " +
                                ((double)_outW * _outH / ((double)image.Width * image.Height) * 100).ToString("0", CultureInfo.InvariantCulture) + "% of the photo's pixels", ui.MutedSmall);
            }

            // ---------------------------------------------------------------- edges
            GUILayout.Label("EDGES", ui.Section);
            GUILayout.BeginHorizontal();
            changed |= Choice(ui, s.Edges == DefishEdges.Extend, "Full frame",
                "Keep the whole defished frame. Along the top and bottom edges the projection reaches a few pixels past the photo (up to about 5 on a 42 MP frame); those repeat its edge, as Fisheye Friend always did.",
                () => s.Edges = DefishEdges.Extend);
            changed |= Choice(ui, s.Edges == DefishEdges.Crop, "Crop to photo",
                "Crop very slightly (about 0.3% each way), keeping the frame's shape, so every pixel comes from the photo.",
                () => s.Edges = DefishEdges.Crop);
            GUILayout.EndHorizontal();

            // ---------------------------------------------------------------- filter
            GUILayout.Label("FILTER", ui.Section);
            GUILayout.BeginHorizontal();
            string description = "";
            foreach (var f in Filters)
            {
                bool sel = f.Kernel == s.Kernel && Math.Abs(f.Radius - s.Radius) < 1e-9;
                if (sel) description = f.Description;
                var fi = f;
                changed |= Choice(ui, sel, f.Name, f.Description, () =>
                {
                    s.Kernel = fi.Kernel;
                    s.Radius = fi.Radius;
                });
            }
            GUILayout.EndHorizontal();
            if (description.Length > 0) GUILayout.Label(description, ui.MutedSmall);

            // ---------------------------------------------------------------- advanced
            if (GUILayout.Button(_advancedOpen ? "Advanced  ^" : "Advanced  v", ui.Button))
                _defer(() => _advancedOpen = !_advancedOpen);
            if (_advancedOpen)
            {
                GUILayout.Label(new GUIContent("Light space", "Where the filtering happens. Sigmoidal measured best (fewest halos where the projection enlarges)."), ui.Label);
                GUILayout.BeginHorizontal();
                for (int i = 0; i < LightValues.Length; i++)
                {
                    var l = LightValues[i];
                    changed |= Choice(ui, s.Light == l, LightNames[i], "Filter in " + LightNames[i].ToLowerInvariant() + " light", () => s.Light = l);
                }
                GUILayout.EndHorizontal();

                bool ap = ui.Toggle(s.ApertureCorrection, new GUIContent("Aperture correction",
                    "Towards the edges the projection enlarges the photo (up to about 1.8x). There each new pixel covers less of the scene than a " +
                    "photo pixel did, so the photo's own pixel blur is partly undone. Measured +0.1 to +0.3 dB."));
                if (ap != s.ApertureCorrection)
                {
                    s.ApertureCorrection = ap;
                    changed = true;
                }
                bool aa = ui.Toggle(s.Antialias, new GUIContent("Anti-alias the centre",
                    "Widen the filter where the projection shrinks the photo (the centre, by up to 1.19x). Measured 0.1-0.2 dB worse: at so mild a " +
                    "shrink the photo's own blur already prevents aliasing."));
                if (aa != s.Antialias)
                {
                    s.Antialias = aa;
                    changed = true;
                }

                if (GUILayout.Button(new GUIContent("Reset to defaults", "Restore the defaults from " + AppConfig.FileName), ui.Button))
                {
                    _defer(() =>
                    {
                        var d = _cfg.CreateDefaultSettings();
                        s.OutputSize = d.OutputSize;
                        s.Edges = d.Edges;
                        s.Kernel = d.Kernel;
                        s.Radius = d.Radius;
                        s.Light = d.Light;
                        s.Antialias = d.Antialias;
                        s.ApertureCorrection = d.ApertureCorrection;
                        _onChanged();
                    });
                }
            }
            return changed;
        }

        private static bool Choice(UiStyles ui, bool selected, string label, string tip, Action select)
        {
            if (GUILayout.Button(new GUIContent(label, tip), selected ? ui.ButtonSelected : ui.Button) && !selected)
            {
                select();
                return true;
            }
            return false;
        }

        private void PlannedSize(DefishGeometry g, PixelImage image, DefishSettings s)
        {
            int o = image.Metadata != null ? image.Metadata.SourceOrientation : 0;
            if (_sizeW == image.Width && _sizeH == image.Height && _sizeO == o && _sizeFor == s.OutputSize && _edgesFor == s.Edges && _sizeGeometry == g)
                return;
            Defisher.OutputSize(g, image, s, out _outW, out _outH);
            _sizeW = image.Width;
            _sizeH = image.Height;
            _sizeO = o;
            _sizeFor = s.OutputSize;
            _edgesFor = s.Edges;
            _sizeGeometry = g;
        }
    }
}
