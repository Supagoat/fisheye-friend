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
    /// The defish controls shared by interactive and batch modes: the projection and its adjustment (with a live
    /// preview), the output size, what happens at the frame's edges, the resampling filter and the advanced options.
    /// </summary>
    internal sealed class SettingsPanel : IDisposable
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

        private static readonly Color ColorX = UiStyles.Accent;
        private static readonly Color ColorY = UiStyles.Hex(0x62A8E5);
        private static readonly int CurveHint = "FisheyeFriend.DefishCurve".GetHashCode();

        /// <summary>Seconds an adjustment must stay put before the full-quality geometry is made for it.</summary>
        private const float SettleSeconds = 0.6f;

        private readonly AppConfig _cfg;
        private readonly Action<Action> _defer;
        private readonly Action _onChanged;
        private readonly LivePreview _preview = new LivePreview();
        private bool _advancedOpen;
        private bool _linked;
        private bool _editY;                    // unlinked: which curve the graph edits

        // state taken in Tick (from Update), so the controls drawn stay the same within a GUI event
        private string _adjustKey;
        private float _adjustChangedAt;
        private double _aspect = 2.0 / 3.0;
        private DefishGeometry _geometry;       // for the photo and the current adjustment: the full one, else the preview's
        private bool _geometryApprox;           // _geometry is the preview's coarse one
        private bool _folds;

        // the preview's size, measured when it was last drawn
        private int _thumbW, _thumbH;
        private float _thumbBoxH = -1;

        // curve dragging
        private int _dragIndex = -1;
        private float _dragOffset;

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
            _linked = Prefs.GetBool("DefishLinked", true);
        }

        public static string FilterName(DefishSettings s)
        {
            foreach (var f in Filters)
                if (f.Kernel == s.Kernel && Math.Abs(f.Radius - s.Radius) < 1e-9) return f.Name;
            return s.Kernel + " " + s.Radius.ToString("0.#", CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Follows the settings and the open photo: asks for the geometry and the preview they need. Call once per frame
        /// from Update, after <see cref="Projection.Tick"/>.
        /// </summary>
        public void Tick(DefishSettings s, PixelImage image, Projection projection)
        {
            double aspect = image != null ? Defisher.SensorAspect(image) : 2.0 / 3.0;
            string key = LivePreview.KeyFor(aspect, s.Adjust);
            float now = Time.realtimeSinceStartup;
            if (key != _adjustKey)
            {
                _adjustKey = key;
                _adjustChangedAt = now;
            }
            // While a slider or the curve is being dragged, only the quick preview follows; the full geometry (about half
            // a second) is made once the adjustment has stayed put for a moment.
            bool settled = now - _adjustChangedAt > SettleSeconds && GUIUtility.hotControl == 0;
            var full = projection.PeekGeometry(aspect, s.Adjust);
            if (full == null && settled) full = projection.TryGeometry(aspect, s.Adjust);

            _preview.Want(image, aspect, s.Adjust, _thumbW, _thumbH);
            _preview.Tick(projection);

            var coarse = _preview.GeometryKey == key ? _preview.Geometry : null;
            _geometry = full ?? coarse;
            _geometryApprox = full == null && coarse != null;
            _folds = _geometry != null && _geometry.Folds;
            _aspect = aspect;
        }

        /// <summary>
        /// Draws the controls (GUILayout). Returns true if any setting changed. <paramref name="image"/> is the open
        /// (upright) image in interactive mode, or null; <paramref name="projection"/> is the projection being prepared.
        /// </summary>
        public bool Draw(DefishSettings s, UiStyles ui, PixelImage image, Projection projection)
        {
            bool changed = false;
            // What decides which controls are drawn is taken here, before any control can change a setting: IMGUI needs the
            // same controls in the layout and the event pass.
            var g = _geometry;
            bool approx = _geometryApprox;
            bool folds = _folds;
            bool linked = _linked;
            bool adjustable = !projection.Ready || projection.Adjustable;
            bool fixedPairs = projection.Ready && !projection.Adjustable;
            bool sameSize = s.OutputSize == DefishOutputSize.SameAsPhoto;
            bool showSize = image != null && g != null;
            bool shapeDiffers = image != null && g != null && Defisher.ShapeDiffers(g, image);
            var model = projection.Model ?? DefishModel.Hemi;
            double centre = g != null ? Defisher.CentreScale(g, 1500, (int)Math.Round(1500 * _aspect)) : 1.19;
            string centreText = centre.ToString("0.00", CultureInfo.InvariantCulture) + "x";

            // ---------------------------------------------------------------- projection
            GUILayout.Label("PROJECTION", ui.Section);
            GUILayout.Label(projection.Describe(), projection.Failed ? ui.Error : ui.MutedSmall);
            if (shapeDiffers)
                GUILayout.Label("This photo is not 3:2 like the frame the projection was measured on: it is defished stretched to " +
                                "that frame, so the result may be slightly off.", ui.Error);

            // ---------------------------------------------------------------- preview
            if (image != null) DrawPreview(ui, image);

            // ---------------------------------------------------------------- defish adjustment
            GUILayout.Label("DEFISH", ui.Section);
            if (fixedPairs)
                GUILayout.Label("The projection file is a fixed set of point pairs, which can't be adjusted. Set projectionFile to " +
                                "defish-hemi-model.json in " + AppConfig.FileName + " to adjust the defish.", ui.MutedSmall);
            bool wasEnabled = GUI.enabled;
            GUI.enabled = wasEnabled && adjustable;
            changed |= DrawAdjustment(s, ui, image, model, linked, folds);
            GUI.enabled = wasEnabled;
            if (image == null)
                GUILayout.Label("Open a photo (in Interactive mode) to preview these settings on it as you change them.", ui.MutedSmall);

            // ---------------------------------------------------------------- output size
            GUILayout.Label("OUTPUT SIZE", ui.Section);
            GUILayout.BeginHorizontal();
            changed |= Choice(ui, s.OutputSize == DefishOutputSize.FullDetail, "Full detail",
                "Recommended. Enlarged so the centre of the frame isn't reduced (with these settings " + centreText + " each way) and keeps the " +
                "photo's sharpness: measured within 1% of the original for all but the very finest detail.",
                () => s.OutputSize = DefishOutputSize.FullDetail);
            changed |= Choice(ui, s.OutputSize == DefishOutputSize.SameAsPhoto, "Same as photo",
                "The photo's pixel size, as the old Fisheye Friend made it. The projection then reduces the centre of the frame (with these " +
                "settings " + centreText + "), which loses its finest detail: at the same scale it looks softer than the original.",
                () => s.OutputSize = DefishOutputSize.SameAsPhoto);
            GUILayout.EndHorizontal();
            if (sameSize)
                GUILayout.Label("The centre is reduced " + centreText + " and loses its finest detail. Full detail keeps the photo's sharpness.", ui.MutedSmall);
            if (showSize)
            {
                PlannedSize(g, image, s);
                string about = approx ? "≈ " : "";
                GUILayout.Label(image.Width + " x " + image.Height + "  »  " + about + _outW + " x " + _outH + " px", ui.Label);
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
                "Crop slightly (about 0.3% each way at the default defish), keeping the frame's shape, so every pixel comes from the photo.",
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

                if (GUILayout.Button(new GUIContent("Reset to defaults", "Restore the defaults from " + AppConfig.FileName +
                                                                         " (the defish adjustment has its own Reset, above)"), ui.Button))
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

        // ============================================================================ preview

        private void DrawPreview(UiStyles ui, PixelImage image)
        {
            float boxH = _thumbBoxH > 0 ? _thumbBoxH : ui.S(200);
            Rect r = GUILayoutUtility.GetRect(ui.S(10), boxH, GUILayout.ExpandWidth(true));
            if (Event.current.type == EventType.Repaint)
            {
                // the picture has the photo's shape (the defish keeps it); portrait photos are kept from growing too tall
                double ratio = (double)image.Height / image.Width;
                float w = r.width, h = (float)(w * ratio), maxH = r.width * 0.75f;
                if (h > maxH)
                {
                    h = maxH;
                    w = (float)(h / ratio);
                }
                var pr = new Rect(Mathf.Round(r.x + (r.width - w) * 0.5f), Mathf.Round(r.y), Mathf.Round(w), Mathf.Round(h));
                _thumbW = (int)pr.width;
                _thumbH = (int)pr.height;
                _thumbBoxH = h;
                ui.Fill(pr, UiStyles.ViewerBackground);
                if (_preview.HasPicture && _preview.Texture != null) GUI.DrawTexture(pr, _preview.Texture, ScaleMode.StretchToFill);
                else GUI.Label(pr, "Preparing the preview…", new GUIStyle(ui.MutedSmall) { alignment = TextAnchor.MiddleCenter });
            }
            GUILayout.Label(new GUIContent("Preview of the whole defished frame. Defish (Ctrl/Cmd+Enter) to see it at full size.",
                "A quick, reduced rendering that follows the settings as you change them. The full-size result uses the same projection " +
                "and the filter chosen below."), ui.MutedSmall);
            if (_preview.Error != null) GUILayout.Label(_preview.Error, ui.Error);
        }

        // ============================================================================ adjustment

        private bool DrawAdjustment(DefishSettings s, UiStyles ui, PixelImage image, DefishModel model, bool linked, bool folds)
        {
            bool changed = false;
            var a = s.Adjust;
            bool editY = !linked && _editY;

            // ---- strength
            double fov = model.DiagonalFieldOfView(_aspect, a.Strength);
            const string strengthTip =
                "How strongly your fisheye lens bends the picture. 100% is the lens this defish was measured on (177° across the diagonal of a " +
                "3:2 frame); lower it for a milder fisheye, raise it for a stronger one. If vertical lines near the sides still bow outwards " +
                "after defishing, raise it; if they bow inwards, lower it. The result keeps the photo's width either way.";
            double strength = Slider(ui, new GUIContent("Lens strength", strengthTip),
                Pct(a.Strength) + "  ·  " + fov.ToString("0", CultureInfo.InvariantCulture) + "°", a.Strength,
                DefishAdjustment.MinStrength, DefishAdjustment.MaxStrength, 1, 120);
            if (strength != a.Strength)
            {
                a.Strength = strength;
                changed = true;
            }

            // ---- amounts
            const string xTip =
                "How far the defish moves pixels sideways: this straightens vertical lines and takes away the fisheye's sideways stretch towards " +
                "the edges. 0% keeps the photo's left-right positions; past 100% it goes further than the defish was measured to.";
            const string yTip =
                "How far the defish moves pixels up and down: this evens out the vertical stretch and sets how much horizontal lines curve. " +
                "0% keeps the photo's up-down positions; past 100% it goes further than the defish was measured to.";
            double ax = Slider(ui, new GUIContent("Horizontal", xTip), Pct(a.AmountX), a.AmountX, 0, DefishAdjustment.MaxAmount, 1, 64);
            if (ax != a.AmountX)
            {
                a.AmountX = ax;
                if (linked) a.AmountY = ax;
                changed = true;
            }
            double ay = Slider(ui, new GUIContent("Vertical", yTip), Pct(a.AmountY), a.AmountY, 0, DefishAdjustment.MaxAmount, 1, 64);
            if (ay != a.AmountY)
            {
                a.AmountY = ay;
                if (linked) a.AmountX = ay;
                changed = true;
            }
            bool link = ui.Toggle(linked, new GUIContent("Link horizontal and vertical",
                "Move both amounts, and both curves below, together. Linking copies the horizontal settings to the vertical ones."));
            if (link != linked)
            {
                _defer(() =>
                {
                    _linked = link;
                    Prefs.SetBool("DefishLinked", link);
                    if (link)
                    {
                        s.Adjust.AmountY = s.Adjust.AmountX;
                        s.Adjust.ShapeY = (double[])s.Adjust.ShapeX.Clone();
                    }
                    _onChanged();
                });
            }

            // ---- across the frame
            const string curveTip =
                "How much of the correction to apply at each distance from the centre, as a share of the amounts above: the centre, a third " +
                "and two thirds of the way out, and the corners. Drag a point up or down (Alt/Option for exact values, without snapping to " +
                "100%); double-click a point for 100%. Lowering the outer points keeps more of the fisheye look towards the edges.";
            GUILayout.BeginHorizontal();
            GUILayout.Label(new GUIContent("Across the frame", curveTip), ui.Label);
            GUILayout.FlexibleSpace();
            if (!linked)
            {
                if (GUILayout.Button(new GUIContent("Horizontal", "Edit the horizontal curve"), !editY ? ui.ButtonSelected : ui.Button)) _editY = false;
                if (GUILayout.Button(new GUIContent("Vertical", "Edit the vertical curve"), editY ? ui.ButtonSelected : ui.Button)) _editY = true;
            }
            GUILayout.EndHorizontal();
            changed |= DrawCurves(ui, a, linked, editY, image);
            var edit = linked || !editY ? a.ShapeX : a.ShapeY;
            GUILayout.Label((linked ? "Both" : editY ? "Vertical" : "Horizontal") + ":  " + Pct(edit[0]) + "  ·  " + Pct(edit[1]) + "  ·  " +
                            Pct(edit[2]) + "  ·  " + Pct(edit[3]), ui.MutedSmall);

            GUILayout.BeginHorizontal();
            if (GUILayout.Button(new GUIContent("Flatten curves", "Set every point of both curves to 100%: the amounts above apply evenly across the frame"), ui.Button) &&
                !a.ShapesFlat)
            {
                a.ShapeX = new[] { 1.0, 1, 1, 1 };
                a.ShapeY = new[] { 1.0, 1, 1, 1 };
                changed = true;
            }
            if (GUILayout.Button(new GUIContent("Reset defish", "Lens strength, amounts and curves back to the defaults from " + AppConfig.FileName), ui.Button))
            {
                _defer(() =>
                {
                    s.Adjust = _cfg.CreateDefaultSettings().Adjust.Clone();
                    _onChanged();
                });
            }
            GUILayout.EndHorizontal();

            if (folds)
                GUILayout.Label("These settings fold the picture over itself near the edges, so it can't be defished. Lower the horizontal or " +
                                "vertical amount, or the steepest part of the curve.", ui.Error);
            return changed;
        }

        private static string Pct(double v)
        {
            return Math.Round(v * 100).ToString("0", CultureInfo.InvariantCulture) + "%";
        }

        /// <summary>A labelled slider in whole percent, snapping to <paramref name="snap"/> within 2%.</summary>
        private static double Slider(UiStyles ui, GUIContent label, string valueText, double value, double min, double max, double snap, float valueWidth)
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label(label, ui.Label);
            GUILayout.FlexibleSpace();
            GUILayout.Label(valueText, new GUIStyle(ui.Value) { alignment = TextAnchor.MiddleRight }, GUILayout.Width(ui.S(valueWidth)));
            GUILayout.EndHorizontal();
            float p = GUILayout.HorizontalSlider((float)(value * 100), (float)(min * 100), (float)(max * 100), ui.Slider, ui.Thumb, GUILayout.ExpandWidth(true));
            if (Math.Abs(p - value * 100) < 1e-3) return value;        // untouched: keep the value exactly
            p = Mathf.Round(p);
            if (Math.Abs(p - snap * 100) < 2) p = (float)(snap * 100);
            return p / 100.0;
        }

        // ============================================================================ curve editor

        /// <summary>
        /// The shape curves: factor (0-200%) against distance from the centre (centre to corner). The curve being edited is
        /// drawn on top with its four points; the other (when not linked) dimmed behind it.
        /// </summary>
        private bool DrawCurves(UiStyles ui, DefishAdjustment a, bool linked, bool editY, PixelImage image)
        {
            Rect r = GUILayoutUtility.GetRect(ui.S(10), ui.S(124), GUILayout.ExpandWidth(true));
            int id = GUIUtility.GetControlID(CurveHint, FocusType.Passive, r);
            float pad = ui.S(8);
            var plot = new Rect(r.x + pad, r.y + pad, r.width - 2 * pad, r.height - pad - ui.S(22));   // room below for the axis labels
            var e = Event.current;
            var edit = linked || !editY ? a.ShapeX : a.ShapeY;
            bool changed = false;

            switch (e.GetTypeForControl(id))
            {
                case EventType.MouseDown:
                    if (e.button == 0 && GUI.enabled && r.Contains(e.mousePosition))
                    {
                        int i = NearestPoint(plot, e.mousePosition.x);
                        float hy = CurveY(plot, edit[i]);
                        bool onPoint = Mathf.Abs(e.mousePosition.y - hy) <= ui.S(9);
                        _dragIndex = i;
                        _dragOffset = onPoint ? hy - e.mousePosition.y : 0;
                        GUIUtility.hotControl = id;
                        if (onPoint && e.clickCount == 2) changed |= SetPoint(a, linked, editY, i, 1);
                        else if (!onPoint) changed |= SetPoint(a, linked, editY, i, ValueAt(ui, plot, e.mousePosition.y, e.alt));
                        e.Use();
                    }
                    break;
                case EventType.MouseDrag:
                    if (GUIUtility.hotControl == id && _dragIndex >= 0)
                    {
                        changed |= SetPoint(a, linked, editY, _dragIndex, ValueAt(ui, plot, e.mousePosition.y + _dragOffset, e.alt));
                        e.Use();
                    }
                    break;
                case EventType.MouseUp:
                    if (GUIUtility.hotControl == id)
                    {
                        GUIUtility.hotControl = 0;
                        _dragIndex = -1;
                        e.Use();
                    }
                    break;
                case EventType.Repaint:
                    PaintCurves(ui, r, plot, a, linked, editY, image, GUIUtility.hotControl == id ? _dragIndex : -1);
                    break;
            }
            return changed;
        }

        private static float PointX(Rect plot, int i)
        {
            return plot.x + plot.width * i / (DefishAdjustment.ShapePoints - 1);
        }

        private static float CurveY(Rect plot, double v)
        {
            return plot.yMax - (float)(v / DefishAdjustment.MaxShape) * plot.height;
        }

        private static int NearestPoint(Rect plot, float x)
        {
            double f = (x - plot.x) / plot.width * (DefishAdjustment.ShapePoints - 1);
            return Mathf.Clamp((int)Math.Round(f), 0, DefishAdjustment.ShapePoints - 1);
        }

        /// <summary>The factor at screen height <paramref name="y"/>: whole percent, snapping to 100% nearby unless <paramref name="exact"/>.</summary>
        private static double ValueAt(UiStyles ui, Rect plot, float y, bool exact)
        {
            double v = (plot.yMax - y) / plot.height * DefishAdjustment.MaxShape;
            v = Math.Max(0, Math.Min(DefishAdjustment.MaxShape, v));
            v = Math.Round(v * 100) / 100;
            if (!exact && Math.Abs(CurveY(plot, v) - CurveY(plot, 1)) < ui.S(5)) v = 1;
            return v;
        }

        private static bool SetPoint(DefishAdjustment a, bool linked, bool editY, int i, double v)
        {
            var edit = linked || !editY ? a.ShapeX : a.ShapeY;
            bool changed = Math.Abs(edit[i] - v) > 1e-9 || (linked && Math.Abs(a.ShapeY[i] - v) > 1e-9);
            edit[i] = v;
            if (linked) a.ShapeY[i] = v;
            return changed;
        }

        private void PaintCurves(UiStyles ui, Rect r, Rect plot, DefishAdjustment a, bool linked, bool editY, PixelImage image, int dragging)
        {
            float px1 = Mathf.Max(1, Mathf.Round(ui.S(1)));
            ui.Fill(r, UiStyles.Hex(0x18191C));

            // levels: 0, 50, 100 (brighter), 150 and 200%
            for (int k = 0; k <= 4; k++)
            {
                float y = Mathf.Round(CurveY(plot, k * 0.5));
                ui.Fill(new Rect(plot.x, y, plot.width, px1), k == 2 ? UiStyles.Hex(0x4A4E56) : UiStyles.Hex(0x2A2C31));
            }

            // where the frame's edges are: the nearer edge midpoints, and the farther ones (centre to corner = 0 to 1)
            double half = Math.Sqrt(0.25 + 0.25 * _aspect * _aspect);
            double nearRho = 0.5 * _aspect / half, farRho = 0.5 / half;
            bool portrait = image != null && image.Height > image.Width;
            var tiny = new GUIStyle(ui.MutedSmall) { wordWrap = false, alignment = TextAnchor.UpperRight, fontSize = (int)ui.S(10) };
            var tinyLeft = new GUIStyle(tiny) { alignment = TextAnchor.UpperLeft };
            Marker(ui, plot, nearRho, portrait ? "sides" : "top", tiny, px1);
            Marker(ui, plot, farRho, portrait ? "top" : "sides", tiny, px1);
            GUI.Label(new Rect(plot.x + plot.width / 6 - ui.S(20), CurveY(plot, 1) - ui.S(14), ui.S(40), ui.S(14)), "100%",
                new GUIStyle(tinyLeft) { alignment = TextAnchor.UpperCenter });

            // axis labels
            float axisY = plot.yMax + ui.S(5);
            GUI.Label(new Rect(plot.x - ui.S(2), axisY, plot.width * 0.5f, ui.S(14)), "centre", tinyLeft);
            GUI.Label(new Rect(plot.x + plot.width * 0.5f, axisY, plot.width * 0.5f + ui.S(2), ui.S(14)), "corners", tiny);

            // curves: the other one dimmed, then the one being edited
            if (!linked)
            {
                var other = editY ? a.ShapeX : a.ShapeY;
                var c = editY ? ColorX : ColorY;
                PaintCurve(ui, plot, other, new Color(c.r, c.g, c.b, 0.35f), px1);
            }
            var edit = linked || !editY ? a.ShapeX : a.ShapeY;
            var col = linked || !editY ? ColorX : ColorY;
            PaintCurve(ui, plot, edit, col, Mathf.Max(px1, Mathf.Round(ui.S(2))));

            // the points
            float size = Mathf.Round(ui.S(9));
            for (int i = 0; i < DefishAdjustment.ShapePoints; i++)
            {
                float x = PointX(plot, i), y = CurveY(plot, edit[i]);
                float sz = i == dragging ? Mathf.Round(ui.S(12)) : size;
                var box = new Rect(Mathf.Round(x - sz / 2), Mathf.Round(y - sz / 2), sz, sz);
                ui.Fill(new Rect(box.x - px1, box.y - px1, box.width + 2 * px1, box.height + 2 * px1), UiStyles.Hex(0x111214));
                ui.Fill(box, i == dragging ? Color.white : col);
                if (i == dragging)
                {
                    bool left = i == DefishAdjustment.ShapePoints - 1;
                    var lr = left ? new Rect(x - ui.S(54), y - ui.S(22), ui.S(46), ui.S(14)) : new Rect(x + ui.S(8), y - ui.S(22), ui.S(46), ui.S(14));
                    GUI.Label(lr, Pct(edit[i]), new GUIStyle(ui.Value) { fontSize = (int)ui.S(11), alignment = left ? TextAnchor.UpperRight : TextAnchor.UpperLeft });
                }
            }
        }

        private static void Marker(UiStyles ui, Rect plot, double rho, string label, GUIStyle style, float px1)
        {
            if (rho <= 0 || rho >= 1) return;
            float x = Mathf.Round(plot.x + (float)rho * plot.width);
            // dashed
            float dash = Mathf.Max(2, Mathf.Round(ui.S(3)));
            for (float y = plot.y; y < plot.yMax; y += 2 * dash)
                ui.Fill(new Rect(x, y, px1, Mathf.Min(dash, plot.yMax - y)), UiStyles.Hex(0x3A3D44));
            GUI.Label(new Rect(x - ui.S(60) - ui.S(3), plot.y, ui.S(60), ui.S(14)), label, style);
        }

        /// <summary>Draws the curve as a run of thin columns, one per pixel across the plot.</summary>
        private static void PaintCurve(UiStyles ui, Rect plot, double[] shape, Color c, float thickness)
        {
            int n = Mathf.Max(2, Mathf.RoundToInt(plot.width));
            float prev = CurveY(plot, DefishAdjustment.ShapeAt(shape, 0));
            for (int i = 0; i <= n; i++)
            {
                float y = CurveY(plot, DefishAdjustment.ShapeAt(shape, (double)i / n));
                float top = Mathf.Min(prev, y) - thickness * 0.5f, bottom = Mathf.Max(prev, y) + thickness * 0.5f;
                ui.Fill(new Rect(plot.x + plot.width * i / n - 0.5f, top, 1, bottom - top), c);
                prev = y;
            }
        }

        // ============================================================================ helpers

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

        public void Dispose()
        {
            _preview.Dispose();
        }
    }
}
