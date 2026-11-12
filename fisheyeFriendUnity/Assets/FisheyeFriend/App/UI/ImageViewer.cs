using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using FisheyeFriend.Core.Imaging;
using UnityEngine;

namespace FisheyeFriend.App.UI
{
    /// <summary>
    /// GPU copy of one image for the viewer: point-filtered tiles (so images larger than the GPU's
    /// texture limit still display), plus whole-number reductions built on demand for zooming out.
    /// </summary>
    internal sealed class ImageLayer : IDisposable
    {
        private const int TileSize = 2048;

        internal sealed class Tile
        {
            public int X, Y, W, H;
            public Texture2D Tex;
        }

        internal sealed class TileSet
        {
            public int Width, Height;
            public readonly List<Tile> Tiles = new List<Tile>();
        }

        private sealed class PreparedTile
        {
            public int X, Y, W, H;
            public byte[] Bytes;
        }

        public PixelImage Image { get; private set; }
        public TileSet Full { get; private set; }

        private readonly Dictionary<int, TileSet> _reduced = new Dictionary<int, TileSet>();
        private readonly HashSet<int> _building = new HashSet<int>();
        private bool _disposed;

        private ImageLayer() { }

        /// <summary>Prepares tiles on worker threads, then uploads them (call from the main thread).</summary>
        public static async Task<ImageLayer> CreateAsync(PixelImage image)
        {
            var layer = new ImageLayer { Image = image };
            var prepared = await Task.Run(() => PrepareFull(image));
            layer.Full = Upload(prepared, image.Width, image.Height);
            return layer;
        }

        /// <summary>Returns the 1/n reduction if ready; otherwise starts building it and returns null.</summary>
        public TileSet GetReduced(int n)
        {
            if (n <= 1) return Full;
            TileSet set;
            if (_reduced.TryGetValue(n, out set)) return set;
            if (_building.Add(n)) BuildReducedAsync(n);
            return null;
        }

        private async void BuildReducedAsync(int n)
        {
            try
            {
                var img = Image;
                int w = 0, h = 0;
                var prepared = await Task.Run(() =>
                {
                    byte[] rgba = DisplayReducer.ReduceToRgba32(img, n, out w, out h);
                    return SplitTopDown(rgba, w, h);
                });
                if (_disposed) return;
                _reduced[n] = Upload(prepared, w, h);
            }
            catch (Exception e)
            {
                Debug.LogWarning("Could not build zoomed-out preview: " + e.Message);
            }
            finally
            {
                _building.Remove(n);
            }
        }

        private static List<PreparedTile> PrepareFull(PixelImage img)
        {
            var list = new List<PreparedTile>();
            for (int y = 0; y < img.Height; y += TileSize)
                for (int x = 0; x < img.Width; x += TileSize)
                {
                    int w = Math.Min(TileSize, img.Width - x);
                    int h = Math.Min(TileSize, img.Height - y);
                    var bytes = new byte[w * h * 4];
                    img.CopyToRgba32(x, y, w, h, bytes, true); // Unity textures are bottom-up
                    list.Add(new PreparedTile { X = x, Y = y, W = w, H = h, Bytes = bytes });
                }
            return list;
        }

        private static List<PreparedTile> SplitTopDown(byte[] rgba, int width, int height)
        {
            var list = new List<PreparedTile>();
            for (int y = 0; y < height; y += TileSize)
                for (int x = 0; x < width; x += TileSize)
                {
                    int w = Math.Min(TileSize, width - x);
                    int h = Math.Min(TileSize, height - y);
                    var bytes = new byte[w * h * 4];
                    for (int r = 0; r < h; r++)
                        Buffer.BlockCopy(rgba, ((y + r) * width + x) * 4, bytes, (h - 1 - r) * w * 4, w * 4);
                    list.Add(new PreparedTile { X = x, Y = y, W = w, H = h, Bytes = bytes });
                }
            return list;
        }

        private static TileSet Upload(List<PreparedTile> prepared, int width, int height)
        {
            var set = new TileSet { Width = width, Height = height };
            foreach (var p in prepared)
            {
                var tex = new Texture2D(p.W, p.H, TextureFormat.RGBA32, false, false)
                {
                    filterMode = FilterMode.Point,
                    wrapMode = TextureWrapMode.Clamp,
                    hideFlags = HideFlags.HideAndDontSave,
                    anisoLevel = 0,
                };
                tex.LoadRawTextureData(p.Bytes);
                tex.Apply(false, true); // upload and free the CPU copy
                set.Tiles.Add(new Tile { X = p.X, Y = p.Y, W = p.W, H = p.H, Tex = tex });
            }
            return set;
        }

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;
            Destroy(Full);
            foreach (var s in _reduced.Values) Destroy(s);
            _reduced.Clear();
            Full = null;
            Image = null;
        }

        private static void Destroy(TileSet set)
        {
            if (set == null) return;
            foreach (var t in set.Tiles)
                if (t.Tex != null) UnityEngine.Object.Destroy(t.Tex);
            set.Tiles.Clear();
        }
    }

    /// <summary>
    /// Pixel-perfect zoom and pan. Magnification is always a whole number (each image pixel becomes
    /// an exact N x N block of screen pixels) and zooming out uses exact 1/N block averages.
    /// The image is always drawn at whole screen-pixel offsets.
    /// </summary>
    internal sealed class ImageViewer
    {
        /// <summary>Zoom steps: negative = 1/n reduction, positive = n x magnification.</summary>
        private static readonly int[] Steps =
        {
            -32, -24, -16, -12, -10, -8, -6, -5, -4, -3, -2, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64, 96, 128, 192, 256,
        };

        private ImageLayer _layer;
        private int _zoom = 1;
        private double _cx, _cy;          // view centre, in image pixels
        private int _shownW, _shownH;     // size of the image _cx/_cy refer to (kept even if that layer is released)
        private bool _fitPending = true;
        private int _maxZoom = 64;
        private float _scroll;
        private Vector2 _dragLast;
        private int _controlId;
        private Vector2 _viewSizePx;

        public Vector2Int? Hover { get; private set; }
        public ImageLayer Layer => _layer;
        public int ZoomValue => _zoom;

        /// <summary>A saved zoom level and position (see <see cref="Capture"/> / <see cref="Restore"/>).</summary>
        internal struct ViewState
        {
            public int Zoom;
            public double Cx, Cy;
            public bool Valid;
        }

        public ViewState Capture()
        {
            return new ViewState { Zoom = _zoom, Cx = _cx, Cy = _cy, Valid = _layer != null };
        }

        /// <summary>Restores a saved view; call after <see cref="Show"/> of the layer it was captured on.</summary>
        public void Restore(ViewState s)
        {
            if (!s.Valid || _layer == null) return;
            _zoom = s.Zoom;
            _cx = s.Cx;
            _cy = s.Cy;
            _fitPending = false;
            ClampCentre();
        }

        /// <summary>View centre in image pixels (pixel edges at whole numbers).</summary>
        public void GetCentre(out double x, out double y)
        {
            x = _cx;
            y = _cy;
        }

        public void SetCentre(double x, double y)
        {
            _cx = x;
            _cy = y;
            ClampCentre();
        }

        /// <summary>Screen pixels per image pixel at the current zoom.</summary>
        public double ScreenPixelsPerImagePixel => ScaleFactor;

        /// <summary>
        /// GUI position of an image point given in pixel-centre coordinates (pixel (0,0) is centred on (0,0)), as drawn by
        /// the last <see cref="Draw"/> into <paramref name="view"/>.
        /// </summary>
        public bool ImageToGui(Rect view, double x, double y, out Vector2 gui)
        {
            gui = default(Vector2);
            if (_layer == null) return false;
            double ox, oy;
            Origin(out ox, out oy);
            double s = ScaleFactor;
            float ppp = UiStyles.PixelsPerPoint;
            gui = new Vector2((float)(view.x + (ox + (x + 0.5) * s) / ppp), (float)(view.y + (oy + (y + 0.5) * s) / ppp));
            return true;
        }

        /// <summary>The image point (pixel-centre coordinates) under a GUI position; false if it is outside the image.</summary>
        public bool GuiToImage(Rect view, Vector2 gui, out double x, out double y)
        {
            x = y = 0;
            if (_layer == null) return false;
            double ox, oy;
            Origin(out ox, out oy);
            double s = ScaleFactor;
            float ppp = UiStyles.PixelsPerPoint;
            x = ((gui.x - view.x) * ppp - ox) / s - 0.5;
            y = ((gui.y - view.y) * ppp - oy) / s - 0.5;
            return x >= -0.5 && y >= -0.5 && x <= _layer.Image.Width - 0.5 && y <= _layer.Image.Height - 0.5;
        }

        public void SetMaxZoom(int max)
        {
            _maxZoom = Math.Max(1, max);
        }

        /// <summary>
        /// Shows a layer. With keepView the zoom level and the position relative to the image (e.g. "40% across,
        /// 25% down") are preserved, so the same spot stays under the view when switching between the original
        /// and a result, or between results of different settings, even though their pixel sizes differ.
        /// </summary>
        public void Show(ImageLayer layer, bool keepView)
        {
            if (keepView && layer != null && layer.Image != null && _shownW > 0 && _shownH > 0)
            {
                double rx = _cx / _shownW, ry = _cy / _shownH;
                _layer = layer;
                _shownW = layer.Image.Width;
                _shownH = layer.Image.Height;
                _cx = rx * _shownW;
                _cy = ry * _shownH;
                return;
            }
            _layer = layer;
            if (layer != null && layer.Image != null)
            {
                _shownW = layer.Image.Width;
                _shownH = layer.Image.Height;
                _cx = _shownW * 0.5;
                _cy = _shownH * 0.5;
                _fitPending = true;
            }
            else
            {
                _shownW = _shownH = 0;
                _fitPending = true;
            }
        }

        public string ZoomLabel
        {
            get
            {
                if (_zoom >= 1) return (_zoom * 100) + "%";
                int n = -_zoom;
                return "1/" + n + " (" + (100.0 / n).ToString("0.#") + "%)";
            }
        }

        private int Mag => _zoom >= 1 ? _zoom : 1;
        private int Red => _zoom >= 1 ? 1 : -_zoom;
        private double ScaleFactor => (double)Mag / Red;

        public void ZoomIn() { StepZoom(+1, null); }
        public void ZoomOut() { StepZoom(-1, null); }

        public void ActualSize()
        {
            SetZoom(1, null);
        }

        public void Fit()
        {
            _fitPending = true;
        }

        private void FitNow()
        {
            if (_layer == null || _viewSizePx.x < 4 || _viewSizePx.y < 4) return;
            double availW = Math.Max(1, _viewSizePx.x - 16), availH = Math.Max(1, _viewSizePx.y - 16);
            int best = Steps[0];
            foreach (int z in Steps)
            {
                if (z > _maxZoom) break;
                int mag = z >= 1 ? z : 1, red = z >= 1 ? 1 : -z;
                double w = Math.Ceiling(_layer.Image.Width / (double)red) * mag;
                double h = Math.Ceiling(_layer.Image.Height / (double)red) * mag;
                if (w <= availW && h <= availH) best = z;
            }
            _zoom = best;
            _cx = _layer.Image.Width * 0.5;
            _cy = _layer.Image.Height * 0.5;
            _fitPending = false;
        }

        private void StepZoom(int dir, Vector2? anchorPx)
        {
            int idx = Array.IndexOf(Steps, _zoom);
            if (idx < 0) idx = Array.IndexOf(Steps, 1);
            idx = Mathf.Clamp(idx + dir, 0, Steps.Length - 1);
            if (Steps[idx] > _maxZoom) return;
            SetZoom(Steps[idx], anchorPx);
        }

        private void SetZoom(int zoom, Vector2? anchorPx)
        {
            if (_layer == null) { _zoom = zoom; return; }
            double oldScale = ScaleFactor;
            Vector2 a = anchorPx ?? _viewSizePx * 0.5f;
            double ax = a.x - _viewSizePx.x * 0.5, ay = a.y - _viewSizePx.y * 0.5;
            double ix = _cx + ax / oldScale, iy = _cy + ay / oldScale; // image point under the anchor
            _zoom = zoom;
            double s = ScaleFactor;
            _cx = ix - ax / s;
            _cy = iy - ay / s;
            ClampCentre();
            _fitPending = false;
        }

        public void PanPixels(float dx, float dy)
        {
            double s = ScaleFactor;
            _cx -= dx / s;
            _cy -= dy / s;
            ClampCentre();
        }

        private void ClampCentre()
        {
            if (_layer == null) return;
            _cx = Math.Max(0, Math.Min(_layer.Image.Width, _cx));
            _cy = Math.Max(0, Math.Min(_layer.Image.Height, _cy));
        }

        /// <summary>Handles mouse input for the view area. Call once per OnGUI event before drawing.</summary>
        public void HandleInput(Rect view)
        {
            Event e = Event.current;
            _controlId = GUIUtility.GetControlID(FocusType.Passive);
            float ppp = UiStyles.PixelsPerPoint;
            Vector2 mp = e.mousePosition;
            bool inside = view.Contains(mp);
            Vector2 localPx = (mp - view.position) * ppp;

            switch (e.GetTypeForControl(_controlId))
            {
                case EventType.MouseDown:
                    if (inside && _layer != null && (e.button == 0 || e.button == 2))
                    {
                        if (e.clickCount == 2 && e.button == 0)
                        {
                            if (_zoom == 1) Fit();
                            else SetZoom(1, localPx);
                        }
                        GUIUtility.hotControl = _controlId;
                        GUIUtility.keyboardControl = 0;
                        _dragLast = mp;
                        e.Use();
                    }
                    break;
                case EventType.MouseDrag:
                    if (GUIUtility.hotControl == _controlId)
                    {
                        Vector2 d = (mp - _dragLast) * ppp;
                        _dragLast = mp;
                        PanPixels(d.x, d.y);
                        e.Use();
                    }
                    break;
                case EventType.MouseUp:
                    if (GUIUtility.hotControl == _controlId)
                    {
                        GUIUtility.hotControl = 0;
                        e.Use();
                    }
                    break;
                case EventType.ScrollWheel:
                    if (inside && _layer != null)
                    {
                        _scroll += e.delta.y;
                        if (Mathf.Abs(_scroll) >= 1.5f)
                        {
                            StepZoom(_scroll < 0 ? +1 : -1, localPx);
                            _scroll = 0;
                        }
                        e.Use();
                    }
                    break;
            }

            if (e.type == EventType.Repaint || e.type == EventType.MouseMove || e.type == EventType.MouseDrag)
                UpdateHover(inside, localPx);
        }

        private void UpdateHover(bool inside, Vector2 localPx)
        {
            if (!inside || _layer == null)
            {
                Hover = null;
                return;
            }
            double ox, oy;
            Origin(out ox, out oy);
            double s = ScaleFactor;
            int x = (int)Math.Floor((localPx.x - ox) / s);
            int y = (int)Math.Floor((localPx.y - oy) / s);
            if (x < 0 || y < 0 || x >= _layer.Image.Width || y >= _layer.Image.Height) Hover = null;
            else Hover = new Vector2Int(x, y);
        }

        /// <summary>Top-left of the image in view pixels, snapped to whole pixels.</summary>
        private void Origin(out double ox, out double oy)
        {
            // Position in the displayed (possibly reduced) image grid, then whole-pixel snap.
            ox = Math.Floor(_viewSizePx.x * 0.5 - (_cx / Red) * Mag);
            oy = Math.Floor(_viewSizePx.y * 0.5 - (_cy / Red) * Mag);
        }

        public void Draw(Rect view, UiStyles ui)
        {
            float ppp = UiStyles.PixelsPerPoint;
            _viewSizePx = new Vector2(Mathf.Round(view.width * ppp), Mathf.Round(view.height * ppp));
            ui.Fill(view, UiStyles.ViewerBackground);
            if (_layer == null || _layer.Full == null) return;
            if (_fitPending) FitNow();
            if (Event.current.type != EventType.Repaint) return;

            ImageLayer.TileSet set = Red == 1 ? _layer.Full : _layer.GetReduced(Red);
            bool fallback = set == null; // reduction still being built: draw full tiles scaled for now
            if (fallback) set = _layer.Full;

            double ox, oy;
            Origin(out ox, out oy);
            double tileScale = fallback ? Mag / (double)Red : Mag;
            bool alpha = _layer.Image.HasAlpha;

            GUI.BeginClip(view);
            foreach (var t in set.Tiles)
            {
                double x = ox + t.X * tileScale, y = oy + t.Y * tileScale;
                double w = t.W * tileScale, h = t.H * tileScale;
                if (x >= _viewSizePx.x || y >= _viewSizePx.y || x + w <= 0 || y + h <= 0) continue;
                var r = new Rect((float)(x / ppp), (float)(y / ppp), (float)(w / ppp), (float)(h / ppp));
                if (alpha)
                {
                    float cell = 8f / ppp;
                    GUI.DrawTextureWithTexCoords(r, ui.Checker, new Rect(0, 0, r.width / (2 * cell), r.height / (2 * cell)));
                }
                GUI.DrawTexture(r, t.Tex, ScaleMode.StretchToFill, true);
            }
            GUI.EndClip();
        }
    }
}
