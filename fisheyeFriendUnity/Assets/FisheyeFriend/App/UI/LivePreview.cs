using System;
using System.Globalization;
using System.Threading.Tasks;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using UnityEngine;

namespace FisheyeFriend.App.UI
{
    /// <summary>
    /// The small defished picture in the settings panel, kept up to date as the defish is adjusted. A reduced copy of the
    /// open photo is made once; each change then makes a coarse geometry and warps that copy on a worker (some tens of
    /// milliseconds). Only one picture is made at a time: changes that arrive meanwhile are folded into the next one.
    ///
    /// Everything the panel reads (<see cref="Texture"/>, <see cref="Geometry"/>...) changes only in <see cref="Tick"/>,
    /// called from Update, so it stays the same between the layout and repaint passes of one GUI event.
    /// </summary>
    internal sealed class LivePreview : IDisposable
    {
        /// <summary>Longest side of the reduced copy of the photo: enough for a sharp picture up to about 900 pixels wide.</summary>
        private const int ReducedSide = 1400;

        private sealed class Reduced
        {
            public byte[] Pixels;
            public int Width, Height;
        }

        private sealed class Result
        {
            public DefishGeometry Geometry;
            public string GeometryKey;
            public byte[] Pixels;       // null when there was no photo or size to draw
            public int Width, Height;
        }

        // the wish (set by Want)
        private PixelImage _wantImage;
        private double _wantAspect = 2.0 / 3.0;
        private DefishAdjustment _wantAdjust = new DefishAdjustment();
        private int _wantW, _wantH;

        // the photo's reduced copy
        private PixelImage _image;
        private int _imageSerial;
        private Task<Reduced> _reducing;
        private Reduced _reduced;

        // the picture being made, and the last one made
        private Task<Result> _job;
        private string _jobKey, _doneKey;
        private int _jobSerial;

        public Texture2D Texture { get; private set; }
        /// <summary>True when <see cref="Texture"/> shows the photo now open (perhaps with slightly older settings).</summary>
        public bool HasPicture { get; private set; }
        /// <summary>True while the picture is behind the latest settings.</summary>
        public bool Behind { get; private set; }
        /// <summary>The coarse geometry of the latest finished picture (made even without a photo, for the fold check).</summary>
        public DefishGeometry Geometry { get; private set; }
        /// <summary>Which frame shape and adjustment <see cref="Geometry"/> is for (see <see cref="KeyFor"/>).</summary>
        public string GeometryKey { get; private set; }
        public string Error { get; private set; }

        public static string KeyFor(double aspect, DefishAdjustment adj)
        {
            var a = (adj ?? new DefishAdjustment()).Clone();
            a.Validate();
            return aspect.ToString("0.00000", CultureInfo.InvariantCulture) + "|" + a.Key;
        }

        /// <summary>What the picture should show: the photo (or null), its lens frame shape, the adjustment and the size to draw.</summary>
        public void Want(PixelImage image, double aspect, DefishAdjustment adjustment, int width, int height)
        {
            _wantImage = image;
            _wantAspect = aspect;
            _wantAdjust = (adjustment ?? new DefishAdjustment()).Clone();
            _wantAdjust.Validate();
            _wantW = Math.Max(0, Math.Min(2048, width));
            _wantH = Math.Max(0, Math.Min(2048, height));
        }

        /// <summary>Takes in finished work and starts the next; call once per frame from Update.</summary>
        public void Tick(Projection projection)
        {
            // a different photo: start its reduced copy, and stop showing the last one's picture
            if (!ReferenceEquals(_wantImage, _image))
            {
                _image = _wantImage;
                _imageSerial++;
                _reduced = null;
                HasPicture = false;
                var img = _image;
                _reducing = img == null ? null : Task.Run(() => Reduce(img));
            }
            if (_reducing != null && _reducing.IsCompleted)
            {
                if (_reducing.Status == TaskStatus.RanToCompletion) _reduced = _reducing.Result;
                else Error = "Could not prepare the preview: " + Message(_reducing.Exception);
                _reducing = null;
            }

            if (_job != null && _job.IsCompleted)
            {
                if (_job.Status == TaskStatus.RanToCompletion)
                {
                    var r = _job.Result;
                    Geometry = r.Geometry;
                    GeometryKey = r.GeometryKey;
                    if (r.Pixels != null && _jobSerial == _imageSerial)
                    {
                        Upload(r);
                        HasPicture = true;
                    }
                    Error = null;
                }
                else
                {
                    Error = "Could not make the preview: " + Message(_job.Exception);
                }
                _doneKey = _jobKey;    // a failed wish is not retried until it changes
                _job = null;
            }

            if (!projection.Ready || projection.Failed)
            {
                Behind = false;
                return;
            }
            string gkey = KeyFor(_wantAspect, _wantAdjust);
            bool draw = _reduced != null && _wantW > 0 && _wantH > 0;
            string key = gkey + "|" + _imageSerial + "|" + (draw ? _wantW + "x" + _wantH : "-");
            Behind = key != _doneKey;
            if (_job != null || !Behind) return;

            // start the next picture
            var reduced = draw ? _reduced : null;
            int turn = reduced != null ? Defisher.SensorTurn(_image) : 0;
            int w = _wantW, h = _wantH;
            double aspect = _wantAspect;
            var adj = _wantAdjust.Clone();
            var reuse = GeometryKey == gkey ? Geometry : null;
            _jobKey = key;
            _jobSerial = _imageSerial;
            _job = Task.Run(() =>
            {
                var g = reuse ?? projection.CoarseGeometry(aspect, adj);
                var r = new Result { Geometry = g, GeometryKey = gkey };
                if (reduced != null && g != null)
                {
                    r.Pixels = DefishPreview.Render(reduced.Pixels, reduced.Width, reduced.Height, turn, g, w, h, true);
                    r.Width = w;
                    r.Height = h;
                }
                return r;
            });
        }

        private static Reduced Reduce(PixelImage img)
        {
            int side = Math.Max(img.Width, img.Height);
            int n = Math.Max(1, side / ReducedSide);
            int w, h;
            var px = DisplayReducer.ReduceToRgba32(img, n, out w, out h);
            return new Reduced { Pixels = px, Width = w, Height = h };
        }

        private void Upload(Result r)
        {
            if (Texture == null || Texture.width != r.Width || Texture.height != r.Height)
            {
                if (Texture != null) UnityEngine.Object.Destroy(Texture);
                Texture = new Texture2D(r.Width, r.Height, TextureFormat.RGBA32, false, false)
                {
                    filterMode = FilterMode.Bilinear,
                    wrapMode = TextureWrapMode.Clamp,
                    hideFlags = HideFlags.HideAndDontSave,
                };
            }
            Texture.LoadRawTextureData(r.Pixels);
            Texture.Apply(false, false);
        }

        private static string Message(AggregateException e)
        {
            var b = e != null ? e.GetBaseException() : null;
            return b != null ? b.Message : "cancelled";
        }

        public void Dispose()
        {
            if (Texture != null) UnityEngine.Object.Destroy(Texture);
            Texture = null;
            HasPicture = false;
        }
    }
}
