using System;
using System.Threading;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>How large the defished image is.</summary>
    public enum DefishOutputSize
    {
        /// <summary>The photo's own size (as Fisheye Friend always made it). The Hemi projection shrinks the centre of
        /// the frame by about 1.19x here, so some of the centre's detail is averaged away.</summary>
        SameAsPhoto = 0,
        /// <summary>Enlarged so the centre of the frame keeps every pixel of the photo (about 1.19x each way).</summary>
        FullDetail = 1,
    }

    /// <summary>The few output pixels the projection places just outside the photo.</summary>
    public enum DefishEdges
    {
        /// <summary>Keep the full frame; those pixels repeat the photo's edge (as Fisheye Friend did). On a 42 MP frame they reach about 5 pixels past it, along the top and bottom.</summary>
        Extend = 0,
        /// <summary>Crop slightly (about 0.3% each way) to the largest centred rectangle of the frame's shape that the photo covers completely.</summary>
        Crop = 1,
    }

    public sealed class DefishSettings
    {
        public DefishOutputSize OutputSize = DefishOutputSize.SameAsPhoto;
        public DefishEdges Edges = DefishEdges.Extend;
        /// <summary>Warp filter; see <see cref="WarpKernel"/>. Separable Lanczos 4 measured best (Lanczos 3: 0.03 dB less, about 1.5x faster).</summary>
        public WarpKernel Kernel = WarpKernel.TensorLanczos;
        public double Radius = 4;
        /// <summary>Sigmoidal measured best (fewest halos where the projection enlarges).</summary>
        public LightSpace Light = LightSpace.Sigmoidal;
        /// <summary>
        /// Widen the filter where the projection shrinks the photo (the centre, by about 1.19x), so detail averages instead
        /// of aliasing. Measured 0.1-0.2 dB worse at that mild shrink, so off by default.
        /// </summary>
        public bool Antialias;
        /// <summary>
        /// Where the projection enlarges the photo (towards the edges, up to about 1.8x), partly undo the photo's own pixel
        /// blur so each result pixel shows its own, smaller footprint (see <see cref="WarpOptions.Aperture"/>). +0.1-0.3 dB.
        /// </summary>
        public bool ApertureCorrection = true;
        /// <summary>The user's adjustment of the defish model (strength, horizontal and vertical amounts, shape across the frame).</summary>
        public DefishAdjustment Adjust = new DefishAdjustment();

        public DefishSettings Clone()
        {
            var c = (DefishSettings)MemberwiseClone();
            c.Adjust = (Adjust ?? new DefishAdjustment()).Clone();
            return c;
        }

        public void Validate()
        {
            if (!Enum.IsDefined(typeof(DefishOutputSize), OutputSize)) OutputSize = DefishOutputSize.SameAsPhoto;
            if (!Enum.IsDefined(typeof(DefishEdges), Edges)) Edges = DefishEdges.Extend;
            if (!Enum.IsDefined(typeof(WarpKernel), Kernel)) Kernel = WarpKernel.TensorLanczos;
            if (!Enum.IsDefined(typeof(LightSpace), Light)) Light = LightSpace.Sigmoidal;
            if (double.IsNaN(Radius) || Radius < 1) Radius = 3;
            if (Radius > 4) Radius = 4;
            if (Adjust == null) Adjust = new DefishAdjustment();
            Adjust.Validate();
        }

        public bool Same(DefishSettings o)
        {
            return o != null && OutputSize == o.OutputSize && Edges == o.Edges && Kernel == o.Kernel && Math.Abs(Radius - o.Radius) < 1e-9 &&
                   Light == o.Light && Antialias == o.Antialias && ApertureCorrection == o.ApertureCorrection &&
                   (Adjust ?? new DefishAdjustment()).Same(o.Adjust ?? new DefishAdjustment());
        }

        public override string ToString()
        {
            return OutputSize + ", " + Edges + ", " + Kernel + (Kernel == WarpKernel.Bicubic ? "" : " " + Radius) + ", " + Light +
                   (Antialias ? ", AA" : "") + (ApertureCorrection ? ", aperture" : "") +
                   (Adjust != null && !Adjust.IsDefault ? ", adjust " + Adjust.Key : "");
        }
    }

    /// <summary>Normalised part of the defished frame an output shows.</summary>
    public struct DefishCanvas
    {
        public double U0, V0, U1, V1;

        public DefishCanvas(double u0, double v0, double u1, double v1)
        {
            U0 = u0;
            V0 = v0;
            U1 = u1;
            V1 = v1;
        }

        public double Width => U1 - U0;
        public double Height => V1 - V0;
        public static DefishCanvas Full => new DefishCanvas(0, 0, 1, 1);
    }

    /// <summary>
    /// Defishes photos with a <see cref="DefishGeometry"/>. The projection belongs to the lens on the sensor, so it is applied
    /// in the sensor's landscape frame: a portrait photo is turned the way the camera's orientation tag says (or a quarter
    /// turn clockwise), defished, and turned back.
    /// </summary>
    public static class Defisher
    {
        /// <summary>Height / width of the photo as the lens saw it (in the sensor's landscape frame): the shape a geometry is made for.</summary>
        public static double SensorAspect(PixelImage upright)
        {
            bool swap = Orientation.SwapsAxes(SensorTurn(upright));
            double w = swap ? upright.Height : upright.Width, h = swap ? upright.Width : upright.Height;
            return h / w;
        }

        /// <summary>The turn from the upright photo to the landscape frame the projection was measured in.</summary>
        public static int SensorTurn(PixelImage upright)
        {
            if (upright.Width >= upright.Height) return Orientation.Normal;
            int o = upright.Metadata != null ? upright.Metadata.SourceOrientation : 0;
            return Orientation.SwapsAxes(o) ? Orientation.Inverse(o) : Orientation.Rotate90Clockwise;
        }

        /// <summary>Source pixels per output pixel at the centre of the frame, for an output the size of the photo (about 1.19).</summary>
        public static double CentreScale(DefishGeometry g, int width, int height)
        {
            double s, t, su, sv, tu, tv;
            g.Map(0.5, 0.5, out s, out t, out su, out sv, out tu, out tv);
            // pixel units: d(src px)/d(out px) with out = src size -> the normalised Jacobian itself (aspect cancels on the diagonal)
            double j00 = su, j01 = sv * width / height, j10 = tu * height / width, j11 = tv;
            return Math.Sqrt(Math.Abs(j00 * j11 - j01 * j10));
        }

        /// <summary>
        /// The largest centred part of the frame, with the frame's shape, whose outline maps inside the photo (by at least
        /// half a pixel of a photo <paramref name="width"/> wide).
        /// </summary>
        public static DefishCanvas CropCanvas(DefishGeometry g, int width, int height)
        {
            double lo = 0.5, hi = 1.0;
            if (OutlineInside(g, 1.0, width, height)) return DefishCanvas.Full;
            for (int it = 0; it < 40; it++)
            {
                double m = 0.5 * (lo + hi);
                if (OutlineInside(g, m, width, height)) lo = m;
                else hi = m;
            }
            return new DefishCanvas(0.5 - lo / 2, 0.5 - lo / 2, 0.5 + lo / 2, 0.5 + lo / 2);
        }

        private static bool OutlineInside(DefishGeometry g, double c, int width, int height)
        {
            double mx = 0.5 / width, my = 0.5 / height;
            const int N = 400;
            for (int i = 0; i <= N; i++)
            {
                double t = 0.5 - c / 2 + c * i / N;
                for (int k = 0; k < 4; k++)
                {
                    double u = k < 2 ? t : (k == 2 ? 0.5 - c / 2 : 0.5 + c / 2);
                    double v = k < 2 ? (k == 0 ? 0.5 - c / 2 : 0.5 + c / 2) : t;
                    double s, tt, a, b, cc, d;
                    g.Map(u, v, out s, out tt, out a, out b, out cc, out d);
                    if (s < mx || s > 1 - mx || tt < my || tt > 1 - my) return false;
                }
            }
            return true;
        }

        /// <summary>Canvas and output size (in the sensor frame) for a landscape photo <paramref name="w"/> x <paramref name="h"/>.</summary>
        public static void Plan(DefishGeometry g, int w, int h, DefishSettings s, out DefishCanvas canvas, out int outW, out int outH)
        {
            canvas = s.Edges == DefishEdges.Crop ? CropCanvas(g, w, h) : DefishCanvas.Full;
            double k = s.OutputSize == DefishOutputSize.FullDetail ? CentreScale(g, w, h) : 1.0;
            outW = Math.Max(1, (int)Math.Round(w * canvas.Width * k));
            outH = Math.Max(1, (int)Math.Round(h * canvas.Height * k));
        }

        /// <summary>The size of the image <see cref="Process"/> returns for this photo.</summary>
        public static void OutputSize(DefishGeometry g, PixelImage upright, DefishSettings s, out int width, out int height)
        {
            bool swap = Orientation.SwapsAxes(SensorTurn(upright));
            int w = swap ? upright.Height : upright.Width, h = swap ? upright.Width : upright.Height;
            DefishCanvas c;
            int ow, oh;
            Plan(g, w, h, s, out c, out ow, out oh);
            width = swap ? oh : ow;
            height = swap ? ow : oh;
        }

        public static WarpOptions WarpOptionsFor(PixelImage src, DefishSettings s)
        {
            return new WarpOptions
            {
                Kernel = s.Kernel,
                Radius = s.Kernel == WarpKernel.Bicubic ? 2 : s.Radius,
                Light = src.HasAlpha && s.Light == LightSpace.Sigmoidal ? LightSpace.Gamma : s.Light,
                Antialias = s.Antialias,
                Aperture = s.ApertureCorrection && s.Kernel != WarpKernel.RadialLanczos ? 1 : 0,
                ApertureSharpenOnly = true,
                ExtendEdges = s.Edges == DefishEdges.Extend,
                OutputBitDepth = src.BitDepth,
            };
        }

        /// <summary>Defishes an upright photo; the result is upright, with the photo's metadata.</summary>
        public static PixelImage Process(PixelImage upright, DefishGeometry g, DefishSettings settings,
            ProgressTracker progress = null, CancellationToken ct = default)
        {
            if (upright == null) throw new ArgumentNullException(nameof(upright));
            if (g == null) throw new ArgumentNullException(nameof(g));
            if (g.Folds)
                throw new InvalidOperationException("These defish settings fold the picture over itself near the edges. Lower the horizontal or vertical amount, or the curve.");
            var s = settings.Clone();
            s.Validate();
            int turn = SensorTurn(upright);
            var sensor = Orientation.Apply(upright, turn);
            DefishCanvas c;
            int ow, oh;
            Plan(g, sensor.Width, sensor.Height, s, out c, out ow, out oh);
            var map = new DefishWarpMap(g, sensor.Width, sensor.Height, ow, oh, c.U0, c.V0, c.U1, c.V1);
            var r = WarpResampler.Warp(sensor, ow, oh, map, WarpOptionsFor(sensor, s), progress, ct);
            sensor = null;
            ct.ThrowIfCancellationRequested();
            if (turn != Orientation.Normal) r = Orientation.Apply(r, Orientation.Inverse(turn));
            r.Metadata = upright.Metadata != null ? upright.Metadata.Clone() : new ImageMetadata();
            r.SourceFormat = upright.SourceFormat;
            return r;
        }

        /// <summary>
        /// Where the photo point (<paramref name="x"/>, <paramref name="y"/>) (upright, pixel centres at whole numbers) appears
        /// in the result of <see cref="Process"/> with these settings. False if the result does not show it.
        /// </summary>
        public static bool PhotoToResult(DefishGeometry g, PixelImage upright, DefishSettings settings, double x, double y, out double rx, out double ry)
        {
            int turn = SensorTurn(upright);
            bool swap = Orientation.SwapsAxes(turn);
            int w = swap ? upright.Height : upright.Width, h = swap ? upright.Width : upright.Height;
            DefishCanvas c;
            int ow, oh;
            Plan(g, w, h, settings, out c, out ow, out oh);
            double sx, sy, u, v;
            Orientation.MapPoint(turn, upright.Width, upright.Height, x, y, out sx, out sy);
            bool ok = g.Unmap((sx + 0.5) / w, (sy + 0.5) / h, out u, out v);
            double ox = (u - c.U0) / c.Width * ow - 0.5, oy = (v - c.V0) / c.Height * oh - 0.5;
            Orientation.MapPoint(Orientation.Inverse(turn), ow, oh, ox, oy, out rx, out ry);
            return ok && ox >= -0.5 && oy >= -0.5 && ox <= ow - 0.5 && oy <= oh - 0.5;
        }

        /// <summary>The photo point (upright, pixel centres) shown at (<paramref name="rx"/>, <paramref name="ry"/>) of the result of <see cref="Process"/>.</summary>
        public static void ResultToPhoto(DefishGeometry g, PixelImage upright, DefishSettings settings, double rx, double ry, out double x, out double y)
        {
            int turn = SensorTurn(upright);
            bool swap = Orientation.SwapsAxes(turn);
            int w = swap ? upright.Height : upright.Width, h = swap ? upright.Width : upright.Height;
            DefishCanvas c;
            int ow, oh;
            Plan(g, w, h, settings, out c, out ow, out oh);
            int rw = swap ? oh : ow, rh = swap ? ow : oh;   // the result as returned (upright)
            double ox, oy, s, t, a, b, cc, d;
            Orientation.MapPoint(turn, rw, rh, rx, ry, out ox, out oy);
            g.Map(c.U0 + (ox + 0.5) / ow * c.Width, c.V0 + (oy + 0.5) / oh * c.Height, out s, out t, out a, out b, out cc, out d);
            double sx = s * w - 0.5, sy = t * h - 0.5;
            Orientation.MapPoint(Orientation.Inverse(turn), w, h, sx, sy, out x, out y);
        }

        /// <summary>True when the photo's shape differs from the 3:2 frame the projection was measured on by more than 2%.</summary>
        public static bool ShapeDiffers(DefishGeometry g, PixelImage upright)
        {
            var m = g.Mapping;
            if (m == null) return false;        // a model works in the lens' own terms: any frame shape is defished correctly
            double want = m.SourceWidth > 0 && m.SourceHeight > 0 ? (double)m.SourceWidth / m.SourceHeight : 1.5;
            double have = (double)Math.Max(upright.Width, upright.Height) / Math.Min(upright.Width, upright.Height);
            return Math.Abs(have / want - 1) > 0.02;
        }
    }
}
