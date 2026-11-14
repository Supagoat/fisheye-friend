using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;

namespace FisheyeFriend.Core.Testing
{
    /// <summary>
    /// The defish conformance suite, shared by the Unity tests and the command-line bench: the golden points of the
    /// model (computed independently by Tools/TestData/golden.py), the test images, the render cases with their golden
    /// results, and the comparisons. Every defish implementation (the CPU one, the GPU one) is held to it.
    /// </summary>
    public static class Conformance
    {
        /// <summary>Frame width the mapping errors are measured on, in pixels (a 64 MP 3:2 photo is about 9800 wide).</summary>
        public const int ReferenceWidth = 8000;

        // ================================================================== golden points

        public sealed class GoldenCase
        {
            public string Name;
            public double Aspect;
            public DefishAdjustment Adjustment;
            /// <summary>u, v (result point), s, t (photo point), all normalised.</summary>
            public double[][] Points;

            public override string ToString()
            {
                return Name;
            }
        }

        public static List<GoldenCase> LoadGolden(string json)
        {
            var root = Json.Parse(json);
            var list = Json.Get(root, "cases") as List<object>;
            if (list == null) throw new FormatException("no cases");
            var cases = new List<GoldenCase>();
            foreach (var o in list)
            {
                var pts = (List<object>)Json.Get(o, "points");
                var c = new GoldenCase
                {
                    Name = Json.String(o, "name", "?"),
                    Aspect = Json.Number(o, "aspect", 2.0 / 3.0),
                    Adjustment = new DefishAdjustment
                    {
                        Strength = Json.Number(o, "strength", 1),
                        AmountX = Json.Number(o, "amountX", 1),
                        AmountY = Json.Number(o, "amountY", 1),
                        ShapeX = Doubles(Json.Get(o, "shapeX")),
                        ShapeY = Doubles(Json.Get(o, "shapeY")),
                    },
                    Points = new double[pts.Count][],
                };
                for (int i = 0; i < pts.Count; i++) c.Points[i] = Doubles(pts[i]);
                cases.Add(c);
            }
            return cases;
        }

        private static double[] Doubles(object o)
        {
            var l = (List<object>)o;
            var r = new double[l.Count];
            for (int i = 0; i < r.Length; i++) r[i] = (double)l[i];
            return r;
        }

        public struct MappingError
        {
            public double MaxPx, RmsPx;
            public double WorstU, WorstV;
            public int Points;

            public override string ToString()
            {
                var ci = CultureInfo.InvariantCulture;
                return "max " + MaxPx.ToString("0.000000", ci) + " px, rms " + RmsPx.ToString("0.000000", ci) + " px over " + Points +
                       " points (on a frame " + ReferenceWidth + " px wide; worst at u " + WorstU.ToString("0.000", ci) + ", v " + WorstV.ToString("0.000", ci) + ")";
            }
        }

        /// <summary>
        /// How far a mapping's photo points are from the golden ones, in pixels of a frame <see cref="ReferenceWidth"/>
        /// wide. <paramref name="map"/> fills s and t for (u, v), all normalised.
        /// </summary>
        public static MappingError CompareMapping(GoldenCase c, Func<double, double, double[]> map)
        {
            double w = ReferenceWidth, h = ReferenceWidth * c.Aspect, sum = 0, max = -1, wu = 0, wv = 0;
            foreach (var p in c.Points)
            {
                var r = map(p[0], p[1]);
                double dx = (r[0] - p[2]) * w, dy = (r[1] - p[3]) * h, d = Math.Sqrt(dx * dx + dy * dy);
                if (double.IsNaN(d)) d = double.PositiveInfinity;
                sum += d * d;
                if (d > max)
                {
                    max = d;
                    wu = p[0];
                    wv = p[1];
                }
            }
            return new MappingError { MaxPx = max, RmsPx = Math.Sqrt(sum / c.Points.Length), WorstU = wu, WorstV = wv, Points = c.Points.Length };
        }

        // ================================================================== test images

        /// <summary>
        /// A deterministic test picture with the kinds of content photos have: smooth gradients, soft and hard edges,
        /// thin lines, fine periodic detail up to the pixel scale, and a little noise. Integer arithmetic only, so it is
        /// identical on every machine. With 4 channels the alpha varies too (a gradient, and a transparent hole).
        /// </summary>
        public static PixelImage Scene(int w, int h, int channels, int depth, uint seed = 1)
        {
            var img = new PixelImage(w, h, channels, depth);
            int max = depth == 8 ? 255 : 65535;
            uint rng = seed * 2654435761u + 12345u;
            int cc = channels == 4 || channels == 2 ? channels - 1 : channels;
            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    // everything in 1/4096 units of full scale
                    long gx = (long)x * 4096 / Math.Max(1, w - 1), gy = (long)y * 4096 / Math.Max(1, h - 1);
                    for (int c = 0; c < cc; c++)
                    {
                        long v = c == 0 ? 600 + gx * 5 / 8 : c == 1 ? 900 + gy / 2 : 3000 - (gx + gy) / 3;
                        // discs: soft (quadratic falloff) and hard
                        long dx = x - w * 3 / 10, dy = y - h * 2 / 5, r2 = dx * dx + dy * dy, R = Math.Max(4, h / 5);
                        if (r2 < R * R) v += (R * R - r2) * (c == 2 ? 1800 : 900) / (R * R);
                        dx = x - w * 7 / 10; dy = y - h * 3 / 5;
                        if (dx * dx + dy * dy < (long)(h / 7) * (h / 7)) v = c == 1 ? 3600 : 400;
                        // a rectangle with hard edges
                        if (x > w / 10 && x < w / 4 && y > h * 3 / 5 && y < h * 9 / 10) v = 4000 - c * 1300;
                        // thin lines: one vertical, one horizontal, one diagonal
                        if (x == w / 2 || y == h / 6 || x - y == w / 3) v = c == 0 ? 3900 : 200;
                        // fine detail: a triangle wave whose period shrinks to 2 pixels towards the right, in a band
                        if (y > h * 3 / 4 && y < h * 7 / 8)
                        {
                            long period = 2 + (w - 1 - x) * 14 / Math.Max(1, w);
                            long phase = (x % period) * 2 * 1600 / period;
                            v += (phase < 1600 ? phase : 3200 - phase) - 800;
                        }
                        // noise: about +-0.4% of full scale
                        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
                        v += (long)(rng % 33) - 16;
                        v = Math.Max(0, Math.Min(4095, v));
                        Set(img, x, y, c, (int)((v * max + 2047) / 4095));
                    }
                    if (cc < channels)
                    {
                        long a = 4095 - gy * 2 / 5;
                        long hx = x - w / 2, hy = y - h / 2;
                        if (hx * hx + hy * hy < (long)(h / 8) * (h / 8)) a = 0;
                        Set(img, x, y, cc, (int)((a * max + 2047) / 4095));
                    }
                }
            return img;
        }

        /// <summary>Uniform noise in every channel: the hardest case for a resampler (all frequencies at full contrast).</summary>
        public static PixelImage Noise(int w, int h, int channels, int depth, uint seed = 7)
        {
            var img = new PixelImage(w, h, channels, depth);
            uint rng = seed * 2654435761u + 977u;
            int n = w * h * channels;
            for (int i = 0; i < n; i++)
            {
                rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
                if (depth == 8) img.Data8[i] = (byte)(rng >> 24);
                else img.Data16[i] = (ushort)(rng >> 16);
            }
            if (channels == 4)
                for (int i = 3; i < n; i += 4)
                    if (depth == 8) img.Data8[i] = 255; else img.Data16[i] = 65535;
            return img;
        }

        private static void Set(PixelImage img, int x, int y, int c, int v)
        {
            int i = (y * img.Width + x) * img.Channels + c;
            if (img.BitDepth == 8) img.Data8[i] = (byte)v;
            else img.Data16[i] = (ushort)v;
        }

        // ================================================================== render cases

        public sealed class RenderCase
        {
            public string Id;
            public string Description;
            public int Width, Height, Channels, Depth;
            public bool NoiseImage;
            /// <summary>The camera's orientation tag (the image is given upright, as displayed); 0 = none.</summary>
            public int Orientation;
            public DefishSettings Settings;

            /// <summary>The photo, upright, with its orientation tag.</summary>
            public PixelImage MakeImage()
            {
                var img = NoiseImage ? Noise(Width, Height, Channels, Depth) : Scene(Width, Height, Channels, Depth);
                img.Metadata.SourceOrientation = Orientation;
                return img;
            }

            public string GoldenFile => "render-" + Id + ".png";

            /// <summary>The frame shape the geometry is made for (the lens' landscape frame).</summary>
            public double SensorAspect => Defisher.SensorAspect(MakeHeader());

            private PixelImage MakeHeader()
            {
                var img = new PixelImage(Width, Height, 1, 8);
                img.Metadata.SourceOrientation = Orientation;
                return img;
            }

            public override string ToString()
            {
                return Id;
            }
        }

        private static DefishSettings S(Action<DefishSettings> f)
        {
            var s = new DefishSettings { OutputSize = DefishOutputSize.FullDetail };
            f(s);
            s.Validate();
            return s;
        }

        /// <summary>The render cases: every image layout and every option the app offers, at several adjustments.</summary>
        public static readonly RenderCase[] RenderCases =
        {
            new RenderCase { Id = "rgb8-default", Description = "RGB 8-bit, the defaults (full detail, full frame, Lanczos 4, sigmoidal, aperture)",
                Width = 300, Height = 200, Channels = 3, Depth = 8, Settings = S(s => { }) },
            new RenderCase { Id = "rgb8-adjusted-crop", Description = "RGB 8-bit, strength 80%, horizontal 70%, vertical 120%, both curves, cropped",
                Width = 300, Height = 200, Channels = 3, Depth = 8, Settings = S(s =>
                {
                    s.Edges = DefishEdges.Crop;
                    s.Adjust = new DefishAdjustment { Strength = 0.8, AmountX = 0.7, AmountY = 1.2, ShapeX = new[] { 1, 1, 0.7, 0.5 }, ShapeY = new[] { 1, 1.1, 1.2, 1.3 } };
                }) },
            new RenderCase { Id = "rgb8-portrait", Description = "RGB 8-bit portrait (orientation tag 6), strength 110%",
                Width = 200, Height = 300, Channels = 3, Depth = 8, Orientation = 6, Settings = S(s => s.Adjust = new DefishAdjustment { Strength = 1.1 }) },
            new RenderCase { Id = "gray8-16x9-lanczos3-gamma", Description = "Grey 8-bit 16:9, Lanczos 3, gamma light, photo size, vertical 60%",
                Width = 320, Height = 180, Channels = 1, Depth = 8, Settings = S(s =>
                {
                    s.Kernel = WarpKernel.TensorLanczos;
                    s.Radius = 3;
                    s.Light = LightSpace.Gamma;
                    s.OutputSize = DefishOutputSize.SameAsPhoto;
                    s.Adjust = new DefishAdjustment { AmountY = 0.6 };
                }) },
            new RenderCase { Id = "rgba8-bicubic-linear-aa", Description = "RGBA 8-bit, bicubic, linear light, anti-aliasing, no aperture correction",
                Width = 240, Height = 160, Channels = 4, Depth = 8, Settings = S(s =>
                {
                    s.Kernel = WarpKernel.Bicubic;
                    s.Light = LightSpace.Linear;
                    s.Antialias = true;
                    s.ApertureCorrection = false;
                    s.Adjust = new DefishAdjustment { ShapeX = new[] { 1, 1.2, 0.8, 0.6 }, ShapeY = new[] { 1, 1.2, 0.8, 0.6 } };
                }) },
            new RenderCase { Id = "rgb16-strong", Description = "RGB 16-bit, strength 115%, horizontal 120%",
                Width = 300, Height = 200, Channels = 3, Depth = 16, Settings = S(s => s.Adjust = new DefishAdjustment { Strength = 1.15, AmountX = 1.2 }) },
            new RenderCase { Id = "rgba16-crop-sigmoidal", Description = "RGBA 16-bit, cropped, sigmoidal light (premultiplied alpha), strength 60%",
                Width = 210, Height = 140, Channels = 4, Depth = 16, Settings = S(s =>
                {
                    s.Edges = DefishEdges.Crop;
                    s.Adjust = new DefishAdjustment { Strength = 0.6 };
                }) },
            new RenderCase { Id = "rgb8-noise", Description = "RGB 8-bit noise (the worst case for any resampler), the defaults",
                Width = 150, Height = 100, Channels = 3, Depth = 8, NoiseImage = true, Settings = S(s => { }) },
            new RenderCase { Id = "rgb8-none", Description = "RGB 8-bit with no correction (both amounts 0%): the photo itself",
                Width = 300, Height = 200, Channels = 3, Depth = 8, Settings = S(s => s.Adjust = new DefishAdjustment { AmountX = 0, AmountY = 0 }) },
        };

        // ================================================================== image comparison

        public sealed class ImageDiff
        {
            public bool SameLayout;
            public string Layout;
            /// <summary>Largest difference of one sample, in the image's own units (8-bit: 0..255; 16-bit: 0..65535).</summary>
            public int Max;
            public long Differing, Samples;
            public double MeanAbs;
            /// <summary>Number of samples differing by 0, 1, 2, 3, 4 and more.</summary>
            public long[] Histogram = new long[5];
            public int MaxValue;

            public double IdenticalFraction => Samples > 0 ? 1 - (double)Differing / Samples : 0;

            public override string ToString()
            {
                if (!SameLayout) return "different images: " + Layout;
                var ci = CultureInfo.InvariantCulture;
                var sb = new StringBuilder();
                sb.Append(Layout).Append(": ").Append((IdenticalFraction * 100).ToString("0.0000", ci)).Append("% of samples identical, max difference ")
                  .Append(Max).Append(" of ").Append(MaxValue).Append(", mean ").Append(MeanAbs.ToString("0.000000", ci)).Append(" (by difference 0/1/2/3/4+: ");
                for (int i = 0; i < Histogram.Length; i++) sb.Append(i == 0 ? "" : "/").Append(Histogram[i]);
                return sb.Append(")").ToString();
            }
        }

        public static ImageDiff Compare(PixelImage a, PixelImage b)
        {
            var d = new ImageDiff
            {
                SameLayout = a.Width == b.Width && a.Height == b.Height && a.Channels == b.Channels && a.BitDepth == b.BitDepth,
                Layout = a.Width + "x" + a.Height + "x" + a.Channels + " " + a.BitDepth + "-bit" +
                         (a.Width == b.Width && a.Height == b.Height && a.Channels == b.Channels && a.BitDepth == b.BitDepth ? "" :
                             " vs " + b.Width + "x" + b.Height + "x" + b.Channels + " " + b.BitDepth + "-bit"),
                MaxValue = a.MaxValue,
            };
            if (!d.SameLayout) return d;
            long n = (long)a.Width * a.Height * a.Channels;
            double sum = 0;
            for (long i = 0; i < n; i++)
            {
                int va = a.BitDepth == 8 ? a.Data8[i] : a.Data16[i];
                int vb = b.BitDepth == 8 ? b.Data8[i] : b.Data16[i];
                int e = Math.Abs(va - vb);
                if (e > d.Max) d.Max = e;
                if (e != 0) d.Differing++;
                d.Histogram[Math.Min(4, e)]++;
                sum += e;
            }
            d.Samples = n;
            d.MeanAbs = sum / n;
            return d;
        }
    }
}
