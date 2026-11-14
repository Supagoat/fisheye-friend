using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Color;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;
using System.Threading.Tasks;

namespace FisheyeFriendBench
{
    /// <summary>
    /// Command-line bench for the Fisheye Friend core library (compile it together with Assets/FisheyeFriend/Core).
    /// Images are read and written as PNG. Commands:
    ///   selftest &lt;pairs.json&gt; [model.json]            correctness checks of the projection, the model, the warp and the defish pipeline
    ///   analyze &lt;pairs.json&gt; [W H]                    local scale of the projection across the frame
    ///   grid &lt;pairs.json&gt;                             interpolation error of the grid for several grid sizes
    ///   exp-prep &lt;pairs&gt; &lt;scene.png&gt; &lt;prefix&gt; &lt;W&gt; &lt;H&gt; &lt;blur&gt;  ground truth and fisheye capture of a scene
    ///   exp-run &lt;pairs&gt; &lt;capture.png&gt; &lt;prefix&gt; [methods]   defish a capture with each method (see score.py)
    ///   ff-replica &lt;pairs&gt; &lt;photo.png&gt; &lt;out.png&gt; [fix]   the old Fisheye Friend CPU renderer, exactly
    ///   defish &lt;pairs&gt; &lt;photo.png&gt; &lt;out.png&gt; [radius] [full] [crop] [noap] [aa]   the app's pipeline, timed
    ///   prof &lt;pairs&gt; &lt;W&gt; &lt;H&gt;                         timing of the map and the warp
    /// </summary>
    public static class Program
    {
        public static int Main(string[] args)
        {
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            System.Threading.Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
            try
            {
                switch (args.Length > 0 ? args[0] : "")
                {
                    case "analyze": return Analyze(args);
                    case "grid": return GridCmd(args);
                    case "exp-prep": return ExpPrep(args);
                    case "exp-run": return ExpRun(args);
                    case "ff-replica": return ReplicaCmd(args);
                    case "defish": return DefishCmd(args);
                    case "model-defish": return ModelDefishCmd(args);
                    case "prof": return ProfCmd(args);
                    case "selftest": return SelfTest(args);
                    case "jac": return JacCmd(args);
                    case "tps-dump": return TpsDumpCmd(args);
                    case "model-check": return ModelCheckCmd(args);
                    default:
                        Console.WriteLine("commands: selftest | analyze | grid | exp-prep | exp-run | ff-replica | defish | prof (see the comment in Bench.cs)");
                        return 1;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                return 2;
            }
        }

        private static DefishGeometry s_geo;
        private static DefishGeometry Geo(string pairs)
        {
            return s_geo ?? (s_geo = new DefishGeometry(HemiMapping.Parse(File.ReadAllText(pairs))));
        }

        private static double ToLin(double v) { return v <= 0.04045 ? v / 12.92 : Math.Pow((v + 0.055) / 1.055, 2.4); }
        private static double FromLin(double l) { l = Math.Max(0, Math.Min(1, l)); return l <= 0.0031308 ? l * 12.92 : 1.055 * Math.Pow(l, 1 / 2.4) - 0.055; }

        /// <summary>Solves geometry(u, v) = (s, t) for the defished point of a photo point (Newton, from (s, t)).</summary>
        private static bool Invert(DefishGeometry g, double s, double t, ref double u, ref double v)
        {
            for (int it = 0; it < 30; it++)
            {
                double ms, mt, a, b, c, d;
                g.Map(u, v, out ms, out mt, out a, out b, out c, out d);
                double es = ms - s, et = mt - t;
                if (Math.Abs(es) < 1e-12 && Math.Abs(et) < 1e-12) return true;
                double det = a * d - b * c;
                if (Math.Abs(det) < 1e-12) return false;
                double du = (d * es - b * et) / det, dv = (a * et - c * es) / det;
                double lim = 0.05;
                if (Math.Abs(du) > lim || Math.Abs(dv) > lim) { double k = lim / Math.Max(Math.Abs(du), Math.Abs(dv)); du *= k; dv *= k; }
                u -= du; v -= dv;
            }
            return false;
        }

        // exp-prep <pairs> <scene4x.png> <outPrefix> <Wo> <Ho> <blur>: ground truth (scene integrated over the defished output
        // pixels) and a simulated fisheye capture (scene through the inverse map, integrated over 4x4 samples per sensor pixel).
        private static int ExpPrep(string[] a)
        {
            var g = Geo(a[1]);
            var scene = PngCodec.Decode(File.ReadAllBytes(a[2]));
            string prefix = a[3];
            int Wo = int.Parse(a[4]), Ho = int.Parse(a[5]);
            double blur = double.Parse(a[6], CultureInfo.InvariantCulture);
            int SW = scene.Width, SH = scene.Height, ch = scene.Channels;
            if (SW != 4 * Wo || SH != 4 * Ho) throw new Exception("scene must be 4x the output");
            var lin = new float[SW * SH * 3];
            for (int i = 0; i < SW * SH; i++)
                for (int c = 0; c < 3; c++)
                    lin[i * 3 + c] = (float)ToLin((scene.BitDepth == 8 ? scene.Data8[i * ch + c] / 255.0 : scene.Data16[i * ch + c] / 65535.0));
            if (!File.Exists(prefix + "_gt.png"))
            {
                var gt = new PixelImage(Wo, Ho, 3, 8);
                Parallel.For(0, Ho, y =>
                {
                    for (int x = 0; x < Wo; x++)
                        for (int c = 0; c < 3; c++)
                        {
                            double acc = 0;
                            for (int j = 0; j < 4; j++) for (int i = 0; i < 4; i++) acc += lin[((4 * y + j) * SW + 4 * x + i) * 3 + c];
                            gt.Data8[(y * Wo + x) * 3 + c] = (byte)Math.Round(FromLin(acc / 16) * 255);
                        }
                });
                File.WriteAllBytes(prefix + "_gt.png", PngCodec.Encode(gt));
            }
            var cap = new double[Wo * Ho * 3];
            int failures = 0;
            Parallel.For(0, Ho, y =>
            {
                double gu = 0.5, gv = (y + 0.5) / Ho;
                for (int x = 0; x < Wo; x++)
                {
                    for (int j = 0; j < 4; j++)
                        for (int i = 0; i < 4; i++)
                        {
                            double s = (x + (i + 0.5) / 4) / Wo, t = (y + (j + 0.5) / 4) / Ho;
                            double u = s, v = t;
                            if (!Invert(g, s, t, ref u, ref v)) System.Threading.Interlocked.Increment(ref failures);
                            // scene pixel centres at whole numbers; bilinear, clamped
                            double px = Math.Max(0, Math.Min(SW - 1.001, u * SW - 0.5)), py = Math.Max(0, Math.Min(SH - 1.001, v * SH - 0.5));
                            int x0 = (int)px, y0 = (int)py;
                            double fx = px - x0, fy = py - y0;
                            for (int c = 0; c < 3; c++)
                            {
                                double v00 = lin[(y0 * SW + x0) * 3 + c], v10 = lin[(y0 * SW + x0 + 1) * 3 + c];
                                double v01 = lin[((y0 + 1) * SW + x0) * 3 + c], v11 = lin[((y0 + 1) * SW + x0 + 1) * 3 + c];
                                cap[(y * Wo + x) * 3 + c] += ((v00 * (1 - fx) + v10 * fx) * (1 - fy) + (v01 * (1 - fx) + v11 * fx) * fy) / 16;
                            }
                        }
                }
            });
            if (blur > 0) GaussianBlur(cap, Wo, Ho, blur);
            var capImg = new PixelImage(Wo, Ho, 3, 8);
            for (int i = 0; i < cap.Length; i++) capImg.Data8[i] = (byte)Math.Round(FromLin(cap[i]) * 255);
            File.WriteAllBytes(prefix + "_cap.png", PngCodec.Encode(capImg));
            Console.WriteLine(Path.GetFileName(prefix) + ": capture " + Wo + "x" + Ho + (failures > 0 ? ", " + failures + " samples did not invert (outside the map)" : ""));
            return 0;
        }

        private static void GaussianBlur(double[] img, int w, int h, double sigma)
        {
            int r = (int)Math.Ceiling(3 * sigma);
            var k = new double[2 * r + 1];
            double sum = 0;
            for (int i = -r; i <= r; i++) { k[i + r] = Math.Exp(-i * i / (2 * sigma * sigma)); sum += k[i + r]; }
            for (int i = 0; i < k.Length; i++) k[i] /= sum;
            var tmp = new double[img.Length];
            Parallel.For(0, h, y => { for (int x = 0; x < w; x++) for (int c = 0; c < 3; c++) { double acc = 0; for (int i = -r; i <= r; i++) acc += k[i + r] * img[(y * w + Math.Max(0, Math.Min(w - 1, x + i))) * 3 + c]; tmp[(y * w + x) * 3 + c] = acc; } });
            Parallel.For(0, h, y => { for (int x = 0; x < w; x++) for (int c = 0; c < 3; c++) { double acc = 0; for (int i = -r; i <= r; i++) acc += k[i + r] * tmp[(Math.Max(0, Math.Min(h - 1, y + i)) * w + x) * 3 + c]; img[(y * w + x) * 3 + c] = acc; } });
        }

        /// <summary>Fisheye Friend's CPU renderer, line for line (radial Lanczos 3, 8-bit gamma, clamped edges), on a given map.</summary>
        private static PixelImage FisheyeFriendOriginal(PixelImage src, WarpMap map, int dw, int dh, bool fixOffset)
        {
            int sw = src.Width, sh = src.Height, ch = src.Channels;
            var dst = new PixelImage(dw, dh, 3, 8);
            const double A = 3.0;
            Func<double, double> sinc = t => Math.Abs(t) < 1e-8 ? 1.0 : Math.Sin(Math.PI * t) / (Math.PI * t);
            Parallel.For(0, dh, y =>
            {
                for (int x = 0; x < dw; x++)
                {
                    double sx, sy, a, b, c, d;
                    map.Map(x, y, out sx, out sy, out a, out b, out c, out d);   // = map*sw - 0.5, as in the original
                    double off = fixOffset ? 0 : 0.5;
                    int xmin = (int)Math.Floor(sx - A + 1), xmax = (int)Math.Ceiling(sx + A);
                    int ymin = (int)Math.Floor(sy - A + 1), ymax = (int)Math.Ceiling(sy + A);
                    double W = 0, R = 0, G = 0, B = 0;
                    for (int yy = ymin; yy <= ymax; yy++)
                        for (int xx = xmin; xx <= xmax; xx++)
                        {
                            double dx = (xx + off) - sx, dy = (yy + off) - sy;
                            double r = Math.Sqrt(dx * dx + dy * dy);
                            if (r >= A) continue;
                            double w = sinc(r) * sinc(r / A);
                            int px = Math.Max(0, Math.Min(sw - 1, xx)), py = Math.Max(0, Math.Min(sh - 1, yy));
                            int o = (py * sw + px) * ch;
                            W += w;
                            R += w * src.Data8[o] / 255.0; G += w * src.Data8[o + 1] / 255.0; B += w * src.Data8[o + 2] / 255.0;
                        }
                    int q = (y * dw + x) * 3;
                    dst.Data8[q] = (byte)Math.Max(0, Math.Min(255, Math.Round(R / W * 255)));
                    dst.Data8[q + 1] = (byte)Math.Max(0, Math.Min(255, Math.Round(G / W * 255)));
                    dst.Data8[q + 2] = (byte)Math.Max(0, Math.Min(255, Math.Round(B / W * 255)));
                }
            });
            return dst;
        }

        // exp-run <pairs> <capture.png> <outPrefix> [methods,comma]: defish the capture with each method.
        private static int ExpRun(string[] a)
        {
            var g = Geo(a[1]);
            var cap = PngCodec.Decode(File.ReadAllBytes(a[2]));
            string prefix = a[3];
            var only = a.Length > 4 ? new HashSet<string>(a[4].Split(',')) : null;
            int W = cap.Width, H = cap.Height;
            var map = new DefishWarpMap(g, W, H, W, H);
            var times = new List<string>();
            Action<string, Func<PixelImage>> run = (name, f) =>
            {
                if (only != null && !only.Contains(name)) return;
                var sw = Stopwatch.StartNew();
                var r = f();
                times.Add(name + " " + sw.ElapsedMilliseconds);
                File.WriteAllBytes(prefix + "_" + name + ".png", PngCodec.Encode(r));
            };
            Func<WarpKernel, double, LightSpace, bool, PixelImage> warp = (k, rad, l, aa) =>
                WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = k, Radius = rad, Light = l, Antialias = aa, ExtendEdges = true, OutputBitDepth = 8 });
            string maskPath = Path.Combine(Path.GetDirectoryName(prefix), "mask_" + W + "x" + H + ".png");
            if (!File.Exists(maskPath))
            {
                var mask = new PixelImage(W, H, 1, 8);
                for (int y = 0; y < H; y++)
                    for (int x = 0; x < W; x++)
                    {
                        double sx, sy, q0, q1, q2, q3;
                        map.Map(x, y, out sx, out sy, out q0, out q1, out q2, out q3);
                        mask.Data8[y * W + x] = (byte)(sx >= 6 && sx <= W - 7 && sy >= 6 && sy <= H - 7 ? 255 : 0);
                    }
                File.WriteAllBytes(maskPath, PngCodec.Encode(mask));
            }
            run("ff", () => FisheyeFriendOriginal(cap, map, W, H, false));
            run("ff_fixed", () => FisheyeFriendOriginal(cap, map, W, H, true));
            run("t3s", () => warp(WarpKernel.TensorLanczos, 3, LightSpace.Sigmoidal, true));
            run("t3g", () => warp(WarpKernel.TensorLanczos, 3, LightSpace.Gamma, true));
            run("t3l", () => warp(WarpKernel.TensorLanczos, 3, LightSpace.Linear, true));
            run("t3s_noaa", () => warp(WarpKernel.TensorLanczos, 3, LightSpace.Sigmoidal, false));
            run("t2s", () => warp(WarpKernel.TensorLanczos, 2, LightSpace.Sigmoidal, true));
            run("t4s", () => warp(WarpKernel.TensorLanczos, 4, LightSpace.Sigmoidal, true));
            run("r3s", () => warp(WarpKernel.RadialLanczos, 3, LightSpace.Sigmoidal, true));
            run("bcs", () => warp(WarpKernel.Bicubic, 2, LightSpace.Sigmoidal, true));
            Func<double, double, LightSpace, PixelImage> warpK = (rad, k, l) =>
                WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = WarpKernel.TensorLanczos, Radius = rad, Light = l, Antialias = true, AntialiasStrength = k, ExtendEdges = true, OutputBitDepth = 8 });
            run("t3s_aa25", () => warpK(3, 0.25, LightSpace.Sigmoidal));
            run("t3s_aa50", () => warpK(3, 0.5, LightSpace.Sigmoidal));
            run("t3s_aa75", () => warpK(3, 0.75, LightSpace.Sigmoidal));
            run("t4s_noaa", () => warp(WarpKernel.TensorLanczos, 4, LightSpace.Sigmoidal, false));
            run("t3g_noaa", () => warp(WarpKernel.TensorLanczos, 3, LightSpace.Gamma, false));
            Func<double, double, double, PixelImage> warpAp = (rad, ap, k) =>
                WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = WarpKernel.TensorLanczos, Radius = rad, Light = LightSpace.Sigmoidal, Antialias = k > 0, AntialiasStrength = k, Aperture = ap, ExtendEdges = true, OutputBitDepth = 8 });
            foreach (var gain in new[] { 0.5, 1.0, 1.5, 2.0, 3.0 })
            {
                string gs = ((int)(gain * 100)).ToString(); double gv = gain;
                run("t3s_ap" + gs, () => warpAp(3, gv, 0));
                run("t4s_ap" + gs, () => warpAp(4, gv, 0));
                run("t4s_app" + gs, () => WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = WarpKernel.TensorLanczos, Radius = 4, Light = LightSpace.Sigmoidal, Antialias = false, Aperture = gv, ApertureSharpenOnly = true, ExtendEdges = true, OutputBitDepth = 8 }));
                run("t3s_app" + gs, () => WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = WarpKernel.TensorLanczos, Radius = 3, Light = LightSpace.Sigmoidal, Antialias = false, Aperture = gv, ApertureSharpenOnly = true, ExtendEdges = true, OutputBitDepth = 8 }));
            }
            run("bcs_app100", () => WarpResampler.Warp(cap, W, H, map, new WarpOptions { Kernel = WarpKernel.Bicubic, Radius = 2, Light = LightSpace.Sigmoidal, Antialias = false, Aperture = 1, ApertureSharpenOnly = true, ExtendEdges = true, OutputBitDepth = 8 }));
            Console.WriteLine(Path.GetFileName(prefix) + ": " + string.Join(", ", times) + " ms");
            return 0;
        }

        // defish <pairs> <photo.png|jpg-as-png> <out.png> [radius] [full] [crop]: the app's pipeline at full size, timed
        private static int DefishCmd(string[] a)
        {
            var sw = Stopwatch.StartNew();
            var g = Geo(a[1]);
            Console.WriteLine("geometry " + sw.ElapsedMilliseconds + " ms");
            var src = ImageIO.MakeUpright(PngCodec.Decode(File.ReadAllBytes(a[2])));
            var s = new DefishSettings();
            for (int i = 4; i < a.Length; i++)
            {
                double r;
                if (double.TryParse(a[i], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out r)) s.Radius = r;
                else if (a[i] == "full") s.OutputSize = DefishOutputSize.FullDetail;
                else if (a[i] == "crop") s.Edges = DefishEdges.Crop;
                else if (a[i] == "noap") s.ApertureCorrection = false;
                else if (a[i] == "aa") s.Antialias = true;
            }
            sw.Restart();
            var r2 = Defisher.Process(src, g, s);
            Console.WriteLine("defished " + src.Width + "x" + src.Height + " -> " + r2.Width + "x" + r2.Height + " (" + s + ") in " + sw.ElapsedMilliseconds + " ms");
            File.WriteAllBytes(a[3], PngCodec.Encode(r2));
            return 0;
        }

        // model-defish <photo.png> <out.png> [strength=1] [h=1] [v=1] [shapex=a,b,c,d] [shapey=a,b,c,d] [preview=W] [full] [crop]:
        // the photo defished with the built-in model and an adjustment (or, with preview=W, the panel's quick picture W wide)
        private static int ModelDefishCmd(string[] a)
        {
            var src = ImageIO.MakeUpright(PngCodec.Decode(File.ReadAllBytes(a[1])));
            var s = new DefishSettings { OutputSize = DefishOutputSize.SameAsPhoto };
            int preview = 0;
            var ci = CultureInfo.InvariantCulture;
            for (int i = 3; i < a.Length; i++)
            {
                string k = a[i], v = "";
                int eq = k.IndexOf('=');
                if (eq > 0) { v = k.Substring(eq + 1); k = k.Substring(0, eq); }
                switch (k)
                {
                    case "strength": s.Adjust.Strength = double.Parse(v, ci); break;
                    case "h": s.Adjust.AmountX = double.Parse(v, ci); break;
                    case "v": s.Adjust.AmountY = double.Parse(v, ci); break;
                    case "shapex": s.Adjust.ShapeX = Array.ConvertAll(v.Split(','), x => double.Parse(x, ci)); break;
                    case "shapey": s.Adjust.ShapeY = Array.ConvertAll(v.Split(','), x => double.Parse(x, ci)); break;
                    case "preview": preview = int.Parse(v, ci); break;
                    case "full": s.OutputSize = DefishOutputSize.FullDetail; break;
                    case "crop": s.Edges = DefishEdges.Crop; break;
                }
            }
            s.Validate();
            var sw = Stopwatch.StartNew();
            double aspect = Defisher.SensorAspect(src);
            if (preview > 0)
            {
                int n = Math.Max(1, Math.Max(src.Width, src.Height) / 1400), rw, rh;
                var reduced = DisplayReducer.ReduceToRgba32(src, n, out rw, out rh);
                long tReduce = sw.ElapsedMilliseconds;
                sw.Restart();
                var gq = DefishGeometry.FromModel(DefishModel.Hemi, s.Adjust, aspect, 129);
                long tGeo = sw.ElapsedMilliseconds;
                sw.Restart();
                int ph = (int)Math.Round(preview * (double)src.Height / src.Width);
                var px = DefishPreview.Render(reduced, rw, rh, Defisher.SensorTurn(src), gq, preview, ph);
                Console.WriteLine("preview " + preview + "x" + ph + ": reduce " + tReduce + " ms (once per photo), geometry " + tGeo + " ms, picture " + sw.ElapsedMilliseconds + " ms");
                var img = new PixelImage(preview, ph, 4, 8);
                Buffer.BlockCopy(px, 0, img.Data8, 0, px.Length);
                File.WriteAllBytes(a[2], PngCodec.Encode(img));
                return 0;
            }
            var g = DefishGeometry.FromModel(DefishModel.Hemi, s.Adjust, aspect);
            Console.WriteLine("geometry " + sw.ElapsedMilliseconds + " ms" + (g.Folds ? " - FOLDS" : "") + ", min area scale " + g.MinAreaScale.ToString("0.000"));
            sw.Restart();
            var r = Defisher.Process(src, g, s);
            Console.WriteLine("defished " + src.Width + "x" + src.Height + " -> " + r.Width + "x" + r.Height + " (" + s + ") in " + sw.ElapsedMilliseconds + " ms");
            File.WriteAllBytes(a[2], PngCodec.Encode(r));
            return 0;
        }

        private sealed class IdentityMap : WarpMap
        {
            public override void Map(double x, double y, out double sx, out double sy, out double j00, out double j01, out double j10, out double j11)
            {
                sx = x; sy = y; j00 = 1; j01 = 0; j10 = 0; j11 = 1;
            }
        }

        private static int _fails;

        private static void Check(bool ok, string what)
        {
            Console.WriteLine((ok ? "PASS  " : "FAIL  ") + what);
            if (!ok) _fails++;
        }

        private static PixelImage Noise(int w, int h, int ch, int depth, int seed)
        {
            var img = new PixelImage(w, h, ch, depth);
            var rnd = new Random(seed);
            if (depth == 8) rnd.NextBytes(img.Data8);
            else for (int i = 0; i < img.Data16.Length; i++) img.Data16[i] = (ushort)rnd.Next(65536);
            if (ch == 4)
                for (int i = 3; i < w * h * 4; i += 4)
                    if (depth == 8) img.Data8[i] = 255; else img.Data16[i] = 65535;
            return img;
        }

        private static int MaxDiff(PixelImage a, PixelImage b, int channels)
        {
            int m = 0;
            for (int i = 0; i < a.Width * a.Height; i++)
                for (int c = 0; c < channels; c++)
                {
                    int va = a.BitDepth == 8 ? a.Data8[i * a.Channels + c] : a.Data16[i * a.Channels + c];
                    int vb = b.BitDepth == 8 ? b.Data8[i * b.Channels + c] : b.Data16[i * b.Channels + c];
                    m = Math.Max(m, Math.Abs(va - vb));
                }
            return m;
        }

        // selftest <pairs>
        private static int SelfTest(string[] a)
        {
            _fails = 0;
            var mapping = HemiMapping.Parse(File.ReadAllText(a[1]));
            Check(mapping.Pairs.Count == 596 && mapping.SourceWidth == 6000 && mapping.SourceHeight == 4000, "mapping: 596 pairs on 6000 x 4000 (" + mapping.Pairs.Count + ")");
            var g = new DefishGeometry(mapping);

            // the spline goes through the measured points (lambda smooths a little)
            ThinPlateSpline tx, ty;
            mapping.FitInverse(1e-5, out tx, out ty);
            double worst = 0;
            foreach (var p in mapping.Pairs)
                worst = Math.Max(worst, Math.Max(Math.Abs(tx.Eval(p.DstX, p.DstY) - p.SrcX) * 6000, Math.Abs(ty.Eval(p.DstX, p.DstY) - p.SrcY) * 4000));
            Check(worst < 2, "spline residual at the 596 points " + worst.ToString("0.000") + " px");

            double ge = g.MaxGridError(7952, 5304, 20000, 7);
            Check(ge < 0.01, "grid error at 7952 x 5304: " + ge.ToString("0.0000") + " px");

            // derivatives match finite differences
            double dErr = 0;
            var rnd = new Random(5);
            for (int k = 0; k < 2000; k++)
            {
                double u = rnd.NextDouble(), v = rnd.NextDouble(), s0, t0, su, sv, tu, tv, s1, t1, s2, t2, q0, q1, q2, q3;
                g.Map(u, v, out s0, out t0, out su, out sv, out tu, out tv);
                const double h = 1e-6;
                g.Map(u + h, v, out s1, out t1, out q0, out q1, out q2, out q3);
                g.Map(u, v + h, out s2, out t2, out q0, out q1, out q2, out q3);
                dErr = Math.Max(dErr, Math.Max(Math.Max(Math.Abs((s1 - s0) / h - su), Math.Abs((t1 - t0) / h - tu)), Math.Max(Math.Abs((s2 - s0) / h - sv), Math.Abs((t2 - t0) / h - tv))));
            }
            Check(dErr < 1e-3, "Jacobian vs finite differences " + dErr.ToString("0.0e0"));

            double ue = 0;
            for (int k = 0; k < 2000; k++)
            {
                double u = 0.02 + 0.96 * rnd.NextDouble(), v = 0.02 + 0.96 * rnd.NextDouble(), s0, t0, q0, q1, q2, q3, uu, vv;
                g.Map(u, v, out s0, out t0, out q0, out q1, out q2, out q3);
                bool ok = g.Unmap(s0, t0, out uu, out vv);
                ue = Math.Max(ue, ok ? Math.Max(Math.Abs(uu - u), Math.Abs(vv - v)) : 1);
            }
            Check(ue < 1e-8, "Unmap(Map) round trip " + ue.ToString("0.0e0"));

            // aperture filter: none at 1:1, sharpening where enlarging, smoothing where shrinking
            Check(Math.Abs(WarpResampler.ApertureCoefficient(1)) < 1e-6 && WarpResampler.ApertureCoefficient(0.55) > 0.05 &&
                  WarpResampler.ApertureCoefficient(1.19) < 0, "aperture coefficients " + WarpResampler.ApertureCoefficient(0.55).ToString("0.000") + " / " +
                  WarpResampler.ApertureCoefficient(1.19).ToString("0.000"));

            // identity map is exact for every kernel, light space and depth, with and without aperture correction
            foreach (var depth in new[] { 8, 16 })
                foreach (var ch in new[] { 1, 3, 4 })
                {
                    var img = Noise(37, 23, ch, depth, ch * depth);
                    int m = 0;
                    foreach (var k in new[] { WarpKernel.TensorLanczos, WarpKernel.Bicubic })
                        foreach (var l in new[] { LightSpace.Gamma, LightSpace.Sigmoidal, LightSpace.Linear })
                            foreach (var ap in new[] { 0.0, 1.0 })
                            {
                                var r = WarpResampler.Warp(img, 37, 23, new IdentityMap(), new WarpOptions { Kernel = k, Radius = k == WarpKernel.Bicubic ? 2 : 4, Light = l, Aperture = ap, ApertureSharpenOnly = true, Antialias = false, ExtendEdges = true });
                                int md = MaxDiff(img, r, ch);
                                if (md > m && Environment.GetEnvironmentVariable("FFDEBUG") != null) Console.WriteLine("   " + k + " " + l + " ap " + ap + ": " + md);
                                m = Math.Max(m, md);
                            }
                    Check(m <= (depth == 8 ? 0 : 1), "identity warp " + ch + " ch " + depth + "-bit: max diff " + m);
                }

            // fast (RGB) and general (RGBA) paths agree
            {
                var rgb = Noise(160, 107, 3, 8, 11);
                var rgba = new PixelImage(160, 107, 4, 8);
                for (int i = 0; i < 160 * 107; i++)
                {
                    for (int c = 0; c < 3; c++) rgba.Data8[i * 4 + c] = rgb.Data8[i * 3 + c];
                    rgba.Data8[i * 4 + 3] = 255;
                }
                var s = new DefishSettings { Light = LightSpace.Gamma };
                var r1 = Defisher.Process(rgb, g, s);
                var r2 = Defisher.Process(rgba, g, s);
                Check(MaxDiff(r1, r2, 3) <= 1, "RGB fast path vs RGBA path: max diff " + MaxDiff(r1, r2, 3));
            }

            // flat field stays flat, in every light space
            {
                int m = 0;
                foreach (var l in new[] { LightSpace.Gamma, LightSpace.Sigmoidal, LightSpace.Linear })
                {
                    var flat = new PixelImage(150, 100, 3, 8);
                    for (int i = 0; i < flat.Data8.Length; i += 3) { flat.Data8[i] = 200; flat.Data8[i + 1] = 90; flat.Data8[i + 2] = 17; }
                    var r = Defisher.Process(flat, g, new DefishSettings { Light = l });
                    for (int i = 0; i < r.Data8.Length; i += 3)
                        m = Math.Max(m, Math.Max(Math.Abs(r.Data8[i] - 200), Math.Max(Math.Abs(r.Data8[i + 1] - 90), Math.Abs(r.Data8[i + 2] - 17))));
                }
                Check(m <= 1, "flat field: max deviation " + m);
            }

            // sizes
            {
                var land = Noise(300, 200, 3, 8, 21);
                int w, h;
                Defisher.OutputSize(g, land, new DefishSettings(), out w, out h);
                Check(w == 300 && h == 200, "same-as-photo size " + w + "x" + h);
                Defisher.OutputSize(g, land, new DefishSettings { OutputSize = DefishOutputSize.FullDetail }, out w, out h);
                double k = Defisher.CentreScale(g, 300, 200);
                Check(Math.Abs(k - 1.19) < 0.01 && w == (int)Math.Round(300 * k) && h == (int)Math.Round(200 * k), "full-detail size " + w + "x" + h + " (scale " + k.ToString("0.0000") + ")");
                var c = Defisher.CropCanvas(g, 7952, 5304);
                Check(c.Width > 0.995 && c.Width < 1, "crop keeps " + (c.Width * 100).ToString("0.000") + "% of each side");
                // every pixel of a cropped output reads from inside the photo
                int W = 7952, H = 5304, ow, oh;
                DefishCanvas cv;
                var cs = new DefishSettings { Edges = DefishEdges.Crop };
                Defisher.Plan(g, W, H, cs, out cv, out ow, out oh);
                var map = new DefishWarpMap(g, W, H, ow, oh, cv.U0, cv.V0, cv.U1, cv.V1);
                double minIn = 1e9;
                for (int i = 0; i < ow; i += 1)
                    foreach (var y in new[] { 0, oh - 1 })
                    {
                        double sx, sy, q0, q1, q2, q3;
                        map.Map(i, y, out sx, out sy, out q0, out q1, out q2, out q3);
                        minIn = Math.Min(minIn, Math.Min(Math.Min(sx, W - 1 - sx), Math.Min(sy, H - 1 - sy)));
                    }
                for (int j = 0; j < oh; j += 1)
                    foreach (var x in new[] { 0, ow - 1 })
                    {
                        double sx, sy, q0, q1, q2, q3;
                        map.Map(x, j, out sx, out sy, out q0, out q1, out q2, out q3);
                        minIn = Math.Min(minIn, Math.Min(Math.Min(sx, W - 1 - sx), Math.Min(sy, H - 1 - sy)));
                    }
                Check(minIn >= -1e-6, "cropped output stays inside the photo (closest " + minIn.ToString("0.000") + " px from the edge centres)");
                Defisher.Plan(g, W, H, new DefishSettings(), out cv, out ow, out oh);
                map = new DefishWarpMap(g, W, H, ow, oh);
                minIn = 1e9;
                int firstOut = -1, lastOut = -1, rowsOut = 0;
                for (int i = 0; i < ow; i++)
                {
                    double sx, sy, q0, q1, q2, q3;
                    map.Map(i, 0, out sx, out sy, out q0, out q1, out q2, out q3);
                    minIn = Math.Min(minIn, sy);
                    if (sy < -0.5) { if (firstOut < 0) firstOut = i; lastOut = i; }
                }
                for (int j = 0; j < 20; j++)
                {
                    double sx, sy, q0, q1, q2, q3;
                    map.Map(ow / 2, j, out sx, out sy, out q0, out q1, out q2, out q3);
                    if (sy < -0.5) rowsOut++;
                }
                Check(minIn < 0, "uncropped output reaches past the photo's top edge (" + minIn.ToString("0.00") + " px; columns " + firstOut + "-" + lastOut + " of " + ow + ", " + rowsOut + " rows at the centre)");
            }

            // portraits are defished in the sensor frame: the same as turning, defishing and turning back
            {
                var land = Noise(240, 160, 3, 8, 31);
                var s = new DefishSettings();
                var rl = Defisher.Process(land, g, s);
                foreach (var o in new[] { Orientation.Rotate90Clockwise, Orientation.Rotate90CounterClockwise })
                {
                    var portrait = Orientation.Apply(land, o);              // the photo as displayed
                    portrait.Metadata.SourceOrientation = o;                // the camera's tag
                    var rp = Defisher.Process(portrait, g, s);
                    var expect = Orientation.Apply(rl, o);
                    Check(rp.Width == expect.Width && rp.Height == expect.Height && MaxDiff(rp, expect, 3) == 0, "portrait (tag " + o + ") equals turned landscape result");
                }
                var untagged = Orientation.Apply(land, Orientation.Rotate90CounterClockwise);
                untagged.Metadata.SourceOrientation = 0;
                var ru = Defisher.Process(untagged, g, s);
                Check(MaxDiff(ru, Orientation.Apply(rl, Orientation.Rotate90CounterClockwise), 3) == 0, "untagged portrait is defished turned a quarter clockwise");
                var r16 = Defisher.Process(Noise(90, 60, 3, 16, 41), g, s);
                Check(r16.BitDepth == 16, "16-bit stays 16-bit");
            }

            // view mapping round trips, landscape and portrait, in each output mode
            {
                double worstRt = 0;
                foreach (var o in new[] { 0, Orientation.Rotate90Clockwise, Orientation.Rotate90CounterClockwise })
                    foreach (var os in new[] { DefishOutputSize.SameAsPhoto, DefishOutputSize.FullDetail })
                        foreach (var ed in new[] { DefishEdges.Extend, DefishEdges.Crop })
                        {
                            var img = o == 0 ? new PixelImage(600, 400, 3, 8) : new PixelImage(400, 600, 3, 8);
                            img.Metadata.SourceOrientation = o;
                            var s = new DefishSettings { OutputSize = os, Edges = ed };
                            for (int k = 0; k < 200; k++)
                            {
                                double x = 20 + rnd.NextDouble() * (img.Width - 40), y = 20 + rnd.NextDouble() * (img.Height - 40), rx, ry, bx, by;
                                if (!Defisher.PhotoToResult(g, img, s, x, y, out rx, out ry)) continue;
                                Defisher.ResultToPhoto(g, img, s, rx, ry, out bx, out by);
                                worstRt = Math.Max(worstRt, Math.Max(Math.Abs(bx - x), Math.Abs(by - y)));
                            }
                        }
                Check(worstRt < 1e-4, "photo -> result -> photo round trip " + worstRt.ToString("0.0e0") + " px");
                // the mapping agrees with the warp: a dot in the photo lands where PhotoToResult says
                var dot = new PixelImage(600, 400, 1, 8);
                double dx = 452.3, dy = 97.6;
                for (int y = 0; y < 400; y++)
                    for (int x = 0; x < 600; x++)
                    {
                        double r2 = (x - dx) * (x - dx) + (y - dy) * (y - dy);
                        dot.Data8[y * 600 + x] = (byte)Math.Round(250 * Math.Exp(-r2 / 18));
                    }
                var rd = Defisher.Process(dot, g, new DefishSettings { Light = LightSpace.Gamma });
                double sw = 0, mx = 0, my = 0;
                for (int y = 0; y < rd.Height; y++)
                    for (int x = 0; x < rd.Width; x++)
                    {
                        double v = rd.Data8[y * rd.Width + x];
                        if (v < 20) continue;
                        sw += v; mx += v * x; my += v * y;
                    }
                double ex, ey;
                Defisher.PhotoToResult(g, dot, new DefishSettings(), dx, dy, out ex, out ey);
                Check(Math.Abs(mx / sw - ex) < 0.3 && Math.Abs(my / sw - ey) < 0.3, "dot lands at the mapped point (" + (mx / sw).ToString("0.00") + ", " + (my / sw).ToString("0.00") + ") vs (" + ex.ToString("0.00") + ", " + ey.ToString("0.00") + ")");
            }

            // the mapping file's other spelling (pixels instead of normalised) parses the same
            {
                var alt = HemiMapping.Parse("{\"source\":{\"width\":100,\"height\":50},\"target\":{\"width\":100,\"height\":50},\"pairs\":[" +
                                            "{\"src\":{\"x\":10,\"y\":5},\"dst\":{\"xn\":0.1,\"yn\":0.1}},{\"src\":{\"x\":90,\"y\":5},\"dst\":{\"x\":90,\"y\":5}}," +
                                            "{\"src\":{\"x\":50,\"y\":45},\"dst\":{\"x\":50,\"y\":45}}]}");
                Check(alt.Pairs.Count == 3 && Math.Abs(alt.Pairs[0].SrcX - 0.1) < 1e-12 && Math.Abs(alt.Pairs[2].DstY - 0.9) < 1e-12, "pixel-unit pairs parse");
                bool threw = false;
                try { HemiMapping.Parse("{\"pairs\": [1, 2"); } catch (FormatException) { threw = true; }
                Check(threw, "malformed mapping is reported");
            }

            ModelSelfTest(g, a.Length > 2 ? a[2] : null);

            Console.WriteLine(_fails == 0 ? "ALL PASSED" : _fails + " FAILED");
            return _fails == 0 ? 0 : 1;
        }

        /// <summary>The defish model and its adjustment: against the point pairs' splines, and its promises.</summary>
        private static void ModelSelfTest(DefishGeometry splines, string modelPath)
        {
            const double A = 2.0 / 3.0;
            var model = modelPath != null ? DefishModel.Parse(File.ReadAllText(modelPath)) : DefishModel.Hemi;

            // the shipped file is the built-in model, and JSON round-trips exactly
            {
                var hemi = DefishModel.Hemi;
                var back = DefishModel.Parse(model.ToJson());
                Check(SameModel(model, hemi), "model file " + (modelPath != null ? Path.GetFileName(modelPath) : "(built in)") + " has the built-in coefficients");
                Check(SameModel(model, back), "model JSON round trip");
                bool threw = false;
                try { DefishModel.Parse("{\"format\": \"fisheye-friend-defish-model\", \"lens\": {}}"); } catch (FormatException) { threw = true; }
                Check(threw, "incomplete model file is reported");
            }

            var g = DefishGeometry.FromModel(model, new DefishAdjustment(), A);
            // matches the point pairs (the hand-measured defish) over the frame
            {
                double sum = 0, mx = 0;
                int n = 0;
                for (int j = 0; j < 60; j++)
                    for (int i = 0; i < 90; i++)
                    {
                        double u = (i + 0.5) / 90, v = (j + 0.5) / 60, s1, t1, s2, t2, q0, q1, q2, q3;
                        g.Map(u, v, out s1, out t1, out q0, out q1, out q2, out q3);
                        splines.Map(u, v, out s2, out t2, out q0, out q1, out q2, out q3);
                        double d = Math.Sqrt(Math.Pow((s1 - s2) * 6000, 2) + Math.Pow((t1 - t2) * 4000, 2));
                        sum += d * d;
                        n++;
                        mx = Math.Max(mx, d);
                    }
                double rms = Math.Sqrt(sum / n);
                Check(rms < 1.2 && mx < 3, "model vs point-pair splines on 6000 x 4000: rms " + rms.ToString("0.00") + " px, max " + mx.ToString("0.00") + " px");
            }
            double ge = g.MaxGridError(7952, 5304, 20000, 7);
            Check(ge < 0.01, "model grid error at 7952 x 5304: " + ge.ToString("0.0000") + " px");
            {
                var rnd = new Random(9);
                double ue = 0;
                for (int k = 0; k < 2000; k++)
                {
                    double u = 0.02 + 0.96 * rnd.NextDouble(), v = 0.02 + 0.96 * rnd.NextDouble(), s0, t0, q0, q1, q2, q3, uu, vv;
                    g.Map(u, v, out s0, out t0, out q0, out q1, out q2, out q3);
                    bool ok = g.Unmap(s0, t0, out uu, out vv);
                    ue = Math.Max(ue, ok ? Math.Max(Math.Abs(uu - u), Math.Abs(vv - v)) : 1);
                }
                Check(ue < 1e-8, "model Unmap(Map) round trip " + ue.ToString("0.0e0"));
            }

            // symmetric left-right and top-bottom, and vertical lines come out straight (x depends on the column only)
            {
                var mm = new DefishModelMap(model, new DefishAdjustment { Strength = 0.9, AmountX = 1.2, AmountY = 0.7, ShapeX = new[] { 1, 1.1, 0.8, 0.6 } }, A);
                double sym = 0;
                var rnd = new Random(4);
                for (int k = 0; k < 500; k++)
                {
                    double u = rnd.NextDouble(), v = rnd.NextDouble(), s, t, s2, t2;
                    mm.Reset(); mm.Map(u, v, out s, out t);
                    mm.Reset(); mm.Map(1 - u, v, out s2, out t2);
                    sym = Math.Max(sym, Math.Max(Math.Abs(s + s2 - 1), Math.Abs(t - t2)));
                    mm.Reset(); mm.Map(u, 1 - v, out s2, out t2);
                    sym = Math.Max(sym, Math.Max(Math.Abs(s - s2), Math.Abs(t + t2 - 1)));
                }
                Check(sym < 1e-9, "adjusted model is symmetric: " + sym.ToString("0.0e0"));
            }

            // no adjustment at all (both amounts 0) leaves the photo as it is, at any strength and curve
            {
                double idm = 0;
                foreach (var st in new[] { 0.5, 1.0, 1.15 })
                {
                    var mm = new DefishModelMap(model, new DefishAdjustment { Strength = st, AmountX = 0, AmountY = 0, ShapeY = new[] { 2.0, 0, 1, 2 } }, A);
                    for (int j = 0; j <= 10; j++)
                        for (int i = 0; i <= 10; i++)
                        {
                            double s, t;
                            mm.Map(i / 10.0, j / 10.0, out s, out t);
                            idm = Math.Max(idm, Math.Max(Math.Abs(s - i / 10.0), Math.Abs(t - j / 10.0)));
                        }
                }
                Check(idm < 1e-12, "amounts at 0 leave the photo unchanged: " + idm.ToString("0.0e0"));
            }

            // at any strength the result's left and right edges read from the photo's edges (the width is kept; the measured
            // defish itself reaches about 1.7 px past them)
            {
                double we = 0;
                foreach (var st in new[] { 0.5, 0.75, 1.0, 1.15 })
                {
                    var mm = new DefishModelMap(model, new DefishAdjustment { Strength = st }, A);
                    double s, t;
                    mm.Map(1, 0.5, out s, out t);
                    we = Math.Max(we, Math.Abs(s - 1) * 6000);
                    mm.Reset();
                    mm.Map(0, 0.5, out s, out t);
                    we = Math.Max(we, Math.Abs(s) * 6000);
                }
                Check(we < 3, "every strength keeps the width: the edges read within " + we.ToString("0.0") + " px of the photo's edges (6000 wide)");
            }

            // a stronger lens sees more; the strength changes the defish the way it should (more correction at the edges)
            {
                double f0 = model.DiagonalFieldOfView(A, 0.5), f1 = model.DiagonalFieldOfView(A), f2 = model.DiagonalFieldOfView(A, 1.15);
                Check(f0 < f1 && f1 < f2 && Math.Abs(f1 - 177) < 1, "field of view " + f0.ToString("0") + "° / " + f1.ToString("0") + "° / " + f2.ToString("0") + "° at 50 / 100 / 115%");
                // how far the defish moves the frame's corner and a point two thirds of the way along the diagonal
                double prev = 0;
                bool mono = true;
                string seen = "";
                foreach (var st in new[] { 0.5, 0.8, 1.0, 1.15 })
                {
                    var mm = new DefishModelMap(model, new DefishAdjustment { Strength = st }, A);
                    double s, t;
                    mm.Map(1.0 / 6, 1.0 / 6, out s, out t);
                    double d = Math.Sqrt(Math.Pow((s - 1.0 / 6) * 6000, 2) + Math.Pow((t - 1.0 / 6) * 4000, 2));
                    mono &= d > prev;
                    prev = d;
                    seen += (seen.Length > 0 ? " / " : "") + d.ToString("0");
                }
                Check(mono, "a stronger lens gets more correction: a point two thirds out moves " + seen + " px at 50 / 80 / 100 / 115%");
            }

            // no folds anywhere in the range the interface offers with flat curves; folds are found when they happen
            {
                bool any = false;
                double worst = double.MaxValue;
                foreach (var st in new[] { DefishAdjustment.MinStrength, 1.0, DefishAdjustment.MaxStrength })
                    foreach (var ax in new[] { 0.0, 0.5, 1.0, DefishAdjustment.MaxAmount })
                        foreach (var ay in new[] { 0.0, 1.0, DefishAdjustment.MaxAmount })
                        {
                            var gq = DefishGeometry.FromModel(model, new DefishAdjustment { Strength = st, AmountX = ax, AmountY = ay }, A, 257);
                            any |= gq.Folds;
                            worst = Math.Min(worst, gq.MinAreaScale);
                        }
                Check(!any, "no folds over strength 50-115%, amounts 0-150% (smallest area scale " + worst.ToString("0.000") + ")");
                // ordinary curves are safe; a zig-zag between 0 and 200% does fold, and is found
                bool curved = false;
                double least = double.MaxValue;
                foreach (var shape in new[] { new[] { 1, 1, 0.5, 0 }, new[] { 0.5, 0.7, 1, 1.4 }, new[] { 1, 1, 1.5, 1.5 }, new[] { 2.0, 0, 2, 0 } })
                {
                    var gc = DefishGeometry.FromModel(model, new DefishAdjustment { ShapeX = shape, ShapeY = shape }, A, 257);
                    curved |= gc.Folds;
                    least = Math.Min(least, gc.MinAreaScale);
                }
                Check(!curved, "no folds with ordinary curves at 100% (smallest area scale " + least.ToString("0.000") + ")");
                var zig = DefishGeometry.FromModel(model, new DefishAdjustment { ShapeX = new[] { 0.0, 2, 0, 2 }, ShapeY = new[] { 0.0, 2, 0, 2 } }, A, 257);
                Check(zig.Folds, "a zig-zag curve (0 / 200 / 0 / 200%) is found to fold (area scale " + zig.MinAreaScale.ToString("0.000") + ")");
                var fold = DefishGeometry.FromModel(model, new DefishAdjustment { AmountX = 1.5, AmountY = 1.5, ShapeX = new[] { 0, 0, 2.0, 2 } }, A, 257);
                Console.WriteLine("      (steepest curve at 150%: area scale " + fold.MinAreaScale.ToString("0.000") + (fold.Folds ? ", folds" : "") + ")");
            }

            // the shape curves pass through their points and never overshoot them
            {
                var sh = new[] { 1.0, 1.0, 0.5, 0.0 };
                bool through = Math.Abs(DefishAdjustment.ShapeAt(sh, 0) - 1) < 1e-12 && Math.Abs(DefishAdjustment.ShapeAt(sh, 1.0 / 3) - 1) < 1e-12 &&
                               Math.Abs(DefishAdjustment.ShapeAt(sh, 2.0 / 3) - 0.5) < 1e-12 && Math.Abs(DefishAdjustment.ShapeAt(sh, 1)) < 1e-12;
                bool mono = true;
                double last = 2;
                for (int i = 0; i <= 300; i++)
                {
                    double v = DefishAdjustment.ShapeAt(sh, i / 300.0);
                    mono &= v <= last + 1e-12 && v >= -1e-12 && v <= 1 + 1e-12;
                    last = v;
                }
                Check(through && mono, "shape curve passes through its points without overshoot");
                var bad = new DefishAdjustment { Strength = 9, AmountX = -1, AmountY = double.NaN, ShapeX = new[] { 5.0 }, ShapeY = null };
                bad.Validate();
                Check(bad.Strength == DefishAdjustment.MaxStrength && bad.AmountX == 0 && bad.AmountY == 1 && bad.ShapeX.Length == 4 && bad.ShapeX[0] == 2 &&
                      bad.ShapeX[3] == 1 && bad.ShapeY.Length == 4, "out-of-range adjustments are clamped");
                var k1 = new DefishAdjustment { AmountX = 0.9 };
                Check(!k1.Same(new DefishAdjustment()) && k1.Same(k1.Clone()) && new DefishAdjustment().IsDefault && !k1.IsDefault, "adjustment identity (Key, Same, IsDefault)");
                var s1 = new DefishSettings();
                var s2 = s1.Clone();
                s2.Adjust.ShapeY[2] = 0.7;
                Check(!s1.Same(s2) && s1.Adjust.ShapeY[2] == 1, "settings clone deeply and tell adjustments apart");
            }

            // the panel's preview: with no correction it reproduces the photo; a portrait is the turned landscape picture
            {
                var src = Noise(120, 80, 4, 8, 51);
                for (int i = 3; i < src.Data8.Length; i += 4) src.Data8[i] = 255;
                var id = DefishGeometry.FromModel(model, new DefishAdjustment { AmountX = 0, AmountY = 0 }, A, 129);
                var pic = DefishPreview.Render(src.Data8, 120, 80, Orientation.Normal, id, 120, 80);
                int md = 0;
                for (int i = 0; i < pic.Length; i++) md = Math.Max(md, Math.Abs(pic[i] - src.Data8[i]));
                Check(md <= 1, "preview with no correction reproduces the photo (max diff " + md + ")");

                var gq = DefishGeometry.FromModel(model, new DefishAdjustment(), A, 129);
                var land = DefishPreview.Render(src.Data8, 120, 80, Orientation.Normal, gq, 120, 80);
                var portrait = Orientation.Apply(src, Orientation.Rotate90Clockwise);   // 80 x 120, as displayed
                var picP = DefishPreview.Render(portrait.Data8, 80, 120, Orientation.Rotate90CounterClockwise, gq, 80, 120);
                var landImg = new PixelImage(120, 80, 4, 8);
                Buffer.BlockCopy(land, 0, landImg.Data8, 0, land.Length);
                var expect = Orientation.Apply(landImg, Orientation.Rotate90Clockwise);
                int mp = 0;
                for (int i = 0; i < picP.Length; i++) mp = Math.Max(mp, Math.Abs(picP[i] - expect.Data8[i]));
                Check(mp <= 1, "portrait preview is the turned landscape preview (max diff " + mp + ")");
                var flipped = DefishPreview.Render(src.Data8, 120, 80, Orientation.Normal, gq, 120, 80, true);
                bool rows = true;
                for (int x = 0; x < 120 * 4; x++) rows &= flipped[x] == land[79 * 120 * 4 + x];
                Check(rows, "preview bottom-up rows for the texture");
                var sw = Stopwatch.StartNew();
                var big = Noise(1400, 933, 4, 8, 52);
                var gq2 = DefishGeometry.FromModel(model, new DefishAdjustment { Strength = 0.9 }, A, 129);
                DefishPreview.Render(big.Data8, 1400, 933, Orientation.Normal, gq2, 660, 440, true);
                Console.WriteLine("      (preview geometry + 660 x 440 picture: " + sw.ElapsedMilliseconds + " ms)");
            }

            // the defisher uses the model: full detail follows the strength, and a portrait is the turned landscape
            {
                var land = Noise(240, 160, 3, 8, 61);
                var s = new DefishSettings { Adjust = new DefishAdjustment { Strength = 0.8, AmountY = 0.6 } };
                var gm = DefishGeometry.FromModel(model, s.Adjust, A, 257);
                var rl = Defisher.Process(land, gm, s);
                var portrait = Orientation.Apply(land, Orientation.Rotate90Clockwise);
                portrait.Metadata.SourceOrientation = Orientation.Rotate90Clockwise;
                Check(Math.Abs(Defisher.SensorAspect(portrait) - A) < 1e-12, "a portrait's lens frame is landscape");
                var rp = Defisher.Process(portrait, gm, s);
                Check(MaxDiff(rp, Orientation.Apply(rl, Orientation.Rotate90Clockwise), 3) == 0, "adjusted portrait equals the turned landscape result");
                double k0 = Defisher.CentreScale(DefishGeometry.FromModel(model, new DefishAdjustment { Strength = 0.5 }, A, 257), 600, 400);
                double k1 = Defisher.CentreScale(g, 600, 400);
                double k2 = Defisher.CentreScale(DefishGeometry.FromModel(model, new DefishAdjustment { Strength = 1.15 }, A, 257), 600, 400);
                Check(k0 < k1 && k1 < k2 && Math.Abs(k1 - 1.19) < 0.02, "centre scale grows with strength: " + k0.ToString("0.000") + " / " + k1.ToString("0.000") + " / " + k2.ToString("0.000"));
                var folding = new DefishGeometry[] { DefishGeometry.FromModel(model, new DefishAdjustment { AmountX = 1.5, AmountY = 1.5, ShapeX = new[] { 0, 0, 2.0, 2 } }, A, 257) };
                if (folding[0].Folds)
                {
                    bool threw = false;
                    try { Defisher.Process(land, folding[0], s); } catch (InvalidOperationException) { threw = true; }
                    Check(threw, "folding settings are refused, not rendered");
                }
            }
        }

        private static bool SameModel(DefishModel a, DefishModel b)
        {
            bool same = a.Theta.Length == b.Theta.Length && a.Horizontal.Length == b.Horizontal.Length &&
                        a.Vertical.GetLength(0) == b.Vertical.GetLength(0) && a.Vertical.GetLength(1) == b.Vertical.GetLength(1) &&
                        a.Scale == b.Scale && a.Panini == b.Panini && a.PaniniVertical == b.PaniniVertical;
            if (!same) return false;
            for (int i = 0; i < a.Theta.Length; i++) same &= a.Theta[i] == b.Theta[i];
            for (int i = 0; i < a.Horizontal.Length; i++) same &= a.Horizontal[i] == b.Horizontal[i];
            for (int i = 0; i < a.Vertical.GetLength(0); i++)
                for (int j = 0; j < a.Vertical.GetLength(1); j++) same &= a.Vertical[i, j] == b.Vertical[i, j];
            return same;
        }

        // model-check <pairs> <model.json> <ref_points.json>: the C# model against the Python reference, the splines, and its speed
        private static int ModelCheckCmd(string[] a)
        {
            var model = DefishModel.Parse(File.ReadAllText(a[2]));
            var refs = Json.Parse(File.ReadAllText(a[3])) as List<object>;
            double worst = 0;
            foreach (var o in refs)
            {
                var r = (List<object>)o;
                double x = (double)r[0], y = (double)r[1];
                var adj = new DefishAdjustment { Strength = (double)r[2], AmountX = (double)r[3], AmountY = (double)r[4] };
                var mm = new DefishModelMap(model, adj, 2.0 / 3.0);
                double s, t;
                mm.Map(x + 0.5, y / (2.0 / 3.0) + 0.5, out s, out t);
                double ex = (double)r[5], ey = (double)r[6];
                worst = Math.Max(worst, Math.Max(Math.Abs((s - 0.5) - ex), Math.Abs((t - 0.5) * 2 / 3 - ey)) * 6000);
            }
            Console.WriteLine("C# vs Python reference: " + refs.Count + " points, max difference " + worst.ToString("0.000000") + " px (6000 wide)");
            var sw = Stopwatch.StartNew();
            var g = DefishGeometry.FromModel(model, new DefishAdjustment(), 2.0 / 3.0);
            Console.WriteLine("model grid built in " + sw.ElapsedMilliseconds + " ms (" + g.GridNodes + " nodes), grid error " + g.MaxGridError(7952, 5304, 20000, 3).ToString("0.0000") +
                              " px on 7952 x 5304, min area scale " + g.MinAreaScale.ToString("0.000"));
            sw.Restart();
            var quick = DefishGeometry.FromModel(model, new DefishAdjustment(), 2.0 / 3.0, 129);
            Console.WriteLine("preview grid (129 across) in " + sw.ElapsedMilliseconds + " ms, grid error " + quick.MaxGridError(7952, 5304, 5000, 3).ToString("0.00") + " px");
            var tps = Geo(a[1]);
            double sum = 0, mx = 0; int n = 0; double wx = 0, wy = 0;
            for (int j = 0; j < 80; j++)
                for (int i = 0; i < 120; i++)
                {
                    double u = (i + 0.5) / 120, v = (j + 0.5) / 80, s1, t1, s2, t2;
                    g.MapExact(u, v, out s1, out t1);
                    tps.MapExact(u, v, out s2, out t2);
                    double d = Math.Sqrt(Math.Pow((s1 - s2) * 6000, 2) + Math.Pow((t1 - t2) * 4000, 2));
                    sum += d * d; n++;
                    if (d > mx) { mx = d; wx = u; wy = v; }
                }
            Console.WriteLine("model vs splines over the frame (photo px, 6000 x 4000): rms " + Math.Sqrt(sum / n).ToString("0.00") + ", max " + mx.ToString("0.00") +
                              " at (" + (wx * 6000).ToString("0") + ", " + (wy * 4000).ToString("0") + ")");
            foreach (var st in new[] { 0.5, 0.7, 0.85, 1.0, 1.1, 1.15 })
            {
                var gs = DefishGeometry.FromModel(model, new DefishAdjustment { Strength = st }, 2.0 / 3.0, 257);
                Console.WriteLine("  strength " + st.ToString("0.00") + ": diagonal FOV " + model.DiagonalFieldOfView(2.0 / 3.0, st).ToString("0") + " deg, centre scale " +
                                  Defisher.CentreScale(gs, 6000, 4000).ToString("0.000") + ", min area scale " + gs.MinAreaScale.ToString("0.000"));
            }
            return 0;
        }

        // tps-dump <pairs> <nx> <ny> <out.csv>: the exact spline map (output -> source, normalised) on an nx x ny grid of the output
        private static int TpsDumpCmd(string[] a)
        {
            var g = Geo(a[1]);
            int nx = int.Parse(a[2]), ny = int.Parse(a[3]);
            var sb = new System.Text.StringBuilder("u,v,s,t\n");
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                {
                    double u = (i + 0.5) / nx, v = (j + 0.5) / ny, s, t;
                    g.MapExact(u, v, out s, out t);
                    sb.Append(u.ToString("R", CultureInfo.InvariantCulture)).Append(',').Append(v.ToString("R", CultureInfo.InvariantCulture)).Append(',')
                      .Append(s.ToString("R", CultureInfo.InvariantCulture)).Append(',').Append(t.ToString("R", CultureInfo.InvariantCulture)).Append('\n');
                }
            File.WriteAllText(a[4], sb.ToString());
            return 0;
        }

        // jac <pairs> <W> <H> <x> <y>: source position and Jacobian d(src px)/d(out px) at an output pixel (same-size output)
        private static int JacCmd(string[] a)
        {
            var g = Geo(a[1]);
            int W = int.Parse(a[2]), H = int.Parse(a[3]);
            var map = new DefishWarpMap(g, W, H, W, H);
            double sx, sy, j00, j01, j10, j11;
            map.Map(double.Parse(a[4], CultureInfo.InvariantCulture), double.Parse(a[5], CultureInfo.InvariantCulture), out sx, out sy, out j00, out j01, out j10, out j11);
            Console.WriteLine(string.Format(CultureInfo.InvariantCulture, "src ({0:0.00}, {1:0.00})  J = [[{2:0.0000}, {3:0.0000}], [{4:0.0000}, {5:0.0000}]]", sx, sy, j00, j01, j10, j11));
            return 0;
        }

        // prof <pairs> <W> <H>: cost of the map evaluation alone, and of the warp on a synthetic image, best of 3
        private static int ProfCmd(string[] a)
        {
            var g = Geo(a[1]);
            int W = int.Parse(a[2]), H = int.Parse(a[3]);
            var map = new DefishWarpMap(g, W, H, W, H);
            var img = new PixelImage(W, H, 3, 8);
            new Random(3).NextBytes(img.Data8);
            for (int rep = 0; rep < 3; rep++)
            {
                var sw = Stopwatch.StartNew();
                double sink = 0;
                Parallel.For(0, H, y =>
                {
                    double loc = 0;
                    for (int x = 0; x < W; x++)
                    {
                        double sx, sy, j0, j1, j2, j3;
                        map.Map(x, y, out sx, out sy, out j0, out j1, out j2, out j3);
                        loc += sx + j3;
                    }
                    lock (map) sink += loc;
                });
                long tMap = sw.ElapsedMilliseconds;
                sw.Restart();
                var s = new DefishSettings();
                WarpResampler.Warp(img, W, H, map, Defisher.WarpOptionsFor(img, s));
                long t4 = sw.ElapsedMilliseconds;
                sw.Restart();
                s.Radius = 3;
                WarpResampler.Warp(img, W, H, map, Defisher.WarpOptionsFor(img, s));
                long t3 = sw.ElapsedMilliseconds;
                sw.Restart();
                s.Radius = 4; s.ApertureCorrection = false;
                WarpResampler.Warp(img, W, H, map, Defisher.WarpOptionsFor(img, s));
                long t4n = sw.ElapsedMilliseconds;
                Console.WriteLine("map " + tMap + " ms, L4+ap " + t4 + ", L3+ap " + t3 + ", L4 " + t4n + " ms  (" + sink.ToString("0") + ")");
            }
            return 0;
        }

        // ff-replica <pairs> <photo.png> <out.png> [fixOffset]: Fisheye Friend's CPU renderer at full size
        private static int ReplicaCmd(string[] a)
        {
            var g = Geo(a[1]);
            var src = PngCodec.Decode(File.ReadAllBytes(a[2]));
            var sw = Stopwatch.StartNew();
            var r = FisheyeFriendOriginal(src, new DefishWarpMap(g, src.Width, src.Height, src.Width, src.Height), src.Width, src.Height, a.Length > 4);
            Console.WriteLine("rendered " + src.Width + "x" + src.Height + " in " + sw.ElapsedMilliseconds + " ms");
            File.WriteAllBytes(a[3], PngCodec.Encode(r));
            return 0;
        }

        private static int GridCmd(string[] a)
        {
            var map = HemiMapping.Parse(File.ReadAllText(a[1]));
            foreach (int n in new[] { 129, 257, 513, 1025 })
            {
                var sw = Stopwatch.StartNew();
                var g = new DefishGeometry(map, 1e-5, n);
                long ms = sw.ElapsedMilliseconds;
                Console.WriteLine($"{n,5} nodes across ({g.GridNodes} nodes, built in {ms} ms): max error {g.MaxGridError(7952, 5304):0.00000} px on a 7952 x 5304 photo");
            }
            return 0;
        }

        private static int Analyze(string[] a)
        {
            var map = HemiMapping.Parse(File.ReadAllText(a[1]));
            int W = a.Length > 3 ? int.Parse(a[2]) : 7952, H = a.Length > 3 ? int.Parse(a[3]) : 5304;
            var sw = Stopwatch.StartNew();
            ThinPlateSpline tx, ty;
            map.FitInverse(1e-5, out tx, out ty);
            Console.WriteLine(map.Pairs.Count + " pairs, fitted in " + sw.ElapsedMilliseconds + " ms (" + map.SourceWidth + "x" + map.SourceHeight + " -> " + map.TargetWidth + "x" + map.TargetHeight + ")");
            double minSx = 1, maxSx = 0, minSy = 1, maxSy = 0, minDx = 1, maxDx = 0, minDy = 1, maxDy = 0;
            foreach (var p in map.Pairs)
            {
                minSx = Math.Min(minSx, p.SrcX); maxSx = Math.Max(maxSx, p.SrcX); minSy = Math.Min(minSy, p.SrcY); maxSy = Math.Max(maxSy, p.SrcY);
                minDx = Math.Min(minDx, p.DstX); maxDx = Math.Max(maxDx, p.DstX); minDy = Math.Min(minDy, p.DstY); maxDy = Math.Max(maxDy, p.DstY);
            }
            Console.WriteLine($"pairs cover src x {minSx:0.###}-{maxSx:0.###} y {minSy:0.###}-{maxSy:0.###};  dst x {minDx:0.###}-{maxDx:0.###} y {minDy:0.###}-{maxDy:0.###}");
            // where the output frame lands in the photo
            foreach (var pt in new[] { new[] { 0.0, 0.0 }, new[] { 0.5, 0.0 }, new[] { 1.0, 0.0 }, new[] { 0.0, 0.5 }, new[] { 0.5, 0.5 }, new[] { 1.0, 0.5 }, new[] { 0.0, 1.0 }, new[] { 0.5, 1.0 }, new[] { 1.0, 1.0 } })
                Console.WriteLine($"  output ({pt[0]:0.0},{pt[1]:0.0}) <- photo ({tx.Eval(pt[0], pt[1]):0.0000},{ty.Eval(pt[0], pt[1]):0.0000})");
            // local scale: singular values of d(src px)/d(out px), output = same size as source
            Console.WriteLine("local scale (source pixels per output pixel: min/max singular value) over the frame:");
            double gmin = 9, gmax = 0;
            for (int j = 0; j <= 4; j++)
            {
                var line = "";
                for (int i = 0; i <= 6; i++)
                {
                    double xn = i / 6.0, yn = j / 4.0, a00, a01, a10, a11;
                    tx.Eval(xn, yn, out a00, out a01);
                    ty.Eval(xn, yn, out a10, out a11);
                    // pixel units: src px = srcN*W, out px = outN*W (same size)
                    double j00 = a00, j01 = a01 * W / H, j10 = a10 * H / W, j11 = a11;
                    double t = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11, d = j00 * j11 - j01 * j10;
                    double s1 = Math.Sqrt(Math.Max(0, (t + Math.Sqrt(Math.Max(0, t * t - 4 * d * d))) / 2)), s2 = Math.Abs(d) / Math.Max(s1, 1e-12);
                    gmin = Math.Min(gmin, s2); gmax = Math.Max(gmax, s1);
                    line += $"  {s2:0.00}/{s1:0.00}";
                }
                Console.WriteLine(line);
            }
            Console.WriteLine($"overall: {gmin:0.000} .. {gmax:0.000}");
            // full-frame scan: does every output pixel map inside the photo?
            int outside = 0, n = 0;
            for (int y = 0; y <= 200; y++)
                for (int x = 0; x <= 300; x++)
                {
                    double sx = tx.Eval(x / 300.0, y / 200.0), sy = ty.Eval(x / 300.0, y / 200.0);
                    n++;
                    if (sx < 0 || sx > 1 || sy < 0 || sy > 1) outside++;
                }
            Console.WriteLine("output points mapping outside the photo: " + outside + " of " + n);
            return 0;
        }
    }
}
