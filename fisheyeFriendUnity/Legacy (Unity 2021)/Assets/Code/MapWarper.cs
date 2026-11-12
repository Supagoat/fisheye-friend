using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using UnityEngine;
using FisheyeFriend;

namespace FisheyeFriend {
    public class MapWarper {
        
        public Bitmap output;
        Renderer renderer;
        static void Main (string[] args) {
            if (args.Length < 3) {
                Console.WriteLine("Usage: MapDefish <mapping.json> <sourceImage> <output.png> [--lambda=1e-5] [--outW=W --outH=H] [--clampEdge]");
                return;
            }
            /*
            string mapPath = args[0];
            string srcPath = args[1];
            string outPath = args[2];

            double lambda = 1e-5;
            int? outWArg = null, outHArg = null;
            bool clampEdge = false;
            */
        }

        public bool IsRenderComplete () {
            return renderer.renderComplete;
        }

        public void WarpImage (WarpConfig config) {

            /* for (int i = 3; i < args.Length; i++) {
                 string s = args[i];
                 if (s.StartsWith("--lambda=")) lambda = double.Parse(s.Substring(9), CultureInfo.InvariantCulture);
                 else if (s.StartsWith("--outW=")) outWArg = int.Parse(s.Substring(7));
                 else if (s.StartsWith("--outH=")) outHArg = int.Parse(s.Substring(7));
                 else if (s == "--clampEdge") clampEdge = true;
             }*/
            throw new Exception("Deprecated: Use WarpImage(image) instead.");
            Mapping map = Mapping.Read(config.mapPath);

            using Bitmap src = new Bitmap(config.srcPath);
            int sw = src.Width, sh = src.Height;
            int outW = config.outWArg ?? sw;
            int outH = config.outHArg ?? sh;

            Debug.Log($"Fitting TPS (normalized) lambda={config.lambda}");
            TPS2D tps = TPS2D.FitInverseNormalized(map, config.lambda);

            Debug.Log($"Rendering {outW}x{outH} with radial Lanczos-3" + (config.clampEdge ? ", clampEdge" : ""));
            using Bitmap dst = new Bitmap(outW, outH, PixelFormat.Format32bppArgb);


            //Renderer.RenderLanczos3Radial(src, dst, tps, config.clampEdge);

            Renderer renderer = new Renderer(src, dst, tps, config.clampEdge);
            renderer.BeginRender();
            Debug.Log("No longer supported path to warp");
            //dst.Save(config.outPath, ImageFormat.Png);
            //Debug.Log($"Wrote {config.outPath}");
        }

        /*
        public async Task<Bitmap> WarpImageAsync(Bitmap src) {
            Bitmap fromNonAsync =  WarpImage(src);
           // using (Graphics g = Graphics.FromImage(dst)) {
           //     g.DrawImage(fromNonAsync, 0, 0, dst.Width, dst.Height);
           // }

        }*/
        
        public void WarpImage (Bitmap src) {
            // Mapping map = Mapping.Read("hemi_pairs_596.json");
            Mapping map = Mapping.Read(((TextAsset)Resources.Load("hemi_pairs_596")).text);
            int sw = src.Width, sh = src.Height;
            int outW = sw;
            int outH = sh;
            TPS2D tps = TPS2D.FitInverseNormalized(map, 1e-5);
            output = new Bitmap(outW, outH, PixelFormat.Format32bppArgb);
            renderer = new Renderer(src, output, tps, true);
            renderer.BeginRender();
        }
        public bool RenderIsComplete () {
            return renderer.CheckCompletion();
        }
    }



    public class WarpConfig {
        public string mapPath;
        public string srcPath;
        public string outPath;

        public double lambda = 1e-5; //TODO: Make configurable
        public int? outWArg = null, outHArg = null;
        public bool clampEdge = false;
    }

    public class Mapping {
        public int sW, sH, tW, tH;
        public List<Pair> Pairs = new();

        public class Pair {
            public double sxN, syN, dxN, dyN;
        }

        public static Mapping Read (string text) {
           // string text = File.ReadAllText(path);
            Mapping m = new();

            m.sW = ReadInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
            m.sH = ReadInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");
            m.tW = ReadInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
            m.tH = ReadInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");

            Regex regex = new(
                "\\{\\s*\"id\"\\s*:\\s*\\d+[^}]*?\"src\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([-0-9.]+)[^}]*?\"y\"\\s*:\\s*([-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([-0-9.]+))?[^}]*?\\}\\s*,\\s*\"dst\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([-0-9.]+)[^}]*?\"y\"\\s*:\\s*([-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([-0-9.]+))?",
                RegexOptions.Singleline);

            var matches = regex.Matches(text);
            foreach (Match match in matches) {
                double sx = double.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
                double sy = double.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture);
                string sxn = match.Groups[3].Value, syn = match.Groups[4].Value;
                double dx = double.Parse(match.Groups[5].Value, CultureInfo.InvariantCulture);
                double dy = double.Parse(match.Groups[6].Value, CultureInfo.InvariantCulture);
                string dxn = match.Groups[7].Value, dyn = match.Groups[8].Value;

                Pair p = new() {
                    sxN = sxn != "" ? double.Parse(sxn, CultureInfo.InvariantCulture) : sx / m.sW,
                    syN = syn != "" ? double.Parse(syn, CultureInfo.InvariantCulture) : sy / m.sH,
                    dxN = dxn != "" ? double.Parse(dxn, CultureInfo.InvariantCulture) : dx / m.tW,
                    dyN = dyn != "" ? double.Parse(dyn, CultureInfo.InvariantCulture) : dy / m.tH
                };
                m.Pairs.Add(p);
            }

            if (m.Pairs.Count == 0) throw new IOException("No pairs found in mapping");
            Console.WriteLine($"Loaded {m.Pairs.Count} pairs (normalized).");
            return m;
        }

        private static int ReadInt (string text, string pattern) {
            Match m = Regex.Match(text, pattern);
            if (!m.Success) throw new IOException($"Missing field: {pattern}");
            return int.Parse(m.Groups[1].Value);
        }
    }

    public class TPS2D {
        private readonly ThinPlateSpline tx, ty;
        public TPS2D (ThinPlateSpline tx, ThinPlateSpline ty) { this.tx = tx; this.ty = ty; }

        public static TPS2D FitInverseNormalized (Mapping m, double lambda) {
            int n = m.Pairs.Count;
            double[] xd = new double[n], yd = new double[n], xs = new double[n], ys = new double[n];

            for (int i = 0; i < n; i++) {
                var p = m.Pairs[i];
                xd[i] = p.dxN; yd[i] = p.dyN; xs[i] = p.sxN; ys[i] = p.syN;
            }

            return new TPS2D(
                ThinPlateSpline.Fit(xd, yd, xs, lambda),
                ThinPlateSpline.Fit(xd, yd, ys, lambda)
            );
        }

        public void MapNorm (double xn, double yn, double[] out2) {
            out2[0] = tx.Eval(xn, yn);
            out2[1] = ty.Eval(xn, yn);
        }
        public TPSCoeffs ExportCoeffs (int sw, int sh, int dw, int dh) {
            return new TPSCoeffs {
                a0x = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).A0,
                axx = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).AX,
                ayx = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).AY,

                a0y = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("ty", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).A0,
                axy = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("ty", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).AX,
                ayy = (float)((ThinPlateSpline)typeof(TPS2D)
                    .GetField("ty", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).AY,
                dims = new Vector4(sw, sh, dw, dh),
                count = ((ThinPlateSpline)typeof(TPS2D)
                    .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                    .GetValue(this)).Count
            };
        }

        // Also add strongly-typed accessors for the control arrays:
        public (float[] x, float[] y, float[] w) ExportTxArrays () {
            var tx = (ThinPlateSpline)typeof(TPS2D)
                .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .GetValue(this);
            int n = tx.Count;
            var xx = new float[n]; var yy = new float[n]; var ww = new float[n];
            for (int i = 0; i < n; i++) { xx[i] = (float)tx.X[i]; yy[i] = (float)tx.Y[i]; ww[i] = (float)tx.W[i]; }
            return (xx, yy, ww);
        }
        public (float[] x, float[] y, float[] w) ExportTyArrays () {
            var ty = (ThinPlateSpline)typeof(TPS2D)
                .GetField("ty", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .GetValue(this);
            int n = ty.Count;
            var xx = new float[n]; var yy = new float[n]; var ww = new float[n];
            for (int i = 0; i < n; i++) { xx[i] = (float)ty.X[i]; yy[i] = (float)ty.Y[i]; ww[i] = (float)ty.W[i]; }
            return (xx, yy, ww);
        }
    }

    /*
    public class TPSActor {
       

        private ManualResetEvent Signal { get; set; }


        public MobMoveActor ( ManualResetEvent signal) {
            this.Signal = signal;
        }

        public void ThreadPoolCallback (System.Object threadContext) {
            try {

                Engine.MoveMobs(Mobs, NewMobs, ToDestroyMobs);
            } catch (Exception e) {
                Debug.Log(e);
            } finally {
                Signal.Set();
            }
           
        }

    }
    */

    


    public struct TPSCoeffs {
        public float a0x, axx, ayx;  // for tx
        public float a0y, axy, ayy;  // for ty
        public Vector4 dims;         // sw, sh, dw, dh
        public int count;
    }

    public class ThinPlateSpline {
        private readonly double[] x, y, w;
        private readonly double a0, ax, ay;
        public int Count => x.Length;
        public IReadOnlyList<double> X => x;
        public IReadOnlyList<double> Y => y;
        public IReadOnlyList<double> W => w;
        public double A0 => a0;
        public double AX => ax;
        public double AY => ay;

        private ThinPlateSpline (double[] x, double[] y, double[] w, double a0, double ax, double ay) {
            this.x = x; this.y = y; this.w = w; this.a0 = a0; this.ax = ax; this.ay = ay;
        }

        public static ThinPlateSpline Fit (double[] xi, double[] yi, double[] zi, double lambda) {
            DateTime t = DateTime.Now;
            int n = xi.Length, dim = n + 3;
            double[,] A = new double[dim, dim];
            double[] b = new double[dim];

            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++) {
                    double r = Math.Sqrt((xi[i] - xi[j]) * (xi[i] - xi[j]) + (yi[i] - yi[j]) * (yi[i] - yi[j]));
                    double v = U(r) + (i == j ? lambda : 0);
                    A[i, j] = v;
                    A[j, i] = v;
                }

            for (int i = 0; i < n; i++) {
                A[i, n] = 1; A[n, i] = 1;
                A[i, n + 1] = xi[i]; A[n + 1, i] = xi[i];
                A[i, n + 2] = yi[i]; A[n + 2, i] = yi[i];
                b[i] = zi[i];
            }

            double[] sol = Solve(A, b);
            double[] w = new double[n];
            Array.Copy(sol, w, n);
          //  Debug.Log("TPS fit time: " + (DateTime.Now - t) + " ms");
            return new ThinPlateSpline((double[])xi.Clone(), (double[])yi.Clone(), w, sol[n], sol[n + 1], sol[n + 2]);
        }

        public double Eval (double px, double py) {
            double s = a0 + ax * px + ay * py;
            for (int i = 0; i < x.Length; i++) {
                double dx = px - x[i], dy = py - y[i];
                double r = Math.Sqrt(dx * dx + dy * dy);
                s += w[i] * U(r);
            }
            return s;
        }

        private static double U (double r) => r <= 1e-12 ? 0 : r * r * Math.Log(r);

        private static double[] Solve (double[,] A, double[] b) {
            int n = b.Length;
            double[,] M = (double[,])A.Clone();
            double[] x = new double[n];
            double[] rhs = (double[])b.Clone();

            for (int k = 0; k < n; k++) {
                int piv = k;
                double mx = Math.Abs(M[k, k]);
                for (int i = k + 1; i < n; i++) {
                    double v = Math.Abs(M[i, k]);
                    if (v > mx) { mx = v; piv = i; }
                }
                if (mx < 1e-14) throw new Exception("Singular TPS system");
                if (piv != k) {
                    for (int j = 0; j < n; j++) {
                        (M[k, j], M[piv, j]) = (M[piv, j], M[k, j]);
                    }
                    (rhs[k], rhs[piv]) = (rhs[piv], rhs[k]);
                }
                double akk = M[k, k];
                for (int i = k + 1; i < n; i++) {
                    double f = M[i, k] / akk;
                    rhs[i] -= f * rhs[k];
                    for (int j = k; j < n; j++) M[i, j] -= f * M[k, j];
                }
            }

            for (int i = n - 1; i >= 0; i--) {
                double s = rhs[i];
                for (int j = i + 1; j < n; j++) s -= M[i, j] * x[j];
                x[i] = s / M[i, i];
            }

            return x;
        }
    }

    class RenderActor {
        public int startX;
        public int startY;
        public int endX;
        public int endY;
        public BitmapData srcData;
        public BitmapData dstData;
        public ManualResetEvent reset;
        public TPS2D tps;
        bool clampEdge;

        public RenderActor (int startX, int startY, int endX, int endY, TPS2D tps, bool clampEdge, BitmapData src, BitmapData dst, ManualResetEvent reset) {
            this.startX = startX;
            this.startY = startY;
            this.endX = endX;
            this.endY = endY;
            this.srcData = src;
            this.dstData = dst;
            this.reset = reset;
            this.tps = tps;
            this.clampEdge = clampEdge; 
        }

        public void ThreadPoolCallback (System.Object threadContext) {
            try {
                RenderAPortion();
            } catch (Exception e) {
                Debug.Log(e);
            } finally {
                Debug.Log("Thread complete");
                reset.Set();
            }

        }

        private static double Sinc (double t) => Math.Abs(t) < 1e-8 ? 1.0 : Math.Sin(Math.PI * t) / (Math.PI * t);
        private static double Lanczos (double r, double a) => Sinc(r) * Sinc(r / a);


        public void RenderWhole () {
            int sw = srcData.Width, sh = srcData.Height;
            int dw = dstData.Width, dh = dstData.Height;
     
            long t = DateTime.Now.Millisecond;
            Debug.Log("Beginning render");


            unsafe {
                int* sdata = (int*)srcData.Scan0;
                int* ddata = (int*)dstData.Scan0;
                double a = 3.0;
                double[] map = new double[2];

                for (int y = 0; y < dh; y++) {
                    int row = y * dw;
                    double yn = (y + 0.5) / dh;
                    for (int x = 0; x < dw; x++) {
                        double xn = (x + 0.5) / dw;
                        tps.MapNorm(xn, yn, map);
                        double sx = map[0] * sw - 0.5;
                        double sy = map[1] * sh - 0.5;

                        int xmin = (int)Math.Floor(sx - a + 1);
                        int xmax = (int)Math.Ceiling(sx + a);
                        int ymin = (int)Math.Floor(sy - a + 1);
                        int ymax = (int)Math.Ceiling(sy + a);

                        double W = 0, Aa = 0, Ar = 0, Ag = 0, Ab = 0;

                        for (int yy = ymin; yy <= ymax; yy++) {
                            for (int xx = xmin; xx <= xmax; xx++) {
                                double dx = (xx + 0.5) - sx;
                                double dy = (yy + 0.5) - sy;
                                double r = Math.Sqrt(dx * dx + dy * dy);
                                if (r >= a) continue;
                                double w = Lanczos(r, a);
                                int px = xx, py = yy;
                                if (clampEdge) {
                                    px = Math.Clamp(px, 0, sw - 1);
                                    py = Math.Clamp(py, 0, sh - 1);
                                } else {
                                    if (px < 0 || py < 0 || px >= sw || py >= sh) continue;
                                }

                                int c = sdata[py * sw + px];
                                double a0 = ((c >> 24) & 255) / 255.0;
                                double r0 = ((c >> 16) & 255) / 255.0;
                                double g0 = ((c >> 8) & 255) / 255.0;
                                double b0 = (c & 255) / 255.0;

                                W += w;
                                Aa += w * a0;
                                Ar += w * r0 * a0;
                                Ag += w * g0 * a0;
                                Ab += w * b0 * a0;
                            }
                        }

                        int outARGB;
                        if (W <= 1e-12 || Aa <= 1e-12) {
                            outARGB = 0;
                        } else {
                            double A1 = Aa / W;
                            double R1 = (Ar / W) / (A1 + 1e-12);
                            double G1 = (Ag / W) / (A1 + 1e-12);
                            double B1 = (Ab / W) / (A1 + 1e-12);
                            int ia = (int)Math.Clamp(Math.Round(A1 * 255), 0, 255);
                            int ir = (int)Math.Clamp(Math.Round(R1 * 255), 0, 255);
                            int ig = (int)Math.Clamp(Math.Round(G1 * 255), 0, 255);
                            int ib = (int)Math.Clamp(Math.Round(B1 * 255), 0, 255);
                            outARGB = (ia << 24) | (ir << 16) | (ig << 8) | ib;
                        }

                        ddata[row + x] = outARGB;
                    }
                }
            }
            Debug.Log("lancLoop time: " + (DateTime.Now.Millisecond - t) + " ms");
            //src.UnlockBits(srcData);
           // dst.UnlockBits(dstData);
        }

        public void RenderAPortion () {
            int sw = srcData.Width, sh = srcData.Height;
            int dw = dstData.Width, dh = dstData.Height;

            long t = DateTime.Now.Millisecond;


            unsafe {
                int* sdata = (int*)srcData.Scan0;
                int* ddata = (int*)dstData.Scan0;
                double a = 3.0;
                double[] map = new double[2];

                for (int y = startY; y < endY; y++) {
                    int row = y * dw;
                    double yn = (y + 0.5) / dh;
                    for (int x = startX; x < endX; x++) {
                        double xn = (x + 0.5) / dw;
                        tps.MapNorm(xn, yn, map);
                        double sx = map[0] * sw - 0.5;
                        double sy = map[1] * sh - 0.5;

                        int xmin = (int)Math.Floor(sx - a + 1);
                        int xmax = (int)Math.Ceiling(sx + a);
                        int ymin = (int)Math.Floor(sy - a + 1);
                        int ymax = (int)Math.Ceiling(sy + a);

                        double W = 0, Aa = 0, Ar = 0, Ag = 0, Ab = 0;

                        for (int yy = ymin; yy <= ymax; yy++) {
                            for (int xx = xmin; xx <= xmax; xx++) {
                                double dx = (xx + 0.5) - sx;
                                double dy = (yy + 0.5) - sy;
                                double r = Math.Sqrt(dx * dx + dy * dy);
                                if (r >= a) continue;
                                double w = Lanczos(r, a);
                                int px = xx, py = yy;
                                if (clampEdge) {
                                    px = Math.Clamp(px, 0, sw - 1);
                                    py = Math.Clamp(py, 0, sh - 1);
                                } else {
                                    if (px < 0 || py < 0 || px >= sw || py >= sh) continue;
                                }

                                int c = sdata[py * sw + px];
                                double a0 = ((c >> 24) & 255) / 255.0;
                                double r0 = ((c >> 16) & 255) / 255.0;
                                double g0 = ((c >> 8) & 255) / 255.0;
                                double b0 = (c & 255) / 255.0;

                                W += w;
                                Aa += w * a0;
                                Ar += w * r0 * a0;
                                Ag += w * g0 * a0;
                                Ab += w * b0 * a0;
                            }
                        }

                        int outARGB;
                        if (W <= 1e-12 || Aa <= 1e-12) {
                            outARGB = 0;
                        } else {
                            double A1 = Aa / W;
                            double R1 = (Ar / W) / (A1 + 1e-12);
                            double G1 = (Ag / W) / (A1 + 1e-12);
                            double B1 = (Ab / W) / (A1 + 1e-12);
                            int ia = (int)Math.Clamp(Math.Round(A1 * 255), 0, 255);
                            int ir = (int)Math.Clamp(Math.Round(R1 * 255), 0, 255);
                            int ig = (int)Math.Clamp(Math.Round(G1 * 255), 0, 255);
                            int ib = (int)Math.Clamp(Math.Round(B1 * 255), 0, 255);
                            outARGB = (ia << 24) | (ir << 16) | (ig << 8) | ib;
                        }

                        ddata[row + x] = outARGB;
                    }
                }
            }


        }


        public void RenderPortion () {
            int sw = srcData.Width, sh = srcData.Height;
            int dw = dstData.Width, dh = dstData.Height;

            unsafe {
                int* sdata = (int*)srcData.Scan0;
                int* ddata = (int*)dstData.Scan0;
                double a = 3.0;
                double[] map = new double[2];

              //  Debug.Log($"Thread {startX},{startY} starting render portion {startX}-{endX} x {startY}-{endY}");
 
                for (int y = startY; y < endY; y++) {
                    int row = y * endX;
                    double yn = (y + 0.5) / endY;
                    for (int x = startX; x < endX; x++) {
                        double xn = (x + 0.5) / endX;
                        tps.MapNorm(xn, yn, map);
                        double sx = map[0] * endX - 0.5;
                        double sy = map[1] * endY - 0.5;

                        int xmin = (int)Math.Floor(sx - a + 1);
                        int xmax = (int)Math.Ceiling(sx + a);
                        int ymin = (int)Math.Floor(sy - a + 1);
                        int ymax = (int)Math.Ceiling(sy + a);

                        double W = 0, Aa = 0, Ar = 0, Ag = 0, Ab = 0;

                        for (int yy = ymin; yy <= ymax; yy++) {
                            //Debug.Log($"Thread {startX},{startY} processing pixel {x},{y} sample line {yy}");
                            for (int xx = xmin; xx <= xmax; xx++) {
                                double dx = (xx + 0.5) - sx;
                                double dy = (yy + 0.5) - sy;
                                double r = Math.Sqrt(dx * dx + dy * dy);
                                if (r >= a) continue;
                                double w = Lanczos(r, a);
                                int px = xx, py = yy;
                                if (clampEdge) {
                                    px = Math.Clamp(px, 0, endX - 1);
                                    py = Math.Clamp(py, 0, endY - 1);
                                } else {
                                    if (px < 0 || py < 0 || px >= endX || py >= endY) continue;
                                }

                                int c = sdata[py * endX + px];
                                double a0 = ((c >> 24) & 255) / 255.0;
                                double r0 = ((c >> 16) & 255) / 255.0;
                                double g0 = ((c >> 8) & 255) / 255.0;
                                double b0 = (c & 255) / 255.0;

                                W += w;
                                Aa += w * a0;
                                Ar += w * r0 * a0;
                                Ag += w * g0 * a0;
                                Ab += w * b0 * a0;
                            }
                        }

                        int outARGB;
                        if (W <= 1e-12 || Aa <= 1e-12) {
                            outARGB = 0;
                        } else {
                            double A1 = Aa / W;
                            double R1 = (Ar / W) / (A1 + 1e-12);
                            double G1 = (Ag / W) / (A1 + 1e-12);
                            double B1 = (Ab / W) / (A1 + 1e-12);
                            int ia = (int)Math.Clamp(Math.Round(A1 * 255), 0, 255);
                            int ir = (int)Math.Clamp(Math.Round(R1 * 255), 0, 255);
                            int ig = (int)Math.Clamp(Math.Round(G1 * 255), 0, 255);
                            int ib = (int)Math.Clamp(Math.Round(B1 * 255), 0, 255);
                            outARGB = (ia << 24) | (ir << 16) | (ig << 8) | ib;
                        }

                        ddata[row + x] = outARGB;
                    }
                }
            }

        }

    }
}
class Renderer {
    private Bitmap src;
    private Bitmap dst;
    private BitmapData srcData;
    private BitmapData dstData;
    private TPS2D tps;
    private bool clampEdge;
    ManualResetEvent[] waitOn;
    public bool renderComplete = false;

    public Renderer (Bitmap src, Bitmap dst, TPS2D tps, bool clampEdge) {
        this.src = src;
        this.dst = dst;
        this.tps = tps;
        this.clampEdge = clampEdge;
    }

    //TODO: Split into as many portions as there is CPU.  Consider doing vertical or horizontal slices on each thread
    public void BeginRender () {

        if(waitOn != null) {
            Debug.Log("Attempting a render while one is in progress");
            return;
        }
        renderComplete = false;
        int sw = src.Width, sh = src.Height;
        int dw = dst.Width, dh = dst.Height;
        srcData = src.LockBits(new Rectangle(0, 0, sw, sh), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);
        dstData = dst.LockBits(new Rectangle(0, 0, dw, dh), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
       
        waitOn = new ManualResetEvent []{ new ManualResetEvent(false), new ManualResetEvent(false), new ManualResetEvent(false), new ManualResetEvent(false) };
        RenderActor[] actors = {
                new RenderActor(0, 0, sw/2, sh/2, tps, clampEdge, srcData, dstData, waitOn[0]),
                new RenderActor(0, sh/2, sw/2, sh, tps, clampEdge, srcData, dstData, waitOn[1]),
                new RenderActor(sw/2, 0, sw, sh/2, tps, clampEdge, srcData, dstData, waitOn[2]),
                new RenderActor(sw/2, sh/2, sw, sh, tps, clampEdge, srcData, dstData, waitOn[3]) };


          for (int i = 0; i < actors.Length; i++) {
                ThreadPool.QueueUserWorkItem(actors[i].ThreadPoolCallback);
        }
      //  WaitHandle.WaitAll(waitOn);


    }

    public bool CheckCompletion() {
        foreach(ManualResetEvent threadCheck in waitOn) {
            if(!threadCheck.WaitOne(0)) {
                return false;
            }
        }
        OnRenderComplete();
        return true;
    }

    void OnRenderComplete() {
        src.UnlockBits(srcData);
        dst.UnlockBits(dstData);
        waitOn = null;
        renderComplete = true;
    }
    /*


    public static void RenderLanczos3Radial (Bitmap src, Bitmap dst, TPS2D tps, bool clampEdge) {
        int sw = src.Width, sh = src.Height;
        int dw = dst.Width, dh = dst.Height;
        BitmapData srcData = src.LockBits(new Rectangle(0, 0, sw, sh), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);
        BitmapData dstData = dst.LockBits(new Rectangle(0, 0, dw, dh), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
        long t = DateTime.Now.Millisecond;
        Debug.Log("Beginning render");


        unsafe {
            int* sdata = (int*)srcData.Scan0;
            int* ddata = (int*)dstData.Scan0;
            double a = 3.0;
            double[] map = new double[2];

            for (int y = 0; y < dh; y++) {
                int row = y * dw;
                double yn = (y + 0.5) / dh;
                for (int x = 0; x < dw; x++) {
                    double xn = (x + 0.5) / dw;
                    tps.MapNorm(xn, yn, map);
                    double sx = map[0] * sw - 0.5;
                    double sy = map[1] * sh - 0.5;

                    int xmin = (int)Math.Floor(sx - a + 1);
                    int xmax = (int)Math.Ceiling(sx + a);
                    int ymin = (int)Math.Floor(sy - a + 1);
                    int ymax = (int)Math.Ceiling(sy + a);

                    double W = 0, Aa = 0, Ar = 0, Ag = 0, Ab = 0;

                    for (int yy = ymin; yy <= ymax; yy++) {
                        for (int xx = xmin; xx <= xmax; xx++) {
                            double dx = (xx + 0.5) - sx;
                            double dy = (yy + 0.5) - sy;
                            double r = Math.Sqrt(dx * dx + dy * dy);
                            if (r >= a) continue;
                            double w = Lanczos(r, a);
                            int px = xx, py = yy;
                            if (clampEdge) {
                                px = Math.Clamp(px, 0, sw - 1);
                                py = Math.Clamp(py, 0, sh - 1);
                            } else {
                                if (px < 0 || py < 0 || px >= sw || py >= sh) continue;
                            }

                            int c = sdata[py * sw + px];
                            double a0 = ((c >> 24) & 255) / 255.0;
                            double r0 = ((c >> 16) & 255) / 255.0;
                            double g0 = ((c >> 8) & 255) / 255.0;
                            double b0 = (c & 255) / 255.0;

                            W += w;
                            Aa += w * a0;
                            Ar += w * r0 * a0;
                            Ag += w * g0 * a0;
                            Ab += w * b0 * a0;
                        }
                    }

                    int outARGB;
                    if (W <= 1e-12 || Aa <= 1e-12) {
                        outARGB = 0;
                    } else {
                        double A1 = Aa / W;
                        double R1 = (Ar / W) / (A1 + 1e-12);
                        double G1 = (Ag / W) / (A1 + 1e-12);
                        double B1 = (Ab / W) / (A1 + 1e-12);
                        int ia = (int)Math.Clamp(Math.Round(A1 * 255), 0, 255);
                        int ir = (int)Math.Clamp(Math.Round(R1 * 255), 0, 255);
                        int ig = (int)Math.Clamp(Math.Round(G1 * 255), 0, 255);
                        int ib = (int)Math.Clamp(Math.Round(B1 * 255), 0, 255);
                        outARGB = (ia << 24) | (ir << 16) | (ig << 8) | ib;
                    }

                    ddata[row + x] = outARGB;
                }
            }
        }
        Debug.Log("lancLoop time: " + (DateTime.Now.Millisecond - t) + " ms");
        src.UnlockBits(srcData);
        dst.UnlockBits(dstData);
    }*/
}

