using System;
using System.Threading.Tasks;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Gpu;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriendBench
{
    /// <summary>
    /// FisheyeFriendDefish.compute translated line for line into C#, computing in 32-bit floats as the GPU does (every
    /// literal is a float: C# would otherwise silently promote to double). It runs the shader's kernels on the inputs
    /// GpuPlan prepares, so the GPU design (and those inputs) can be checked against the CPU defish without a GPU. The
    /// real GPU's sin, cos, log and atan are approximations of their own, so it can differ from this in the last bits.
    /// </summary>
    internal sealed class GpuEmulator
    {
        // ======================================================================== MapModel

        private readonly GpuPlan.ModelParams _m;

        public GpuEmulator(GpuPlan.ModelParams m)
        {
            _m = m;
        }

        private float LensAngle(float r)
        {
            float r2 = r * r, p = r, s = 0f;
            for (int i = 0; i < 8; i++)
            {
                s += _m.Theta[i] * p;
                p *= r2;
            }
            return s;
        }

        private float LensAngleDerivative(float r)
        {
            float r2 = r * r, p = 1f, s = 0f;
            for (int i = 0; i < 8; i++)
            {
                s += (2 * i + 1) * _m.Theta[i] * p;
                p *= r2;
            }
            return s;
        }

        private float LensRadius(float theta)
        {
            float t0 = _m.Theta[0];
            float r = theta / t0;
            for (int it = 0; it < 50; it++)
            {
                float f = LensAngle(r) - theta, d = LensAngleDerivative(r);
                if (d <= 1e-9f) d = t0;
                float step = Math.Clamp(f / d, -0.2f, 0.2f);
                r -= step;
                if (r < 0f) r = 0f;
                if (MathF.Abs(step) < 1e-7f) break;
            }
            return r;
        }

        private float ProjectX(float lon, out float derivative)
        {
            float s = MathF.Sin(lon), c = MathF.Cos(lon), u = s * s;
            float den = _m.Panini + c;
            float basis = (_m.Panini + 1f) * s / den;
            float dBasis = (_m.Panini + 1f) * (c * den + s * s) / (den * den);
            float corr = 1f, dCorr = 0f, prev = 1f;
            for (int i = 0; i < 4; i++)
            {
                float h = _m.Horizontal[i];
                corr += h * prev * u;
                dCorr += h * (i + 1) * prev * 2f * s * c;
                prev *= u;
            }
            derivative = dBasis * corr + basis * dCorr;
            return basis * corr;
        }

        private float ProjectY(float lon, float lat, out float dLat)
        {
            float sl = MathF.Sin(lon), cl = MathF.Cos(lon), u = sl * sl;
            float sp = MathF.Sin(lat), cp = MathF.Cos(lat), v = sp * sp;
            float side = (_m.PaniniVertical + 1f) / (_m.PaniniVertical + cl);
            float tn = sp / cp, dTan = 1f / (cp * cp);
            float corr = 1f, dCorr = 0f, prev = 1f;
            for (int i = 0; i < 4; i++)
            {
                float row = 0f, up = 1f;
                for (int j = 0; j < 4; j++)
                {
                    row += _m.Vertical[i * 4 + j] * up;
                    up *= u;
                }
                corr += row * prev * v;
                dCorr += row * (i + 1) * prev * 2f * sp * cp;
                prev *= v;
            }
            dLat = side * (dTan * corr + tn * dCorr);
            return side * tn * corr;
        }

        private static float Tangent(float[] y, int k, float h)
        {
            if (k == 0) return 0f;
            if (k == 3) return (y[3] - y[2]) / h;
            float d0 = (y[k] - y[k - 1]) / h, d1 = (y[k + 1] - y[k]) / h;
            if (d0 * d1 <= 0f) return 0f;
            return 2f / (1f / d0 + 1f / d1);
        }

        private static float ShapeAt(float[] y, float rho)
        {
            float h = 1f / 3f;
            if (rho <= 0f) return y[0];
            if (rho >= 1f) return y[3] + Tangent(y, 3, h) * (rho - 1f);
            int k = Math.Min(2, (int)(rho / h));
            float t = (rho - k * h) / h;
            float m0 = Tangent(y, k, h), m1 = Tangent(y, k + 1, h);
            float t2 = t * t, t3 = t2 * t;
            return (2f * t3 - 3f * t2 + 1f) * y[k] + (t3 - 2f * t2 + t) * h * m0 + (-2f * t3 + 3f * t2) * y[k + 1] + (t3 - t2) * h * m1;
        }

        public void MapPoint(float u, float v, out float s, out float t)
        {
            float x = u - 0.5f, y = (v - 0.5f) * _m.Aspect;
            float lon = MathF.Atan(x / _m.K);
            for (int it = 0; it < 60; it++)
            {
                float d;
                float f = _m.K * ProjectX(lon, out d) - x;
                d *= _m.K;
                if (!(d > 1e-9f)) d = _m.K;
                float step = Math.Clamp(f / d, -0.2f, 0.2f);
                lon = Math.Clamp(lon - step, -1.55f, 1.55f);
                if (MathF.Abs(step) < 1e-7f) break;
            }
            float lat = MathF.Atan(y / (_m.K * 1.2f));
            for (int it2 = 0; it2 < 60; it2++)
            {
                float d;
                float f = _m.K * ProjectY(lon, lat, out d) - y;
                d *= _m.K;
                if (!(d > 1e-9f)) d = _m.K;
                float step = Math.Clamp(f / d, -0.2f, 0.2f);
                lat = Math.Clamp(lat - step, -1.5f, 1.5f);
                if (MathF.Abs(step) < 1e-7f) break;
            }
            float cp = MathF.Cos(lat);
            float X = cp * MathF.Sin(lon), Y = MathF.Sin(lat), Z = cp * MathF.Cos(lon);
            float sinTh = MathF.Sqrt(X * X + Y * Y);
            float theta = MathF.Atan2(sinTh, Z);
            float sx, sy;
            if (theta < 1e-9f)
            {
                float inv = 1f / (_m.Theta[0] * _m.Strength);
                sx = X * inv;
                sy = Y * inv;
            }
            else
            {
                float r = LensRadius(theta) / _m.Strength;
                sx = r * X / sinTh;
                sy = r * Y / sinTh;
            }
            float ax = _m.AmountX, ay = _m.AmountY;
            if (!_m.Flat)
            {
                float rho = MathF.Sqrt(x * x + y * y) * _m.RhoNorm;
                ax *= ShapeAt(_m.ShapeX, rho);
                ay *= ShapeAt(_m.ShapeY, rho);
            }
            sx = x + ax * (sx - x);
            sy = y + ay * (sy - y);
            s = sx + 0.5f;
            t = sy / _m.Aspect + 0.5f;
        }

        /// <summary>The MapModel kernel over a list of (u, v) points.</summary>
        public float[] MapPoints(float[] points)
        {
            var r = new float[points.Length];
            Parallel.For(0, points.Length / 2, i =>
            {
                float s, t;
                MapPoint(points[2 * i], points[2 * i + 1], out s, out t);
                r[2 * i] = s;
                r[2 * i + 1] = t;
            });
            return r;
        }

        // ======================================================================== Warp

        private const int MaxTaps = 32;

        private sealed class WarpState
        {
            public GpuPlan.WarpGrid G;
            public GpuPlan.WarpParams P;
            public uint[] Src;
            public int SrcW, SrcH, OutW, OutH;
        }

        private static uint SampleAt(WarpState w, int index)
        {
            if (w.P.Depth == 8) return (w.Src[index >> 2] >> ((index & 3) * 8)) & 0xFFu;
            return (w.Src[index >> 1] >> ((index & 1) * 16)) & 0xFFFFu;
        }

        private static float Lookup(WarpState w, float t)
        {
            var lut = w.P.WeightLut;
            int n = w.P.WeightLutSize;
            int i = (int)t;
            if (i >= n - 1) return lut[n - 1];
            float f = t - i;
            float a = lut[i];
            return a + (lut[i + 1] - a) * f;
        }

        private static float Weight1D(WarpState w, float t, float inv, float ap)
        {
            float R = w.P.R, S = w.P.LutScale;
            float u = MathF.Abs(t * inv);
            float wt = u >= R ? 0f : Lookup(w, u * S);
            if (ap == 0f) return wt;
            float um = MathF.Abs((t - 1f) * inv), up = MathF.Abs((t + 1f) * inv);
            float wm = um >= R ? 0f : Lookup(w, um * S);
            float wp = up >= R ? 0f : Lookup(w, up * S);
            return (1f + 2f * ap) * wt - ap * (wm + wp);
        }

        private static float Widen(float sq, float k)
        {
            if (sq <= 1f) return 1f;
            float f = 1f + k * (MathF.Sqrt(sq) - 1f);
            return f * f;
        }

        private static float ApertureCoefficient(WarpState w, float scale)
        {
            var p = w.P;
            if (!(scale > 0f)) return 0f;
            float t = (MathF.Log(scale) - p.ApLogMin) / (p.ApLogMax - p.ApLogMin) * (p.ApCount - 1);
            if (t <= 0f) return p.ApTable[0];
            if (t >= p.ApCount - 1) return p.ApTable[p.ApCount - 1];
            int i = (int)t;
            return p.ApTable[i] + (p.ApTable[i + 1] - p.ApTable[i]) * (t - i);
        }

        private static float FromSigmoid(WarpState w, float v)
        {
            var t = w.P.FromSig;
            if (!(v > 0f)) return t[0];
            if (v >= 1f) return t[65535];
            float p = v * 65535f;
            int i = (int)p;
            float f = p - i;
            float a = t[i];
            return a + (t[i + 1] - a) * f;
        }

        private static float SrgbEncode(WarpState w, float l)
        {
            var t = w.P.Encode;
            if (!(l > 0f)) return 0f;
            if (l >= 1f) return 1f;
            float p = l * 65535f;
            int i = (int)p;
            float f = p - i;
            float a = t[i];
            return a + (t[i + 1] - a) * f;
        }

        private static float Clamp01(float v)
        {
            if (!(v > 0f)) return 0f;
            return v > 1f ? 1f : v;
        }

        private static void CatmullRom(float x, float[] w, float[] d)
        {
            float x2 = x * x, x3 = x2 * x;
            w[0] = 0.5f * (-x3 + 2f * x2 - x);
            w[1] = 0.5f * (3f * x3 - 5f * x2 + 2f);
            w[2] = 0.5f * (-3f * x3 + 4f * x2 + x);
            w[3] = 0.5f * (x3 - x2);
            d[0] = 0.5f * (-3f * x2 + 4f * x - 1f);
            d[1] = 0.5f * (9f * x2 - 10f * x);
            d[2] = 0.5f * (-9f * x2 + 8f * x + 1f);
            d[3] = 0.5f * (3f * x2 - 2f * x);
        }

        private static void MapPixel(WarpState w, int x, int y, out int ix, out float fx, out int iy, out float fy,
            out float j00, out float j01, out float j10, out float j11)
        {
            var g = w.G;
            var cu = g.Cols[x];
            var cv = g.Rows[y];
            float[] wa = new float[4], da = new float[4], wb = new float[4], db = new float[4];
            CatmullRom(cu.A, wa, da);
            CatmullRom(cv.A, wb, db);
            var b = g.Nodes[cv.I * g.NodesU + cu.I];
            float sx0 = 0f, sx1 = 0f, sx2 = 0f, sy0 = 0f, sy1 = 0f, sy2 = 0f;
            for (int r = 0; r < 4; r++)
            {
                int o = (cv.I - 1 + r) * g.NodesU + cu.I - 1;
                float rs = 0f, rds = 0f, rt = 0f, rdt = 0f;
                for (int c = 0; c < 4; c++)
                {
                    var n = g.Nodes[o + c];
                    float dx = (float)(n.IX - b.IX) + (n.FX - b.FX);
                    float dy = (float)(n.IY - b.IY) + (n.FY - b.FY);
                    float ww = wa[c], dw = da[c];
                    rs += ww * dx;
                    rds += dw * dx;
                    rt += ww * dy;
                    rdt += dw * dy;
                }
                float wr = wb[r], dr = db[r];
                sx0 += wr * rs;
                sx1 += wr * rds;
                sx2 += dr * rs;
                sy0 += wr * rt;
                sy1 += wr * rdt;
                sy2 += dr * rt;
            }
            float f = b.FX + sx0;
            float fl = MathF.Floor(f);
            ix = b.IX + (int)fl;
            fx = f - fl;
            f = b.FY + sy0;
            fl = MathF.Floor(f);
            iy = b.IY + (int)fl;
            fy = f - fl;
            j00 = sx1 * g.ColStep;
            j01 = sx2 * g.RowStep;
            j10 = sy1 * g.ColStep;
            j11 = sy2 * g.RowStep;
        }

        private static void LoadSamples(WarpState w, int si, int count, float[] v)
        {
            var lut = w.P.ColorLut;
            v[0] = lut[SampleAt(w, si)];
            v[1] = count > 1 ? lut[SampleAt(w, si + 1)] : 0f;
            v[2] = count > 2 ? lut[SampleAt(w, si + 2)] : 0f;
            v[3] = count > 3 ? lut[SampleAt(w, si + 3)] : 0f;
        }

        private static uint Quantize(WarpState w, float v)
        {
            return (uint)(v * w.P.OutMax + 0.5f);
        }

        /// <summary>The Warp kernel for one output pixel; writes its word(s) into <paramref name="outWords"/> at pixel index <paramref name="o"/>.</summary>
        private static void WarpPixel(WarpState w, int x, int y, uint[] outWords, int o)
        {
            var P = w.P;
            int ch = P.Channels;
            bool alpha = P.Alpha;
            int cc = alpha ? ch - 1 : ch;

            int ix, iy;
            float fx, fy, j00, j01, j10, j11;
            MapPixel(w, x, y, out ix, out fx, out iy, out fy, out j00, out j01, out j10, out j11);

            float sxc = (float)ix + fx, syc = (float)iy + fy;
            float nx = MathF.Sqrt(j00 * j00 + j01 * j01), ny = MathF.Sqrt(j10 * j10 + j11 * j11);
            float inX = MathF.Min(sxc + 0.5f, w.SrcW - 0.5f - sxc) / MathF.Max(nx, 1e-9f);
            float inY = MathF.Min(syc + 0.5f, w.SrcH - 0.5f - syc) / MathF.Max(ny, 1e-9f);
            float cover = P.Extend ? 1f : MathF.Min(1f, 0.5f + MathF.Min(inX, inY));

            var val = new float[4];
            if (cover <= 0f)
            {
                for (int c = 0; c < 4; c++) val[c] = alpha ? 0f : P.Bg[c];
            }
            else
            {
                float a = j00 * j00 + j01 * j01, b = j00 * j10 + j01 * j11, d = j10 * j10 + j11 * j11;
                float hx = P.R, hy = P.R;
                if (P.Antialias && (a > 1f || d > 1f))
                {
                    float tr = 0.5f * (a + d), disc = MathF.Sqrt(0.25f * (a - d) * (a - d) + b * b);
                    float l1 = tr + disc, l2 = MathF.Max(tr - disc, 0f);
                    float e1x, e1y;
                    if (MathF.Abs(b) > 1e-12f) { e1x = b; e1y = l1 - a; }
                    else if (a >= d) { e1x = 1f; e1y = 0f; }
                    else { e1x = 0f; e1y = 1f; }
                    float n1 = MathF.Sqrt(e1x * e1x + e1y * e1y);
                    e1x /= n1;
                    e1y /= n1;
                    float e2x = -e1y, e2y = e1x;
                    float c1 = Widen(l1, P.AaStrength), c2 = Widen(l2, P.AaStrength);
                    hx = P.R * MathF.Sqrt(c1 * e1x * e1x + c2 * e2x * e2x);
                    hy = P.R * MathF.Sqrt(c1 * e1y * e1y + c2 * e2y * e2y);
                }
                float invSx = 1f, invSy = 1f, apx = 0f, apy = 0f;
                if (P.Antialias)
                {
                    invSx = 1f / MathF.Sqrt(Widen(a, P.AaStrength));
                    invSy = 1f / MathF.Sqrt(Widen(d, P.AaStrength));
                }
                if (P.Aperture != 0f)
                {
                    apx = P.Aperture * ApertureCoefficient(w, MathF.Sqrt(a));
                    apy = P.Aperture * ApertureCoefficient(w, MathF.Sqrt(d));
                    if (P.SharpenOnly)
                    {
                        if (apx < 0f) apx = 0f;
                        if (apy < 0f) apy = 0f;
                    }
                }

                int ox0 = (int)MathF.Ceiling(fx - hx), ox1 = (int)MathF.Floor(fx + hx);
                int oy0 = (int)MathF.Ceiling(fy - hy), oy1 = (int)MathF.Floor(fy + hy);
                int nTx = Math.Min(ox1 - ox0 + 1, MaxTaps), nTy = oy1 - oy0 + 1;
                var wx = new float[MaxTaps];
                var xo = new int[MaxTaps];
                float sumX = 0f;
                for (int k = 0; k < nTx; k++)
                {
                    wx[k] = Weight1D(w, (float)(ox0 + k) - fx, invSx, apx);
                    sumX += wx[k];
                    xo[k] = Math.Clamp(ix + ox0 + k, 0, w.SrcW - 1);
                }

                var acc = new float[4];
                var col = new float[4];
                float wsum = 0f;
                if (!alpha && (ch == 3 || ch == 1))
                {
                    float sumY = 0f;
                    var h = new float[4];
                    for (int r = 0; r < nTy; r++)
                    {
                        float wyr = Weight1D(w, (float)(oy0 + r) - fy, invSy, apy);
                        if (wyr == 0f) continue;
                        sumY += wyr;
                        int rb = Math.Clamp(iy + oy0 + r, 0, w.SrcH - 1) * w.SrcW;
                        h[0] = h[1] = h[2] = h[3] = 0f;
                        for (int k2 = 0; k2 < nTx; k2++)
                        {
                            LoadSamples(w, (rb + xo[k2]) * ch, ch, col);
                            for (int c = 0; c < 4; c++) h[c] += wx[k2] * col[c];
                        }
                        for (int c = 0; c < 4; c++) acc[c] += wyr * h[c];
                    }
                    wsum = sumX * sumY;
                }
                else
                {
                    var add = new float[4];
                    for (int r = 0; r < nTy; r++)
                    {
                        float wyy = Weight1D(w, (float)(oy0 + r) - fy, invSy, apy);
                        if (wyy == 0f) continue;
                        int rb = Math.Clamp(iy + oy0 + r, 0, w.SrcH - 1) * w.SrcW;
                        for (int k2 = 0; k2 < nTx; k2++)
                        {
                            float wt = wyy * wx[k2];
                            if (wt == 0f) continue;
                            int si = (rb + xo[k2]) * ch;
                            float al = alpha ? SampleAt(w, si + cc) * P.AlphaScale : 1f;
                            LoadSamples(w, si, cc, col);
                            for (int c = 0; c < 4; c++) add[c] = (wt * col[c]) * al;
                            if (alpha) add[cc] = wt * al;
                            for (int c = 0; c < 4; c++) acc[c] += add[c];
                            wsum += wt;
                        }
                    }
                }
                float inv = wsum > 1e-12f ? 1f / wsum : 0f;
                for (int c = 0; c < 4; c++)
                {
                    val[c] = acc[c] * inv;
                    if (cover < 1f) val[c] = alpha ? val[c] * cover : val[c] * cover + P.Bg[c] * (1f - cover);
                }
            }

            if (P.Light == 2)
            {
                if (alpha)
                {
                    float al2 = val[cc];
                    if (!(al2 > 1e-7f))
                    {
                        for (int c = 0; c < 4; c++) val[c] = 0f;
                        val[cc] = al2;
                    }
                    else
                    {
                        float invA = 1f / al2;
                        for (int c = 0; c < cc; c++) val[c] = FromSigmoid(w, val[c] * invA) * al2;
                    }
                }
                else
                {
                    for (int c = 0; c < ch; c++) val[c] = FromSigmoid(w, val[c]);
                }
            }
            bool linear = P.Light != 1;
            var q = new uint[4];
            float invAl = 1f;
            if (alpha)
            {
                float al3 = Clamp01(val[cc]);
                q[cc] = Quantize(w, al3);
                invAl = al3 > 1e-7f ? 1f / al3 : 0f;
            }
            for (int c = 0; c < cc; c++)
            {
                float v = val[c] * invAl;
                v = linear ? SrgbEncode(w, v) : Clamp01(v);
                q[c] = Quantize(w, v);
            }
            if (P.OutDepth == 8) outWords[o] = q[0] | (q[1] << 8) | (q[2] << 16) | (q[3] << 24);
            else
            {
                outWords[2 * o] = q[0] | (q[1] << 16);
                outWords[2 * o + 1] = q[2] | (q[3] << 16);
            }
        }

        /// <summary>A whole defish as the GPU does it: GpuPlan's job, the Warp kernel strip by strip, unpacked and turned back.</summary>
        public static PixelImage Defish(PixelImage upright, DefishGeometry g, DefishSettings s)
        {
            var job = GpuPlan.Job.Prepare(upright, g, s);
            var w = new WarpState { G = job.Grid, P = job.Params, Src = job.Samples, SrcW = job.Sensor.Width, SrcH = job.Sensor.Height, OutW = job.OutW, OutH = job.OutH };
            var output = job.NewOutput();
            int strip = job.StripRows, wpp = job.Params.WordsPerPixel;
            for (int y0 = 0; y0 < job.OutH; y0 += strip)
            {
                int rows = Math.Min(strip, job.OutH - y0);
                var words = new uint[rows * job.OutW * wpp];
                Parallel.For(0, rows, yl =>
                {
                    for (int x = 0; x < job.OutW; x++) WarpPixel(w, x, y0 + yl, words, yl * job.OutW + x);
                });
                GpuPlan.Unpack(words, output, y0, rows);
            }
            return job.Finish(output);
        }

        /// <summary>A geometry made the GPU's way: the model evaluated in floats at every node.</summary>
        public static DefishGeometry Geometry(DefishModel model, DefishAdjustment adj, double aspect, int nodesAcross = 1025)
        {
            var e = new GpuEmulator(GpuPlan.ModelParams.For(model, adj, aspect));
            return GpuPlan.GeometryFromMapped(model, adj, aspect, nodesAcross, e.MapPoints(GpuPlan.GridPoints(aspect, nodesAcross)));
        }
    }
}
