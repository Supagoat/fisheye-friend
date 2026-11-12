using System;
using System.Threading.Tasks;
using FisheyeFriend.Core.Resampling;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>
    /// The defish projection: for a point of the defished image (normalised 0..1 across its width and height, pixel
    /// edges), the point of the fisheye photo to read, in the same normalised units.
    ///
    /// Defined by thin-plate splines through the mapping's point pairs (exactly Fisheye Friend's fit). Evaluating the 596-term
    /// splines for every pixel is what made the original "sharp" render take minutes, so the splines are evaluated
    /// exactly on a fine grid once (normalised, so one grid serves every image size) and interpolated bicubically
    /// between the nodes; see <see cref="MaxGridError"/> for the error that adds.
    /// </summary>
    public sealed class DefishGeometry
    {
        // grid over [U0, U0 + (NU-1) HU] x [V0, V0 + (NV-1) HV], generous margin around the 0..1 frame
        private const double Margin = 0.06;
        private readonly int _nu, _nv;
        private readonly double _u0, _v0, _hu, _hv;
        private readonly float[] _s, _t;          // normalised source x / y at each node (double precision is not needed: ~1e-7)
        private readonly ThinPlateSpline _tx, _ty;

        public HemiMapping Mapping { get; }

        public DefishGeometry(HemiMapping mapping, double lambda = 1e-5, int nodesAcross = 1025)
        {
            Mapping = mapping;
            mapping.FitInverse(lambda, out _tx, out _ty);
            double aspect = mapping.TargetHeight > 0 && mapping.TargetWidth > 0 ? (double)mapping.TargetHeight / mapping.TargetWidth : 2.0 / 3.0;
            _hu = 1.0 / (nodesAcross - 1);
            _hv = _hu / aspect;                     // square cells in pixels for images of the mapping's shape
            _u0 = -Margin;
            _v0 = -Margin;
            _nu = (int)Math.Ceiling((1 + 2 * Margin) / _hu) + 1;
            _nv = (int)Math.Ceiling((1 + 2 * Margin) / _hv) + 1;
            _s = new float[_nu * _nv];
            _t = new float[_nu * _nv];
            Parallel.For(0, _nv, j =>
            {
                double v = _v0 + j * _hv;
                for (int i = 0; i < _nu; i++)
                {
                    double u = _u0 + i * _hu;
                    _s[j * _nu + i] = (float)_tx.Eval(u, v);
                    _t[j * _nu + i] = (float)_ty.Eval(u, v);
                }
            });
        }

        public int GridNodes => _nu * _nv;

        /// <summary>The exact spline (slow: hundreds of terms).</summary>
        public void MapExact(double u, double v, out double s, out double t)
        {
            s = _tx.Eval(u, v);
            t = _ty.Eval(u, v);
        }

        /// <summary>The interpolated map with its derivatives (normalised units).</summary>
        public void Map(double u, double v, out double s, out double t, out double dsdu, out double dsdv, out double dtdu, out double dtdv)
        {
            double fu = (u - _u0) / _hu, fv = (v - _v0) / _hv;
            int i = (int)Math.Floor(fu), j = (int)Math.Floor(fv);
            if (i < 1) i = 1; else if (i > _nu - 3) i = _nu - 3;
            if (j < 1) j = 1; else if (j > _nv - 3) j = _nv - 3;
            double a = fu - i, b = fv - j;            // outside the grid this extrapolates the edge cubic
            double wa0, wa1, wa2, wa3, da0, da1, da2, da3, wb0, wb1, wb2, wb3, db0, db1, db2, db3;
            CatmullRom(a, out wa0, out wa1, out wa2, out wa3, out da0, out da1, out da2, out da3);
            CatmullRom(b, out wb0, out wb1, out wb2, out wb3, out db0, out db1, out db2, out db3);
            double s0 = 0, s1 = 0, s2 = 0, t0 = 0, t1 = 0, t2 = 0;
            for (int r = 0; r < 4; r++)
            {
                int o = (j - 1 + r) * _nu + i - 1;
                double rs = wa0 * _s[o] + wa1 * _s[o + 1] + wa2 * _s[o + 2] + wa3 * _s[o + 3];
                double rds = da0 * _s[o] + da1 * _s[o + 1] + da2 * _s[o + 2] + da3 * _s[o + 3];
                double rt = wa0 * _t[o] + wa1 * _t[o + 1] + wa2 * _t[o + 2] + wa3 * _t[o + 3];
                double rdt = da0 * _t[o] + da1 * _t[o + 1] + da2 * _t[o + 2] + da3 * _t[o + 3];
                double wb = r == 0 ? wb0 : r == 1 ? wb1 : r == 2 ? wb2 : wb3;
                double db = r == 0 ? db0 : r == 1 ? db1 : r == 2 ? db2 : db3;
                s0 += wb * rs; s1 += wb * rds; s2 += db * rs;
                t0 += wb * rt; t1 += wb * rdt; t2 += db * rt;
            }
            s = s0;
            t = t0;
            dsdu = s1 / _hu;
            dtdu = t1 / _hu;
            dsdv = s2 / _hv;
            dtdv = t2 / _hv;
        }

        /// <summary>
        /// The inverse: the point of the defished frame that shows photo point (<paramref name="s"/>, <paramref name="t"/>)
        /// (Newton's method on <see cref="Map"/>). False if it did not converge, e.g. for a point the frame does not show.
        /// </summary>
        public bool Unmap(double s, double t, out double u, out double v)
        {
            u = s;
            v = t;
            for (int it = 0; it < 30; it++)
            {
                double ms, mt, a, b, c, d;
                Map(u, v, out ms, out mt, out a, out b, out c, out d);
                double es = ms - s, et = mt - t;
                if (Math.Abs(es) < 1e-10 && Math.Abs(et) < 1e-10) return true;
                double det = a * d - b * c;
                if (Math.Abs(det) < 1e-12) return false;
                double du = (d * es - b * et) / det, dv = (a * et - c * es) / det;
                double step = Math.Max(Math.Abs(du), Math.Abs(dv));
                if (step > 0.05)
                {
                    du *= 0.05 / step;
                    dv *= 0.05 / step;
                }
                u -= du;
                v -= dv;
            }
            double fs, ft, a2, b2, c2, d2;
            Map(u, v, out fs, out ft, out a2, out b2, out c2, out d2);
            return Math.Abs(fs - s) < 1e-7 && Math.Abs(ft - t) < 1e-7;
        }

        private static void CatmullRom(double x, out double w0, out double w1, out double w2, out double w3,
            out double d0, out double d1, out double d2, out double d3)
        {
            double x2 = x * x, x3 = x2 * x;
            w0 = 0.5 * (-x3 + 2 * x2 - x);
            w1 = 0.5 * (3 * x3 - 5 * x2 + 2);
            w2 = 0.5 * (-3 * x3 + 4 * x2 + x);
            w3 = 0.5 * (x3 - x2);
            d0 = 0.5 * (-3 * x2 + 4 * x - 1);
            d1 = 0.5 * (9 * x2 - 10 * x);
            d2 = 0.5 * (-9 * x2 + 8 * x + 1);
            d3 = 0.5 * (3 * x2 - 2 * x);
        }

        /// <summary>Largest difference between the interpolated and the exact map over the frame, in pixels of a photo <paramref name="width"/> x <paramref name="height"/>.</summary>
        public double MaxGridError(int width, int height, int samples = 20000, int seed = 1)
        {
            var rnd = new Random(seed);
            double worst = 0;
            for (int k = 0; k < samples; k++)
            {
                double u = rnd.NextDouble(), v = rnd.NextDouble(), s, t, es, et, a, b, c, d;
                Map(u, v, out s, out t, out a, out b, out c, out d);
                MapExact(u, v, out es, out et);
                worst = Math.Max(worst, Math.Max(Math.Abs(s - es) * width, Math.Abs(t - et) * height));
            }
            return worst;
        }
    }

    /// <summary>
    /// <see cref="DefishGeometry"/> as a <see cref="WarpMap"/> for a photo <c>srcW x srcH</c> and an output <c>outW x outH</c>
    /// showing the part <c>[u0, u1] x [v0, v1]</c> of the defished frame (normally all of it, 0..1).
    /// </summary>
    public sealed class DefishWarpMap : WarpMap
    {
        private readonly DefishGeometry _g;
        private readonly int _srcW, _srcH;
        private readonly double _u0, _v0, _du, _dv;

        public DefishWarpMap(DefishGeometry g, int srcW, int srcH, int outW, int outH, double u0 = 0, double v0 = 0, double u1 = 1, double v1 = 1)
        {
            _g = g;
            _srcW = srcW;
            _srcH = srcH;
            _u0 = u0;
            _v0 = v0;
            _du = (u1 - u0) / outW;
            _dv = (v1 - v0) / outH;
        }

        public override void Map(double x, double y, out double sx, out double sy, out double j00, out double j01, out double j10, out double j11)
        {
            double s, t, su, sv, tu, tv;
            _g.Map(_u0 + (x + 0.5) * _du, _v0 + (y + 0.5) * _dv, out s, out t, out su, out sv, out tu, out tv);
            sx = s * _srcW - 0.5;
            sy = t * _srcH - 0.5;
            j00 = su * _du * _srcW;
            j01 = sv * _dv * _srcW;
            j10 = tu * _du * _srcH;
            j11 = tv * _dv * _srcH;
        }
    }
}
