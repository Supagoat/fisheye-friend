using System;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>
    /// A 2-D thin-plate spline f(x, y) = a0 + ax x + ay y + sum w_i U(|p - p_i|), U(r) = r^2 ln r, fitted to scattered
    /// values with smoothing <c>lambda</c> (Bookstein 1989). Same formulation as Fisheye Friend's original fit, so the
    /// defish geometry is unchanged; adds the analytic gradient.
    /// </summary>
    public sealed class ThinPlateSpline
    {
        private readonly double[] _x, _y, _w;
        private readonly double _a0, _ax, _ay;

        private ThinPlateSpline(double[] x, double[] y, double[] w, double a0, double ax, double ay)
        {
            _x = x;
            _y = y;
            _w = w;
            _a0 = a0;
            _ax = ax;
            _ay = ay;
        }

        public int Count => _x.Length;

        public static ThinPlateSpline Fit(double[] xi, double[] yi, double[] zi, double lambda)
        {
            int n = xi.Length, dim = n + 3;
            var A = new double[dim, dim];
            var b = new double[dim];
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    double dx = xi[i] - xi[j], dy = yi[i] - yi[j];
                    double v = U(Math.Sqrt(dx * dx + dy * dy)) + (i == j ? lambda : 0);
                    A[i, j] = v;
                    A[j, i] = v;
                }
            for (int i = 0; i < n; i++)
            {
                A[i, n] = A[n, i] = 1;
                A[i, n + 1] = A[n + 1, i] = xi[i];
                A[i, n + 2] = A[n + 2, i] = yi[i];
                b[i] = zi[i];
            }
            var sol = Solve(A, b);
            var w = new double[n];
            Array.Copy(sol, w, n);
            return new ThinPlateSpline((double[])xi.Clone(), (double[])yi.Clone(), w, sol[n], sol[n + 1], sol[n + 2]);
        }

        public double Eval(double px, double py)
        {
            double s = _a0 + _ax * px + _ay * py;
            for (int i = 0; i < _x.Length; i++)
            {
                double dx = px - _x[i], dy = py - _y[i];
                double r2 = dx * dx + dy * dy;
                if (r2 > 1e-24) s += _w[i] * 0.5 * r2 * Math.Log(r2);      // r^2 ln r = r^2 ln(r^2) / 2
            }
            return s;
        }

        /// <summary>Value and gradient.</summary>
        public double Eval(double px, double py, out double dfdx, out double dfdy)
        {
            double s = _a0 + _ax * px + _ay * py;
            double gx = _ax, gy = _ay;
            for (int i = 0; i < _x.Length; i++)
            {
                double dx = px - _x[i], dy = py - _y[i];
                double r2 = dx * dx + dy * dy;
                if (r2 <= 1e-24) continue;
                double l = Math.Log(r2);
                s += _w[i] * 0.5 * r2 * l;
                double g = _w[i] * (l + 1);                                 // d/dx (r^2 ln r) = x (2 ln r + 1)
                gx += g * dx;
                gy += g * dy;
            }
            dfdx = gx;
            dfdy = gy;
            return s;
        }

        private static double U(double r)
        {
            return r <= 1e-12 ? 0 : r * r * Math.Log(r);
        }

        /// <summary>Gaussian elimination with partial pivoting.</summary>
        private static double[] Solve(double[,] A, double[] b)
        {
            int n = b.Length;
            var M = (double[,])A.Clone();
            var rhs = (double[])b.Clone();
            for (int k = 0; k < n; k++)
            {
                int piv = k;
                double mx = Math.Abs(M[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    double v = Math.Abs(M[i, k]);
                    if (v > mx) { mx = v; piv = i; }
                }
                if (mx < 1e-14) throw new InvalidOperationException("The thin-plate spline system is singular (duplicate points with no smoothing?).");
                if (piv != k)
                {
                    for (int j = 0; j < n; j++)
                    {
                        double t = M[k, j];
                        M[k, j] = M[piv, j];
                        M[piv, j] = t;
                    }
                    double tr = rhs[k];
                    rhs[k] = rhs[piv];
                    rhs[piv] = tr;
                }
                double akk = M[k, k];
                for (int i = k + 1; i < n; i++)
                {
                    double f = M[i, k] / akk;
                    if (f == 0) continue;
                    rhs[i] -= f * rhs[k];
                    for (int j = k; j < n; j++) M[i, j] -= f * M[k, j];
                }
            }
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double s = rhs[i];
                for (int j = i + 1; j < n; j++) s -= M[i, j] * x[j];
                x[i] = s / M[i, i];
            }
            return x;
        }
    }
}
