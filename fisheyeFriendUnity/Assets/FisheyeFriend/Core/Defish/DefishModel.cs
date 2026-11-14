using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>
    /// The defish as a physical model: a fisheye lens, and the projection the scene is redrawn with.
    ///
    /// <b>Lens.</b> A point of the photo at distance r from its centre (in photo widths) shows the direction at angle
    /// theta = sum Theta[i] r^(2i+1) from the optical axis, in the direction it lies from the centre (the photo's pixels
    /// are square and the lens is centred). The fitted lens covers 177 degrees across a 3:2 diagonal.
    ///
    /// <b>Projection.</b> A cylindrical projection of the Panini family. A direction at longitude lon (left-right angle)
    /// and latitude lat (up-down angle) is drawn at
    /// <code>x = Scale * (P + 1) sin(lon) / (P + cos(lon)) * (1 + sum Horizontal[i] sin^(2i+2) lon)
    /// y = Scale * (Pv + 1) / (Pv + cos(lon)) * tan(lat) * (1 + sum Vertical[i,j] sin^(2i+2) lat * sin^(2j) lon)</code>
    /// (in frame widths from the centre). x depends on the longitude only, so every vertical line of the scene comes
    /// out straight; y grows towards the sides, so horizontal lines keep a gentle curve, which is what keeps people
    /// at the edges looking natural. P = Panini, Pv = PaniniVertical; the sums are small corrections.
    ///
    /// The built-in model (<see cref="Hemi"/>) was fitted to Fisheye Friend's 596 hand-placed point pairs: it reproduces
    /// them to 0.9 px RMS on a 6000 px wide frame.
    /// </summary>
    public sealed class DefishModel
    {
        public string Name = "Hemi";
        public string Description = "";
        public double[] Theta;
        public double Scale, Panini, PaniniVertical;
        public double[] Horizontal = new double[0];
        public double[,] Vertical = new double[0, 0];

        /// <summary>The model fitted to Fisheye Friend's hand-made Hemi mapping (defish-hemi-model.json holds the same).</summary>
        public static DefishModel Hemi => new DefishModel
        {
            Name = "Hemi",
            Description = "Fitted to the 596 hand-placed point pairs of defish-hemi-596.json (0.9 px RMS).",
            Theta = new[] { 2.331570216256887, 0.24359368632643005, -0.4135582541599703, 4.3519581618197005 },
            Scale = 0.3591980688944844,
            Panini = 1.1247429558441109,
            PaniniVertical = 2.3143462334046765,
            Horizontal = new[] { 0.04343423280488648, -0.015420606868319116 },
            Vertical = new[,]
            {
                { -0.028188129209643746, -0.08957917940386168, 0.6556793921745786, -0.5603986482144145 },
                { -0.364788678123393, 0.6396620530198249, -3.4765387414545397, 2.870224602088231 },
                { 0.2497711231725929, -0.793538658554969, 4.4518517867127905, -3.4005833609102267 },
            },
        };

        // ------------------------------------------------------------------ lens

        /// <summary>Angle from the optical axis (radians) at distance <paramref name="r"/> (photo widths) from the centre.</summary>
        public double LensAngle(double r)
        {
            double r2 = r * r, p = r, s = 0;
            for (int i = 0; i < Theta.Length; i++, p *= r2) s += Theta[i] * p;
            return s;
        }

        public double LensAngleDerivative(double r)
        {
            double r2 = r * r, p = 1, s = 0;
            for (int i = 0; i < Theta.Length; i++, p *= r2) s += (2 * i + 1) * Theta[i] * p;
            return s;
        }

        /// <summary>The distance from the centre (photo widths) at which the lens shows angle <paramref name="theta"/>; <paramref name="guess"/> speeds it up.</summary>
        public double LensRadius(double theta, double guess = -1)
        {
            double r = guess > 0 ? guess : theta / Theta[0];
            for (int it = 0; it < 50; it++)
            {
                double f = LensAngle(r) - theta, d = LensAngleDerivative(r);
                if (d <= 1e-9) d = Theta[0];
                double step = f / d;
                if (step > 0.2) step = 0.2; else if (step < -0.2) step = -0.2;
                r -= step;
                if (r < 0) r = 0;
                if (Math.Abs(step) < 1e-14) break;
            }
            return r;
        }

        /// <summary>The lens' field of view across a frame of the given shape (height / width), in degrees, at a strength.</summary>
        public double DiagonalFieldOfView(double aspect, double strength = 1)
        {
            return 2 * LensAngle(Math.Sqrt(0.25 + 0.25 * aspect * aspect) * strength) * 180 / Math.PI;
        }

        public double HorizontalFieldOfView(double strength = 1)
        {
            return 2 * LensAngle(0.5 * strength) * 180 / Math.PI;
        }

        // ------------------------------------------------------------------ projection

        /// <summary>Horizontal position (without <see cref="Scale"/>) of longitude <paramref name="lon"/>, with its derivative.</summary>
        public double ProjectX(double lon, out double derivative)
        {
            double s = Math.Sin(lon), c = Math.Cos(lon), u = s * s;
            double den = Panini + c;
            double basis = (Panini + 1) * s / den;
            double dBasis = (Panini + 1) * (c * den + s * s) / (den * den);
            double corr = 1, dCorr = 0, prev = 1;          // prev = u^i
            for (int i = 0; i < Horizontal.Length; i++)
            {
                corr += Horizontal[i] * prev * u;
                dCorr += Horizontal[i] * (i + 1) * prev * 2 * s * c;   // d(u^(i+1))/dlon = (i+1) u^i * 2 sin cos
                prev *= u;
            }
            derivative = dBasis * corr + basis * dCorr;
            return basis * corr;
        }

        public double ProjectX(double lon)
        {
            double d;
            return ProjectX(lon, out d);
        }

        /// <summary>Vertical position (without <see cref="Scale"/>) of latitude <paramref name="lat"/> at longitude <paramref name="lon"/>, with d/dlat.</summary>
        public double ProjectY(double lon, double lat, out double dLat)
        {
            double sl = Math.Sin(lon), cl = Math.Cos(lon), u = sl * sl;
            double sp = Math.Sin(lat), cp = Math.Cos(lat), v = sp * sp;
            double side = (PaniniVertical + 1) / (PaniniVertical + cl);
            double tan = sp / cp, dTan = 1 / (cp * cp);
            double corr = 1, dCorr = 0;
            int ni = Vertical.GetLength(0), nj = Vertical.GetLength(1);
            double prev = 1;                                  // v^i
            for (int i = 0; i < ni; i++)
            {
                double row = 0, up = 1;
                for (int j = 0; j < nj; j++, up *= u) row += Vertical[i, j] * up;
                corr += row * prev * v;
                dCorr += row * (i + 1) * prev * 2 * sp * cp;  // d(v^(i+1))/dlat
                prev *= v;
            }
            dLat = side * (dTan * corr + tan * dCorr);
            return side * tan * corr;
        }

        // ------------------------------------------------------------------ files

        public const string Format = "fisheye-friend-defish-model";

        /// <summary>Reads a model file (see defish-hemi-model.json). Throws <see cref="FormatException"/>.</summary>
        public static DefishModel Parse(string json)
        {
            object root;
            try
            {
                root = Json.Parse(json);
            }
            catch (FormatException e)
            {
                throw new FormatException("The defish model is not valid JSON (" + e.Message + ").");
            }
            if (Json.String(root, "format", "") != Format) throw new FormatException("This file is not a Fisheye Friend defish model.");
            var lens = Json.Get(root, "lens");
            var proj = Json.Get(root, "projection");
            var m = new DefishModel
            {
                Name = Json.String(root, "name", "Model"),
                Description = Json.String(root, "description", ""),
                Theta = Numbers(Json.Get(lens, "theta"), "lens.theta"),
                Scale = Json.Number(proj, "scale", double.NaN),
                Panini = Json.Number(proj, "panini", double.NaN),
                PaniniVertical = Json.Number(proj, "paniniVertical", double.NaN),
                Horizontal = Json.Get(proj, "horizontal") != null ? Numbers(Json.Get(proj, "horizontal"), "projection.horizontal") : new double[0],
            };
            var rows = Json.Get(proj, "vertical") as List<object>;
            if (rows != null && rows.Count > 0)
            {
                var first = Numbers(rows[0], "projection.vertical");
                m.Vertical = new double[rows.Count, first.Length];
                for (int i = 0; i < rows.Count; i++)
                {
                    var r = Numbers(rows[i], "projection.vertical");
                    if (r.Length != first.Length) throw new FormatException("projection.vertical rows differ in length.");
                    for (int j = 0; j < r.Length; j++) m.Vertical[i, j] = r[j];
                }
            }
            if (m.Theta.Length == 0 || !(m.Theta[0] > 0)) throw new FormatException("lens.theta must start with a positive number.");
            if (!(m.Scale > 0) || double.IsNaN(m.Panini) || double.IsNaN(m.PaniniVertical) || m.Panini <= -0.9 || m.PaniniVertical <= -0.9)
                throw new FormatException("projection.scale, panini and paniniVertical are missing or out of range.");
            return m;
        }

        private static double[] Numbers(object o, string what)
        {
            var list = o as List<object>;
            if (list == null) throw new FormatException(what + " must be a list of numbers.");
            var r = new double[list.Count];
            for (int i = 0; i < r.Length; i++)
            {
                if (!(list[i] is double)) throw new FormatException(what + " must be a list of numbers.");
                r[i] = (double)list[i];
            }
            return r;
        }

        public string ToJson()
        {
            var ci = CultureInfo.InvariantCulture;
            Func<double, string> n = v => v.ToString("R", ci);
            var sb = new StringBuilder();
            sb.Append("{\n  \"format\": \"").Append(Format).Append("\",\n  \"formatVersion\": 1,\n");
            sb.Append("  \"name\": \"").Append(Escape(Name)).Append("\",\n");
            sb.Append("  \"description\": \"").Append(Escape(Description)).Append("\",\n");
            sb.Append("  \"lens\": { \"theta\": [").Append(string.Join(", ", Array.ConvertAll(Theta, x => n(x)))).Append("] },\n");
            sb.Append("  \"projection\": {\n    \"scale\": ").Append(n(Scale)).Append(",\n    \"panini\": ").Append(n(Panini))
              .Append(",\n    \"paniniVertical\": ").Append(n(PaniniVertical)).Append(",\n    \"horizontal\": [")
              .Append(string.Join(", ", Array.ConvertAll(Horizontal, x => n(x)))).Append("],\n    \"vertical\": [");
            for (int i = 0; i < Vertical.GetLength(0); i++)
            {
                sb.Append(i == 0 ? "\n      [" : ",\n      [");
                for (int j = 0; j < Vertical.GetLength(1); j++) sb.Append(j == 0 ? "" : ", ").Append(n(Vertical[i, j]));
                sb.Append("]");
            }
            sb.Append("\n    ]\n  }\n}\n");
            return sb.ToString();
        }

        private static string Escape(string s)
        {
            return (s ?? "").Replace("\\", "\\\\").Replace("\"", "\\\"").Replace("\n", "\\n");
        }
    }

    /// <summary>
    /// How the user adjusts the defish. All 1 (100%) is the model as fitted.
    /// </summary>
    public sealed class DefishAdjustment
    {
        /// <summary>Points on the shape curves: the centre, a third, two thirds of the way and the corner of the frame.</summary>
        public const int ShapePoints = 4;

        /// <summary>
        /// How strongly the lens bends: its focal length is divided by this. 1 = the lens the model was fitted to (a 177
        /// degree diagonal fisheye on a 3:2 frame); lower for a milder fisheye, higher for a stronger one.
        /// </summary>
        public double Strength = 1;
        /// <summary>How much of the horizontal correction to apply (0 = none: the photo's own left-right shape; 1 = full).</summary>
        public double AmountX = 1;
        /// <summary>How much of the vertical correction to apply.</summary>
        public double AmountY = 1;
        /// <summary>The horizontal amount across the frame, as a factor at <see cref="ShapePoints"/> distances from the centre.</summary>
        public double[] ShapeX = { 1, 1, 1, 1 };
        public double[] ShapeY = { 1, 1, 1, 1 };

        public const double MinStrength = 0.5, MaxStrength = 1.15, MaxAmount = 1.5, MaxShape = 2.0;

        public DefishAdjustment Clone()
        {
            return new DefishAdjustment
            {
                Strength = Strength,
                AmountX = AmountX,
                AmountY = AmountY,
                ShapeX = (double[])ShapeX.Clone(),
                ShapeY = (double[])ShapeY.Clone(),
            };
        }

        public void Validate()
        {
            Strength = Clamp(Strength, MinStrength, MaxStrength, 1);
            AmountX = Clamp(AmountX, 0, MaxAmount, 1);
            AmountY = Clamp(AmountY, 0, MaxAmount, 1);
            ShapeX = ValidShape(ShapeX);
            ShapeY = ValidShape(ShapeY);
        }

        private static double[] ValidShape(double[] s)
        {
            var r = new double[ShapePoints];
            for (int i = 0; i < ShapePoints; i++) r[i] = Clamp(s != null && i < s.Length ? s[i] : 1, 0, MaxShape, 1);
            return r;
        }

        private static double Clamp(double v, double lo, double hi, double fallback)
        {
            if (double.IsNaN(v) || double.IsInfinity(v)) return fallback;
            return v < lo ? lo : v > hi ? hi : v;
        }

        public bool IsDefault
        {
            get
            {
                if (Math.Abs(Strength - 1) > 1e-9 || Math.Abs(AmountX - 1) > 1e-9 || Math.Abs(AmountY - 1) > 1e-9) return false;
                for (int i = 0; i < ShapePoints; i++)
                    if (Math.Abs(ShapeX[i] - 1) > 1e-9 || Math.Abs(ShapeY[i] - 1) > 1e-9) return false;
                return true;
            }
        }

        public bool ShapesFlat
        {
            get
            {
                for (int i = 0; i < ShapePoints; i++)
                    if (Math.Abs(ShapeX[i] - 1) > 1e-9 || Math.Abs(ShapeY[i] - 1) > 1e-9) return false;
                return true;
            }
        }

        public bool Same(DefishAdjustment o)
        {
            return o != null && Key == o.Key;
        }

        /// <summary>Every value, precisely enough to tell two adjustments apart (a cache key).</summary>
        public string Key
        {
            get
            {
                var ci = CultureInfo.InvariantCulture;
                var sb = new StringBuilder();
                sb.Append(Strength.ToString("R", ci)).Append('|').Append(AmountX.ToString("R", ci)).Append('|').Append(AmountY.ToString("R", ci));
                foreach (var v in ShapeX) sb.Append('|').Append(v.ToString("R", ci));
                foreach (var v in ShapeY) sb.Append('|').Append(v.ToString("R", ci));
                return sb.ToString();
            }
        }

        /// <summary>
        /// A shape curve's factor at distance <paramref name="rho"/> from the centre (0 = centre, 1 = corner of the
        /// frame): a smooth monotone-preserving cubic (Fritsch-Carlson) through the points, flat beyond the corner.
        /// </summary>
        public static double ShapeAt(double[] shape, double rho)
        {
            int n = shape.Length;
            if (rho <= 0) return shape[0];
            if (rho >= 1) return shape[n - 1];
            double h = 1.0 / (n - 1);
            int k = Math.Min(n - 2, (int)(rho / h));
            double t = (rho - k * h) / h;
            double m0 = Tangent(shape, k, h), m1 = Tangent(shape, k + 1, h);
            double t2 = t * t, t3 = t2 * t;
            return (2 * t3 - 3 * t2 + 1) * shape[k] + (t3 - 2 * t2 + t) * h * m0 + (-2 * t3 + 3 * t2) * shape[k + 1] + (t3 - t2) * h * m1;
        }

        private static double Tangent(double[] y, int k, double h)
        {
            int n = y.Length;
            if (k == 0) return 0;                                 // flat at the centre: the curve is symmetric about it
            if (k == n - 1) return (y[k] - y[k - 1]) / h;
            double d0 = (y[k] - y[k - 1]) / h, d1 = (y[k + 1] - y[k]) / h;
            if (d0 * d1 <= 0) return 0;
            return 2 / (1 / d0 + 1 / d1);                         // harmonic mean: no overshoot
        }
    }

    /// <summary>
    /// The model with an adjustment, for a frame of a given shape: for each point of the defished image, the point of
    /// the photo to read. Coordinates are normalised (0..1 across the width and the height of each image).
    /// Not thread-safe (it keeps warm starts between neighbouring points); use one per thread.
    /// </summary>
    public sealed class DefishModelMap
    {
        private readonly DefishModel _m;
        private readonly double _aspect, _k, _strength, _ax, _ay, _rhoNorm;
        private readonly double[] _shapeX, _shapeY;
        private readonly bool _flat;
        private double _lon, _lat, _r;                // warm starts

        /// <param name="aspect">The frame's height / width (as the lens saw it).</param>
        public DefishModelMap(DefishModel model, DefishAdjustment adjustment, double aspect)
        {
            _m = model;
            var a = (adjustment ?? new DefishAdjustment()).Clone();
            a.Validate();
            _aspect = aspect;
            _strength = a.Strength;
            _ax = a.AmountX;
            _ay = a.AmountY;
            _shapeX = a.ShapeX;
            _shapeY = a.ShapeY;
            _flat = a.ShapesFlat;
            // Fit the width: the photo's left and right edge midpoints land on the result's edges, as with the fitted model.
            _k = model.Scale * model.ProjectX(model.LensAngle(0.5)) / model.ProjectX(model.LensAngle(0.5 * _strength));
            _rhoNorm = 1 / Math.Sqrt(0.25 + 0.25 * aspect * aspect);
            Reset();
        }

        /// <summary>The photo point (normalised) shown at defished point (<paramref name="u"/>, <paramref name="v"/>) (normalised).</summary>
        public void Map(double u, double v, out double s, out double t)
        {
            double x = u - 0.5, y = (v - 0.5) * _aspect;      // frame widths from the centre
            double sx, sy;
            Physical(x, y, out sx, out sy);
            double ax = _ax, ay = _ay;
            if (!_flat)
            {
                double rho = Math.Sqrt(x * x + y * y) * _rhoNorm;
                ax *= DefishAdjustment.ShapeAt(_shapeX, rho);
                ay *= DefishAdjustment.ShapeAt(_shapeY, rho);
            }
            sx = x + ax * (sx - x);
            sy = y + ay * (sy - y);
            s = sx + 0.5;
            t = sy / _aspect + 0.5;
        }

        /// <summary>Starts the warm starts afresh (call between unrelated points).</summary>
        public void Reset()
        {
            _lon = double.NaN;
            _lat = double.NaN;
            _r = -1;
        }

        /// <summary>The full model: defished point -> photo point, in frame widths from the centres.</summary>
        private void Physical(double x, double y, out double sx, out double sy)
        {
            var m = _m;
            // longitude: k ProjectX(lon) = x (increasing in lon over the range used)
            double lon = double.IsNaN(_lon) ? Math.Atan(x / _k) : _lon;
            for (int it = 0; it < 60; it++)
            {
                double d, f = _k * m.ProjectX(lon, out d) - x;
                d *= _k;
                if (!(d > 1e-9)) d = _k;
                double step = f / d;
                if (step > 0.2) step = 0.2; else if (step < -0.2) step = -0.2;
                lon -= step;
                if (lon > 1.55) lon = 1.55; else if (lon < -1.55) lon = -1.55;
                if (Math.Abs(step) < 1e-13) break;
            }
            _lon = lon;
            // latitude: k ProjectY(lon, lat) = y (increasing in lat)
            double lat = double.IsNaN(_lat) ? Math.Atan(y / (_k * 1.2)) : _lat;
            for (int it = 0; it < 60; it++)
            {
                double d, f = _k * m.ProjectY(lon, lat, out d) - y;
                d *= _k;
                if (!(d > 1e-9)) d = _k;
                double step = f / d;
                if (step > 0.2) step = 0.2; else if (step < -0.2) step = -0.2;
                lat -= step;
                if (lat > 1.5) lat = 1.5; else if (lat < -1.5) lat = -1.5;
                if (Math.Abs(step) < 1e-13) break;
            }
            _lat = lat;
            // the direction, and where the lens puts it
            double cp = Math.Cos(lat);
            double X = cp * Math.Sin(lon), Y = Math.Sin(lat), Z = cp * Math.Cos(lon);
            double sinTh = Math.Sqrt(X * X + Y * Y);
            double theta = Math.Atan2(sinTh, Z);
            if (theta < 1e-9)
            {
                double inv = 1 / (m.Theta[0] * _strength);
                sx = X * inv;
                sy = Y * inv;
                _r = -1;
                return;
            }
            double r = m.LensRadius(theta, _r);
            _r = r;
            r /= _strength;
            sx = r * X / sinTh;
            sy = r * Y / sinTh;
        }
    }
}
