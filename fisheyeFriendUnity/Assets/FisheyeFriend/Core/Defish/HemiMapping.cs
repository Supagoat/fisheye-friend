using System;
using System.Collections.Generic;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>
    /// The defish geometry as measured point pairs: where points of the fisheye photo ("src") end up in the defished image
    /// ("dst"), in coordinates normalised to each image's width and height (0..1, pixel edges). Fisheye Friend's mapping
    /// file format:
    /// <code>{ "source": {"width": 6000, "height": 4000}, "target": {...}, "pairs": [ {"id": 1, "src": {"x":..,"y":..,"xn":..,"yn":..}, "dst": {...}} ] }</code>
    /// (xn/yn are optional; x/y are then divided by the image size).
    /// </summary>
    public sealed class HemiMapping
    {
        public struct Pair
        {
            public double SrcX, SrcY, DstX, DstY;
        }

        public int SourceWidth, SourceHeight, TargetWidth, TargetHeight;
        public readonly List<Pair> Pairs = new List<Pair>();

        public static HemiMapping Parse(string json)
        {
            object root;
            try
            {
                root = Json.Parse(json);
            }
            catch (FormatException e)
            {
                throw new FormatException("The defish mapping is not valid JSON (" + e.Message + ").");
            }
            var m = new HemiMapping
            {
                SourceWidth = (int)Json.Number(Json.Get(root, "source"), "width", 0),
                SourceHeight = (int)Json.Number(Json.Get(root, "source"), "height", 0),
                TargetWidth = (int)Json.Number(Json.Get(root, "target"), "width", 0),
                TargetHeight = (int)Json.Number(Json.Get(root, "target"), "height", 0),
            };
            var pairs = Json.Get(root, "pairs") as List<object>;
            if (pairs == null) throw new FormatException("The defish mapping has no \"pairs\".");
            foreach (var po in pairs)
            {
                var s = Json.Get(po, "src");
                var d = Json.Get(po, "dst");
                if (s == null || d == null) continue;
                var p = new Pair
                {
                    SrcX = Norm(s, "xn", "x", m.SourceWidth),
                    SrcY = Norm(s, "yn", "y", m.SourceHeight),
                    DstX = Norm(d, "xn", "x", m.TargetWidth),
                    DstY = Norm(d, "yn", "y", m.TargetHeight),
                };
                if (double.IsNaN(p.SrcX + p.SrcY + p.DstX + p.DstY)) continue;
                m.Pairs.Add(p);
            }
            if (m.Pairs.Count < 3) throw new FormatException("The defish mapping needs at least 3 point pairs.");
            return m;
        }

        private static double Norm(object o, string nKey, string pKey, int size)
        {
            var n = Json.Get(o, nKey);
            if (n is double) return (double)n;
            var p = Json.Get(o, pKey);
            if (p is double && size > 0) return (double)p / size;
            return double.NaN;
        }

        /// <summary>
        /// The splines from defished (dst) to fisheye (src) normalised coordinates: for each output pixel, where to read the
        /// photo. Fisheye Friend fits them with <paramref name="lambda"/> = 1e-5.
        /// </summary>
        public void FitInverse(double lambda, out ThinPlateSpline toSrcX, out ThinPlateSpline toSrcY)
        {
            int n = Pairs.Count;
            double[] xd = new double[n], yd = new double[n], xs = new double[n], ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xd[i] = Pairs[i].DstX;
                yd[i] = Pairs[i].DstY;
                xs[i] = Pairs[i].SrcX;
                ys[i] = Pairs[i].SrcY;
            }
            toSrcX = ThinPlateSpline.Fit(xd, yd, xs, lambda);
            toSrcY = ThinPlateSpline.Fit(xd, yd, ys, lambda);
        }
    }
}
