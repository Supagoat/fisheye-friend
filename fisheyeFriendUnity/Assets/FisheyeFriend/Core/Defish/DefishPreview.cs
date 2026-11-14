using System;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Defish
{
    /// <summary>
    /// A quick, small picture of the defished frame for the settings panel: a reduced copy of the photo is warped with a
    /// coarse geometry and sampled bilinearly. It shows the whole frame as <see cref="Defisher.Process"/> would (turned
    /// into the lens' landscape frame and back) in a few tens of milliseconds, so it can follow a slider as it moves.
    /// </summary>
    public static class DefishPreview
    {
        /// <summary>
        /// Renders the defished frame of an upright photo.
        /// </summary>
        /// <param name="src">The upright photo, reduced to any size: RGBA32, top row first.</param>
        /// <param name="turn">The photo's turn into the lens' frame (<see cref="Defisher.SensorTurn"/> of the full photo).</param>
        /// <param name="g">A geometry made for the photo's <see cref="Defisher.SensorAspect"/> (a coarse one is plenty).</param>
        /// <param name="width">Width of the picture (upright, like the photo).</param>
        /// <param name="bottomUp">Write the bottom row first, as a Unity texture expects.</param>
        /// <returns>RGBA32, <paramref name="width"/> x <paramref name="height"/>.</returns>
        public static byte[] Render(byte[] src, int srcWidth, int srcHeight, int turn, DefishGeometry g, int width, int height,
            bool bottomUp = false, CancellationToken ct = default)
        {
            if (src == null) throw new ArgumentNullException(nameof(src));
            if (g == null) throw new ArgumentNullException(nameof(g));
            if (srcWidth < 1 || srcHeight < 1 || src.Length < (long)srcWidth * srcHeight * 4) throw new ArgumentException("bad source size");
            if (width < 1 || height < 1) throw new ArgumentOutOfRangeException(nameof(width));
            var dst = new byte[(long)width * height * 4];
            int back = Orientation.Inverse(turn);
            var po = new ParallelOptions { CancellationToken = ct };
            Parallel.For(0, height, po, oy =>
            {
                int row = bottomUp ? height - 1 - oy : oy;
                long o = (long)row * width * 4;
                double b = (oy + 0.5) / height;
                for (int ox = 0; ox < width; ox++, o += 4)
                {
                    double u, v, s, t, su, sv, tu, tv, p, q;
                    Turn(turn, (ox + 0.5) / width, b, out u, out v);          // upright result -> lens frame
                    g.Map(u, v, out s, out t, out su, out sv, out tu, out tv); // -> photo, lens frame
                    Turn(back, s, t, out p, out q);                            // -> upright photo
                    Sample(src, srcWidth, srcHeight, p * srcWidth - 0.5, q * srcHeight - 0.5, dst, o);
                }
            });
            return dst;
        }

        /// <summary>An orientation applied to normalised coordinates (0..1 across the frame).</summary>
        public static void Turn(int orientation, double x, double y, out double ox, out double oy)
        {
            // Orientation.MapPoint works on pixel centres of a w x h image; a square N x N one makes it exact on fractions.
            const double N = 4096;
            double px, py;
            Orientation.MapPoint(orientation, (int)N, (int)N, x * N - 0.5, y * N - 0.5, out px, out py);
            ox = (px + 0.5) / N;
            oy = (py + 0.5) / N;
        }

        private static void Sample(byte[] src, int w, int h, double x, double y, byte[] dst, long o)
        {
            // past the photo: repeat its edge (as "Full frame" does)
            if (x < 0) x = 0; else if (x > w - 1) x = w - 1;
            if (y < 0) y = 0; else if (y > h - 1) y = h - 1;
            int x0 = (int)x, y0 = (int)y;
            int x1 = x0 + 1 < w ? x0 + 1 : x0, y1 = y0 + 1 < h ? y0 + 1 : y0;
            double fx = x - x0, fy = y - y0;
            long i00 = ((long)y0 * w + x0) * 4, i01 = ((long)y0 * w + x1) * 4, i10 = ((long)y1 * w + x0) * 4, i11 = ((long)y1 * w + x1) * 4;
            for (int c = 0; c < 3; c++)
            {
                double top = src[i00 + c] + (src[i01 + c] - src[i00 + c]) * fx;
                double bot = src[i10 + c] + (src[i11 + c] - src[i10 + c]) * fx;
                dst[o + c] = (byte)(top + (bot - top) * fy + 0.5);
            }
            dst[o + 3] = 255;
        }
    }
}
