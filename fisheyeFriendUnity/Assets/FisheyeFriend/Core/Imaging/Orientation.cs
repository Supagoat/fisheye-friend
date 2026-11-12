using System.Threading.Tasks;

namespace FisheyeFriend.Core.Imaging
{
    /// <summary>
    /// EXIF orientation (TIFF tag 0x0112): how a camera says the stored pixels must be turned to look upright.
    /// Cameras store the sensor's pixels as they were read out and record how the camera was held, so a shot
    /// taken with the camera turned 90° is stored sideways with orientation 6 or 8.
    /// </summary>
    public static class Orientation
    {
        public const int Normal = 1;
        public const int MirrorHorizontal = 2;
        public const int Rotate180 = 3;
        public const int MirrorVertical = 4;
        /// <summary>Mirror across the top-left to bottom-right diagonal.</summary>
        public const int Transpose = 5;
        /// <summary>Stored sideways: turn 90° clockwise to display (the usual tag for a camera held grip-up).</summary>
        public const int Rotate90Clockwise = 6;
        /// <summary>Mirror across the top-right to bottom-left diagonal.</summary>
        public const int Transverse = 7;
        public const int Rotate90CounterClockwise = 8;

        public static bool IsValid(int orientation)
        {
            return orientation >= 1 && orientation <= 8;
        }

        /// <summary>True for the orientations that exchange width and height (a quarter turn, possibly mirrored).</summary>
        public static bool SwapsAxes(int orientation)
        {
            return orientation >= 5 && orientation <= 8;
        }

        /// <summary>The orientation that undoes <paramref name="orientation"/>.</summary>
        public static int Inverse(int orientation)
        {
            if (orientation == Rotate90Clockwise) return Rotate90CounterClockwise;
            if (orientation == Rotate90CounterClockwise) return Rotate90Clockwise;
            return orientation;
        }

        /// <summary>What <see cref="Apply"/> does to the stored pixels for this orientation, in words.</summary>
        public static string Describe(int orientation)
        {
            switch (orientation)
            {
                case Normal: return "left as stored";
                case MirrorHorizontal: return "mirrored left-right";
                case Rotate180: return "turned 180°";
                case MirrorVertical: return "flipped top-bottom";
                case Transpose: return "turned 90° and mirrored";
                case Rotate90Clockwise: return "turned 90° clockwise";
                case Transverse: return "turned 90° and mirrored";
                case Rotate90CounterClockwise: return "turned 90° counter-clockwise";
                default: return "no orientation tag";
            }
        }

        /// <summary>
        /// Where the point (<paramref name="x"/>, <paramref name="y"/>) of an image <paramref name="width"/> x
        /// <paramref name="height"/> ends up in <c>Apply(image, orientation)</c> (pixel centres at whole numbers).
        /// </summary>
        public static void MapPoint(int orientation, int width, int height, double x, double y, out double ox, out double oy)
        {
            double w1 = width - 1, h1 = height - 1;
            switch (orientation)
            {
                case MirrorHorizontal: ox = w1 - x; oy = y; break;
                case Rotate180: ox = w1 - x; oy = h1 - y; break;
                case MirrorVertical: ox = x; oy = h1 - y; break;
                case Transpose: ox = y; oy = x; break;
                case Rotate90Clockwise: ox = h1 - y; oy = x; break;
                case Transverse: ox = h1 - y; oy = w1 - x; break;
                case Rotate90CounterClockwise: ox = y; oy = w1 - x; break;
                default: ox = x; oy = y; break;
            }
        }

        /// <summary>
        /// Returns <paramref name="src"/> transformed as <paramref name="orientation"/> prescribes, i.e. the image as it
        /// should be displayed. Orientation 1 or an invalid value returns <paramref name="src"/> itself.
        /// </summary>
        public static PixelImage Apply(PixelImage src, int orientation)
        {
            if (src == null || !IsValid(orientation) || orientation == Normal) return src;
            int w = src.Width, h = src.Height, ch = src.Channels;
            bool swap = SwapsAxes(orientation);
            int ow = swap ? h : w, oh = swap ? w : h;
            var dst = new PixelImage(ow, oh, ch, src.BitDepth)
            {
                Metadata = src.Metadata != null ? src.Metadata.Clone() : new ImageMetadata(),
                SourceFormat = src.SourceFormat,
            };

            // Source pixel for destination (u, v):  x = ax*u + bx*v + cx,  y = ay*u + by*v + cy.
            int ax, bx, cx, ay, by, cy;
            switch (orientation)
            {
                case MirrorHorizontal: ax = -1; bx = 0; cx = w - 1; ay = 0; by = 1; cy = 0; break;
                case Rotate180: ax = -1; bx = 0; cx = w - 1; ay = 0; by = -1; cy = h - 1; break;
                case MirrorVertical: ax = 1; bx = 0; cx = 0; ay = 0; by = -1; cy = h - 1; break;
                case Transpose: ax = 0; bx = 1; cx = 0; ay = 1; by = 0; cy = 0; break;
                case Rotate90Clockwise: ax = 0; bx = 1; cx = 0; ay = -1; by = 0; cy = h - 1; break;
                case Transverse: ax = 0; bx = -1; cx = w - 1; ay = -1; by = 0; cy = h - 1; break;
                default: ax = 0; bx = -1; cx = w - 1; ay = 1; by = 0; cy = 0; break; // Rotate90CounterClockwise
            }
            int step = (ay * w + ax) * ch;
            Parallel.For(0, oh, v =>
            {
                int s = ((by * v + cy) * w + (bx * v + cx)) * ch;
                int d = v * ow * ch;
                if (src.BitDepth == 8)
                {
                    byte[] a = src.Data8, b = dst.Data8;
                    for (int u = 0; u < ow; u++, s += step, d += ch)
                        for (int c = 0; c < ch; c++) b[d + c] = a[s + c];
                }
                else
                {
                    ushort[] a = src.Data16, b = dst.Data16;
                    for (int u = 0; u < ow; u++, s += step, d += ch)
                        for (int c = 0; c < ch; c++) b[d + c] = a[s + c];
                }
            });
            return dst;
        }
    }
}
