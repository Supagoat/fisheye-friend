using System;
using System.Threading.Tasks;
using FisheyeFriend.Core.Color;

namespace FisheyeFriend.Core.Imaging
{
    /// <summary>
    /// Builds RGBA32 previews for the viewer. Zooming out uses whole-number reductions: each screen
    /// pixel is the exact average (in linear light) of an n x n block of image pixels, so a 1/2 or 1/3
    /// view is still pixel-accurate rather than a GPU-filtered approximation.
    /// </summary>
    public static class DisplayReducer
    {
        /// <summary>Returns a top-down RGBA32 buffer of size ceil(W/n) x ceil(H/n).</summary>
        public static byte[] ReduceToRgba32(PixelImage img, int n, out int width, out int height)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
            int w = (img.Width + n - 1) / n;
            int h = (img.Height + n - 1) / n;
            width = w;
            height = h;
            var dst = new byte[(long)w * h * 4];
            if (n == 1)
            {
                img.CopyToRgba32(0, 0, img.Width, img.Height, dst, false);
                return dst;
            }
            int ch = img.Channels;
            bool alpha = img.HasAlpha;
            bool gray = img.IsGray;
            float[] lut = img.BitDepth == 8 ? Srgb.Decode8 : Srgb.Decode16;
            float aScale = 1f / img.MaxValue;
            int srcW = img.Width, srcH = img.Height;

            Parallel.For(0, h, oy =>
            {
                int y0 = oy * n, y1 = Math.Min(srcH, y0 + n);
                for (int ox = 0; ox < w; ox++)
                {
                    int x0 = ox * n, x1 = Math.Min(srcW, x0 + n);
                    float r = 0, g = 0, b = 0, a = 0;
                    int count = 0;
                    for (int y = y0; y < y1; y++)
                    {
                        int s = (y * srcW + x0) * ch;
                        for (int x = x0; x < x1; x++, s += ch)
                        {
                            float av = 1f;
                            if (img.BitDepth == 8)
                            {
                                var d = img.Data8;
                                if (gray) { r += lut[d[s]]; if (alpha) av = d[s + 1] * aScale; }
                                else { r += lut[d[s]]; g += lut[d[s + 1]]; b += lut[d[s + 2]]; if (alpha) av = d[s + 3] * aScale; }
                            }
                            else
                            {
                                var d = img.Data16;
                                if (gray) { r += lut[d[s]]; if (alpha) av = d[s + 1] * aScale; }
                                else { r += lut[d[s]]; g += lut[d[s + 1]]; b += lut[d[s + 2]]; if (alpha) av = d[s + 3] * aScale; }
                            }
                            a += av;
                            count++;
                        }
                    }
                    float inv = 1f / count;
                    int o = (oy * w + ox) * 4;
                    byte rr = (byte)(Srgb.Encode(r * inv) * 255f + 0.5f);
                    if (gray)
                    {
                        dst[o] = dst[o + 1] = dst[o + 2] = rr;
                    }
                    else
                    {
                        dst[o] = rr;
                        dst[o + 1] = (byte)(Srgb.Encode(g * inv) * 255f + 0.5f);
                        dst[o + 2] = (byte)(Srgb.Encode(b * inv) * 255f + 0.5f);
                    }
                    dst[o + 3] = (byte)(a * inv * 255f + 0.5f);
                }
            });
            return dst;
        }
    }
}
