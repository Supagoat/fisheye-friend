using System.Threading.Tasks;

namespace FisheyeFriend.Core.Imaging
{
    public static class Transforms
    {
        /// <summary>Swaps rows and columns (used to run horizontal algorithms vertically).</summary>
        public static PixelImage Transpose(PixelImage src)
        {
            int w = src.Width, h = src.Height, ch = src.Channels;
            var dst = new PixelImage(h, w, ch, src.BitDepth)
            {
                Metadata = src.Metadata != null ? src.Metadata.Clone() : new ImageMetadata(),
                SourceFormat = src.SourceFormat,
            };
            Parallel.For(0, h, y =>
            {
                for (int x = 0; x < w; x++)
                {
                    int s = (y * w + x) * ch, d = (x * h + y) * ch;
                    for (int c = 0; c < ch; c++)
                    {
                        if (src.BitDepth == 8) dst.Data8[d + c] = src.Data8[s + c];
                        else dst.Data16[d + c] = src.Data16[s + c];
                    }
                }
            });
            return dst;
        }
    }
}
