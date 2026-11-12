using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Imaging;
using UnityEngine;

namespace FisheyeFriend.App.IO
{
    /// <summary>
    /// Loads JPEG and PNG files into an upright <see cref="PixelImage"/>.
    /// PNG uses the core decoder (16-bit preserved) on a worker thread. JPEG pixels are decoded with
    /// Unity's ImageConversion, which must run on the main thread, so <see cref="LoadAsync"/> must be
    /// awaited from the main thread. EXIF and ICC data are kept for writing back out.
    /// Images are turned upright from their EXIF orientation tag (Unity's decoder ignores it), the tag is
    /// reset to 1, and the original value is kept in Metadata.SourceOrientation for the squeeze-direction choice.
    /// </summary>
    internal static class ImageLoader
    {
        public static async Task<PixelImage> LoadAsync(string path, CancellationToken ct)
        {
            byte[] bytes = await Task.Run(() => File.ReadAllBytes(path), ct);
            return await DecodeAsync(bytes, ct);
        }

        public static async Task<PixelImage> DecodeAsync(byte[] bytes, CancellationToken ct)
        {
            var format = ImageIO.DetectFormat(bytes);
            if (format == ImageFileFormat.Png)
                return await Task.Run(() => ImageIO.MakeUpright(PngCodec.Decode(bytes)), ct);

            if (format == ImageFileFormat.Jpeg)
            {
                JpegInfo info = JpegInfo.Read(bytes);
                if (info.IsCmyk)
                    throw new NotSupportedException("CMYK JPEGs are not supported; convert the file to RGB first.");
                if (info.Width > SystemInfo.maxTextureSize || info.Height > SystemInfo.maxTextureSize)
                    throw new NotSupportedException("This JPEG (" + info.Width + " x " + info.Height +
                                                    ") is larger than the " + SystemInfo.maxTextureSize +
                                                    " px Unity can decode on this GPU.");
                int w, h;
                Color32[] px = DecodeJpegOnMainThread(bytes, out w, out h);
                ct.ThrowIfCancellationRequested();
                // Unity's decoder does not apply the orientation tag. Should a future version do so, the decoded size
                // comes back with width and height exchanged relative to the frame header; then only the tag is reset.
                bool decoderTurned = info.Width != info.Height && w == info.Height && h == info.Width;
                return await Task.Run(() => ImageIO.MakeUpright(FromColor32(px, w, h, info), decoderTurned), ct);
            }
            throw new InvalidDataException("Not a JPEG or PNG file.");
        }

        private static Color32[] DecodeJpegOnMainThread(byte[] bytes, out int w, out int h)
        {
            var tex = new Texture2D(2, 2, TextureFormat.RGBA32, false);
            try
            {
                if (!ImageConversion.LoadImage(tex, bytes, false))
                    throw new InvalidDataException("Unity could not decode this JPEG (arithmetic-coded or damaged?).");
                w = tex.width;
                h = tex.height;
                return tex.GetPixels32();
            }
            finally
            {
                if (Application.isPlaying) UnityEngine.Object.Destroy(tex);
                else UnityEngine.Object.DestroyImmediate(tex);
            }
        }

        private static PixelImage FromColor32(Color32[] px, int w, int h, JpegInfo info)
        {
            // Unity textures are bottom-up; PixelImage is top-down.
            var data = new byte[(long)w * h * 3];
            Parallel.For(0, h, y =>
            {
                int s = (h - 1 - y) * w;
                int d = y * w * 3;
                for (int x = 0; x < w; x++, s++, d += 3)
                {
                    Color32 c = px[s];
                    data[d] = c.r;
                    data[d + 1] = c.g;
                    data[d + 2] = c.b;
                }
            });
            return ImageIO.AttachJpegMetadata(new PixelImage(w, h, 3, data), info);
        }
    }
}
