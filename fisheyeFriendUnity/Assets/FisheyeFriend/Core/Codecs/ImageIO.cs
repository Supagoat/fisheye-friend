using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Codecs
{
    /// <summary>Format detection and encoding to the output formats (PNG / max-quality JPEG).</summary>
    public static class ImageIO
    {
        public static readonly string[] SupportedExtensions = { ".jpg", ".jpeg", ".png" };

        public static ImageFileFormat DetectFormat(byte[] data)
        {
            if (PngCodec.IsPng(data)) return ImageFileFormat.Png;
            if (JpegInfo.IsJpeg(data)) return ImageFileFormat.Jpeg;
            return ImageFileFormat.Unknown;
        }

        public static ImageFileFormat FormatFromExtension(string path)
        {
            string ext = (Path.GetExtension(path) ?? "").ToLowerInvariant();
            if (ext == ".png") return ImageFileFormat.Png;
            if (ext == ".jpg" || ext == ".jpeg") return ImageFileFormat.Jpeg;
            return ImageFileFormat.Unknown;
        }

        public static bool IsSupportedPath(string path)
        {
            return FormatFromExtension(path) != ImageFileFormat.Unknown;
        }

        public static string DefaultExtension(ImageFileFormat f)
        {
            return f == ImageFileFormat.Png ? ".png" : ".jpg";
        }

        /// <summary>
        /// Turns a freshly decoded image upright according to its EXIF orientation tag (JPEG APP1 or PNG eXIf),
        /// resets the tag to 1 so the metadata matches the pixels, and records the original tag in
        /// <see cref="ImageMetadata.SourceOrientation"/>. Returns the upright image (the same instance when no
        /// turn was needed). Pass <paramref name="pixelsAlreadyUpright"/> when the decoder has already applied the
        /// tag; then only the metadata is updated.
        /// </summary>
        public static PixelImage MakeUpright(PixelImage img, bool pixelsAlreadyUpright = false)
        {
            if (img == null) throw new ArgumentNullException(nameof(img));
            if (img.Metadata == null) img.Metadata = new ImageMetadata();
            int o = img.Metadata.Exif != null ? ExifPatcher.ReadOrientation(img.Metadata.Exif) : 0;
            if (o <= Orientation.Normal)
            {
                img.Metadata.SourceOrientation = o;
                return img;
            }
            var upright = pixelsAlreadyUpright ? img : Orientation.Apply(img, o);
            var meta = upright.Metadata;
            meta.Exif = ExifPatcher.WithOrientation(meta.Exif, Orientation.Normal);
            meta.SourceOrientation = o;
            ResetXmpOrientation(meta.PngChunks);
            return upright;
        }

        /// <summary>Sets tiff:Orientation to 1 in uncompressed XMP packets carried in PNG iTXt chunks.</summary>
        private static void ResetXmpOrientation(List<PngChunk> chunks)
        {
            if (chunks == null) return;
            const string key = "XML:com.adobe.xmp";
            foreach (var c in chunks)
            {
                if (c == null || c.Type != "iTXt" || c.Data == null || c.Data.Length < key.Length + 5) continue;
                if (Encoding.ASCII.GetString(c.Data, 0, key.Length) != key || c.Data[key.Length] != 0) continue;
                int p = key.Length + 1;
                if (c.Data[p] != 0) continue; // compressed: leave alone
                p += 2;
                for (int n = 0; n < 2; n++) // language tag, translated keyword
                {
                    while (p < c.Data.Length && c.Data[p] != 0) p++;
                    p++;
                }
                if (p > c.Data.Length) continue;
                string xmp = Encoding.UTF8.GetString(c.Data, p, c.Data.Length - p);
                string fixedXmp = Regex.Replace(xmp, "(tiff:Orientation\\s*=\\s*[\"'])[2-8]([\"'])", "${1}1${2}");
                fixedXmp = Regex.Replace(fixedXmp, "(<tiff:Orientation>\\s*)[2-8](\\s*</tiff:Orientation>)", "${1}1${2}");
                if (fixedXmp == xmp) continue;
                byte[] body = Encoding.UTF8.GetBytes(fixedXmp);
                var d = new byte[p + body.Length];
                Buffer.BlockCopy(c.Data, 0, d, 0, p);
                Buffer.BlockCopy(body, 0, d, p, body.Length);
                c.Data = d;
            }
        }

        /// <summary>
        /// Encodes an image. JPEG is always written at maximum quality (Q100, 4:4:4, optimised Huffman).
        /// EXIF pixel-dimension tags are updated to the image's size and the embedded EXIF thumbnail (a stale,
        /// squeezed copy) is unlinked; ICC/EXIF/ancillary data is otherwise preserved.
        /// </summary>
        public static byte[] Encode(PixelImage img, ImageFileFormat format)
        {
            if (img == null) throw new ArgumentNullException(nameof(img));
            var original = img.Metadata;
            if (original != null && original.Exif != null)
            {
                // Encode from a shallow wrapper so the caller's metadata is not modified.
                var meta = original.Clone();
                meta.Exif = ExifPatcher.WithoutThumbnail(ExifPatcher.WithDimensions(meta.Exif, img.Width, img.Height));
                img.Metadata = meta;
                try
                {
                    return EncodeRaw(img, format);
                }
                finally
                {
                    img.Metadata = original;
                }
            }
            return EncodeRaw(img, format);
        }

        private static byte[] EncodeRaw(PixelImage img, ImageFileFormat format)
        {
            switch (format)
            {
                case ImageFileFormat.Png:
                    return PngCodec.Encode(img);
                case ImageFileFormat.Jpeg:
                    return JpegEncoder.Encode(img, new JpegEncodeOptions { Quality = 100 });
                default:
                    throw new ArgumentException("Unsupported output format " + format + ".");
            }
        }

        /// <summary>
        /// Completes a JPEG decoded by an external decoder (which returns RGB): attaches EXIF/ICC from
        /// the file and collapses grayscale JPEGs back to one channel so they are saved as grayscale.
        /// </summary>
        public static PixelImage AttachJpegMetadata(PixelImage decodedRgb, JpegInfo info)
        {
            if (decodedRgb == null) throw new ArgumentNullException(nameof(decodedRgb));
            PixelImage img = decodedRgb;
            if (info != null && info.IsGrayscale && decodedRgb.Channels >= 3 && decodedRgb.BitDepth == 8)
            {
                int n = decodedRgb.Width * decodedRgb.Height;
                var g = new byte[n];
                var s = decodedRgb.Data8;
                int ch = decodedRgb.Channels;
                for (int i = 0; i < n; i++) g[i] = s[i * ch];
                img = new PixelImage(decodedRgb.Width, decodedRgb.Height, 1, g);
            }
            img.SourceFormat = ImageFileFormat.Jpeg;
            img.Metadata = new ImageMetadata();
            if (info != null)
            {
                img.Metadata.Exif = info.Exif;
                img.Metadata.IccProfile = info.IccProfile;
            }
            return img;
        }
    }
}
