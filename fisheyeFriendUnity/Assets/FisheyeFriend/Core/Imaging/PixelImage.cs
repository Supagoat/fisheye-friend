using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace FisheyeFriend.Core.Imaging
{
    public enum ImageFileFormat
    {
        Unknown = 0,
        Png = 1,
        Jpeg = 2,
    }

    /// <summary>
    /// An interleaved, top-down raster with 1-4 channels (gray, gray+alpha, RGB, RGBA)
    /// at 8 or 16 bits per sample. Exactly one of <see cref="Data8"/> / <see cref="Data16"/> is set.
    /// </summary>
    public sealed class PixelImage
    {
        public int Width { get; }
        public int Height { get; }
        public int Channels { get; }
        public int BitDepth { get; }
        public byte[] Data8 { get; }
        public ushort[] Data16 { get; }

        /// <summary>Metadata carried from the source file (ICC profile, EXIF, PNG ancillary chunks).</summary>
        public ImageMetadata Metadata { get; set; } = new ImageMetadata();

        /// <summary>The file format this image was decoded from, if any.</summary>
        public ImageFileFormat SourceFormat { get; set; }

        public PixelImage(int width, int height, int channels, int bitDepth)
        {
            Validate(width, height, channels, bitDepth);
            Width = width;
            Height = height;
            Channels = channels;
            BitDepth = bitDepth;
            long n = (long)width * height * channels;
            if (bitDepth == 8) Data8 = new byte[n];
            else Data16 = new ushort[n];
        }

        public PixelImage(int width, int height, int channels, byte[] data)
        {
            Validate(width, height, channels, 8);
            if (data == null || data.LongLength != (long)width * height * channels)
                throw new ArgumentException("Pixel buffer size does not match the dimensions.");
            Width = width;
            Height = height;
            Channels = channels;
            BitDepth = 8;
            Data8 = data;
        }

        public PixelImage(int width, int height, int channels, ushort[] data)
        {
            Validate(width, height, channels, 16);
            if (data == null || data.LongLength != (long)width * height * channels)
                throw new ArgumentException("Pixel buffer size does not match the dimensions.");
            Width = width;
            Height = height;
            Channels = channels;
            BitDepth = 16;
            Data16 = data;
        }

        private static void Validate(int width, int height, int channels, int bitDepth)
        {
            if (width <= 0 || height <= 0) throw new ArgumentException("Image dimensions must be positive.");
            if (channels < 1 || channels > 4) throw new ArgumentException("Channels must be 1-4.");
            if (bitDepth != 8 && bitDepth != 16) throw new ArgumentException("Bit depth must be 8 or 16.");
            if ((long)width * height * channels > int.MaxValue)
                throw new ArgumentException("Image is too large (" + width + " x " + height + ").");
        }

        public bool HasAlpha => Channels == 2 || Channels == 4;
        public bool IsGray => Channels <= 2;
        public int ColorChannels => HasAlpha ? Channels - 1 : Channels;
        /// <summary>Samples per row.</summary>
        public int Stride => Width * Channels;
        public int MaxValue => BitDepth == 16 ? 65535 : 255;

        public string Describe()
        {
            string model;
            switch (Channels)
            {
                case 1: model = "Gray"; break;
                case 2: model = "Gray+Alpha"; break;
                case 3: model = "RGB"; break;
                default: model = "RGBA"; break;
            }
            return Width + " x " + Height + " " + model + " " + BitDepth + "-bit";
        }

        /// <summary>A new, blank image with the same channel layout, depth, format and a copy of the metadata.</summary>
        public PixelImage CreateLike(int width, int height)
        {
            return CreateLike(width, height, BitDepth);
        }

        /// <summary>A new, blank image with the same channel layout, format and a copy of the metadata, at <paramref name="bitDepth"/> (8 or 16).</summary>
        public PixelImage CreateLike(int width, int height, int bitDepth)
        {
            var img = new PixelImage(width, height, Channels, bitDepth);
            img.Metadata = Metadata != null ? Metadata.Clone() : new ImageMetadata();
            img.SourceFormat = SourceFormat;
            return img;
        }

        /// <summary>Returns this image at 8 bits per sample (rounded), or itself if already 8-bit.</summary>
        public PixelImage To8Bit()
        {
            if (BitDepth == 8) return this;
            var src = Data16;
            var dst = new byte[src.Length];
            for (int i = 0; i < src.Length; i++)
                dst[i] = (byte)((src[i] * 255 + 32767) / 65535);
            var img = new PixelImage(Width, Height, Channels, dst);
            img.Metadata = Metadata != null ? Metadata.Clone() : new ImageMetadata();
            img.SourceFormat = SourceFormat;
            return img;
        }

        /// <summary>
        /// Writes a top-down RGBA32 copy of rows [y0, y0+rows) and columns [x0, x0+cols) into dst
        /// (cols*rows*4 bytes). Gray is expanded, missing alpha is opaque, 16-bit is rounded to 8-bit.
        /// Used for display.
        /// </summary>
        public void CopyToRgba32(int x0, int y0, int cols, int rows, byte[] dst, bool flipVertically)
        {
            if (x0 < 0 || y0 < 0 || cols <= 0 || rows <= 0 || x0 + cols > Width || y0 + rows > Height)
                throw new ArgumentOutOfRangeException("Region is outside the image.");
            if (dst.Length < cols * rows * 4) throw new ArgumentException("Destination is too small.");
            int ch = Channels;
            Parallel.For(0, rows, r =>
            {
                int sy = y0 + r;
                int dy = flipVertically ? rows - 1 - r : r;
                int d = dy * cols * 4;
                int s = (sy * Width + x0) * ch;
                if (BitDepth == 8)
                {
                    var a = Data8;
                    for (int x = 0; x < cols; x++, s += ch, d += 4)
                    {
                        switch (ch)
                        {
                            case 1: dst[d] = dst[d + 1] = dst[d + 2] = a[s]; dst[d + 3] = 255; break;
                            case 2: dst[d] = dst[d + 1] = dst[d + 2] = a[s]; dst[d + 3] = a[s + 1]; break;
                            case 3: dst[d] = a[s]; dst[d + 1] = a[s + 1]; dst[d + 2] = a[s + 2]; dst[d + 3] = 255; break;
                            default: dst[d] = a[s]; dst[d + 1] = a[s + 1]; dst[d + 2] = a[s + 2]; dst[d + 3] = a[s + 3]; break;
                        }
                    }
                }
                else
                {
                    var a = Data16;
                    for (int x = 0; x < cols; x++, s += ch, d += 4)
                    {
                        switch (ch)
                        {
                            case 1: dst[d] = dst[d + 1] = dst[d + 2] = To8(a[s]); dst[d + 3] = 255; break;
                            case 2: dst[d] = dst[d + 1] = dst[d + 2] = To8(a[s]); dst[d + 3] = To8(a[s + 1]); break;
                            case 3: dst[d] = To8(a[s]); dst[d + 1] = To8(a[s + 1]); dst[d + 2] = To8(a[s + 2]); dst[d + 3] = 255; break;
                            default: dst[d] = To8(a[s]); dst[d + 1] = To8(a[s + 1]); dst[d + 2] = To8(a[s + 2]); dst[d + 3] = To8(a[s + 3]); break;
                        }
                    }
                }
            });
        }

        /// <summary>Reads one pixel as 0..MaxValue samples (for the cursor readout).</summary>
        public void GetPixel(int x, int y, int[] samples)
        {
            int s = (y * Width + x) * Channels;
            for (int c = 0; c < Channels; c++)
                samples[c] = BitDepth == 8 ? Data8[s + c] : Data16[s + c];
        }

        private static byte To8(ushort v)
        {
            return (byte)((v * 255 + 32767) / 65535);
        }
    }

    /// <summary>A PNG ancillary chunk carried through unchanged (text, pHYs, gAMA, ...).</summary>
    public sealed class PngChunk
    {
        public string Type;
        public byte[] Data;

        public PngChunk(string type, byte[] data)
        {
            Type = type;
            Data = data;
        }
    }

    /// <summary>Metadata preserved from the source file into the output.</summary>
    public sealed class ImageMetadata
    {
        /// <summary>Raw ICC colour profile, if the source had one.</summary>
        public byte[] IccProfile;

        /// <summary>EXIF block as a TIFF structure (starts with "II*\0" or "MM\0*"), without the "Exif\0\0" prefix.</summary>
        public byte[] Exif;

        /// <summary>PNG ancillary chunks carried through when writing PNG (sRGB, gAMA, cHRM, pHYs, text).</summary>
        public List<PngChunk> PngChunks = new List<PngChunk>();

        /// <summary>
        /// The EXIF orientation (1-8) the source file was tagged with, or 0 if it had none. Set by
        /// <see cref="Codecs.ImageIO.MakeUpright"/>, which turns the pixels upright and resets the tag in
        /// <see cref="Exif"/> to 1; this field remembers how the camera was held.
        /// </summary>
        public int SourceOrientation;

        public ImageMetadata Clone()
        {
            var m = new ImageMetadata
            {
                IccProfile = IccProfile != null ? (byte[])IccProfile.Clone() : null,
                Exif = Exif != null ? (byte[])Exif.Clone() : null,
                SourceOrientation = SourceOrientation,
            };
            if (PngChunks != null)
                foreach (var c in PngChunks)
                    m.PngChunks.Add(new PngChunk(c.Type, (byte[])c.Data.Clone()));
            return m;
        }
    }
}
