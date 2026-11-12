using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Text;
using System.Threading.Tasks;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Codecs
{
    public sealed class PngEncodeOptions
    {
        /// <summary>DEFLATE effort. PNG is lossless either way; Optimal only makes files smaller.</summary>
        public CompressionLevel Compression = CompressionLevel.Optimal;
    }

    /// <summary>
    /// Self-contained PNG reader/writer.
    /// Reads every standard PNG: gray, gray+alpha, RGB, RGBA and palette, bit depths 1-16, tRNS
    /// transparency and Adam7 interlacing. 16-bit images stay 16-bit; palette and low-bit-depth images
    /// become 8-bit. Writes 8- or 16-bit gray/gray+alpha/RGB/RGBA with adaptive row filters.
    /// ICC profiles (iCCP), EXIF (eXIf) and common ancillary chunks are preserved.
    /// </summary>
    public static class PngCodec
    {
        private static readonly byte[] Signature = { 0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A };

        /// <summary>Ancillary chunks copied from input to output unchanged.</summary>
        private static readonly HashSet<string> PassthroughChunks = new HashSet<string>
        {
            "sRGB", "gAMA", "cHRM", "pHYs", "tEXt", "zTXt", "iTXt", "sBIT", "cICP",
        };

        public static bool IsPng(byte[] data)
        {
            if (data == null || data.Length < 8) return false;
            for (int i = 0; i < 8; i++)
                if (data[i] != Signature[i]) return false;
            return true;
        }

        // =========================================================================================
        // Decoding
        // =========================================================================================

        public static PixelImage Decode(byte[] data)
        {
            if (!IsPng(data)) throw new InvalidDataException("Not a PNG file.");
            int pos = 8;
            int width = 0, height = 0, bitDepth = 0, colorType = -1, interlace = 0;
            byte[] palette = null, trns = null;
            var idat = new MemoryStream();
            var meta = new ImageMetadata();
            bool sawIhdr = false;

            while (pos + 8 <= data.Length)
            {
                int len = ReadInt32BE(data, pos);
                string type = Encoding.ASCII.GetString(data, pos + 4, 4);
                int dataStart = pos + 8;
                if (len < 0 || dataStart + len > data.Length) throw new InvalidDataException("PNG chunk " + type + " is truncated.");
                pos = dataStart + len + 4; // skip CRC

                switch (type)
                {
                    case "IHDR":
                        if (len < 13) throw new InvalidDataException("Bad IHDR chunk.");
                        width = ReadInt32BE(data, dataStart);
                        height = ReadInt32BE(data, dataStart + 4);
                        bitDepth = data[dataStart + 8];
                        colorType = data[dataStart + 9];
                        interlace = data[dataStart + 12];
                        sawIhdr = true;
                        break;
                    case "PLTE":
                        palette = Slice(data, dataStart, len);
                        break;
                    case "tRNS":
                        trns = Slice(data, dataStart, len);
                        break;
                    case "IDAT":
                        idat.Write(data, dataStart, len);
                        break;
                    case "iCCP":
                        meta.IccProfile = ReadIccp(data, dataStart, len);
                        break;
                    case "eXIf":
                        meta.Exif = Slice(data, dataStart, len);
                        break;
                    case "IEND":
                        pos = data.Length;
                        break;
                    default:
                        if (PassthroughChunks.Contains(type)) meta.PngChunks.Add(new PngChunk(type, Slice(data, dataStart, len)));
                        break;
                }
            }

            if (!sawIhdr) throw new InvalidDataException("PNG has no IHDR chunk.");
            if (width <= 0 || height <= 0) throw new InvalidDataException("PNG has invalid dimensions.");
            int samples = SamplesPerPixel(colorType);
            if (samples == 0) throw new InvalidDataException("Unsupported PNG colour type " + colorType + ".");
            if (!ValidDepth(colorType, bitDepth)) throw new InvalidDataException("Invalid PNG bit depth " + bitDepth + ".");
            if (colorType == 3 && palette == null) throw new InvalidDataException("Palette PNG has no PLTE chunk.");
            if (interlace > 1) throw new InvalidDataException("Unknown PNG interlace method.");

            int bitsPerPixel = samples * bitDepth;
            int bpp = Math.Max(1, bitsPerPixel / 8);

            // Output layout
            bool hasKey = trns != null && (colorType == 0 || colorType == 2);
            bool paletteAlpha = colorType == 3 && trns != null && trns.Length > 0;
            int outCh;
            switch (colorType)
            {
                case 0: outCh = hasKey ? 2 : 1; break;
                case 2: outCh = hasKey ? 4 : 3; break;
                case 3: outCh = paletteAlpha ? 4 : 3; break;
                case 4: outCh = 2; break;
                default: outCh = 4; break;
            }
            int outDepth = bitDepth == 16 ? 16 : 8;
            var img = new PixelImage(width, height, outCh, outDepth) { Metadata = meta, SourceFormat = ImageFileFormat.Png };

            byte[] raw;
            var idatBytes = idat.ToArray();
            var unpacker = new Unpacker(img, colorType, bitDepth, palette, trns);
            if (interlace == 0)
            {
                long rowBytes = ((long)width * bitsPerPixel + 7) / 8;
                raw = Zlib.Decompress(idatBytes, 0, idatBytes.Length, (rowBytes + 1) * height);
                Unfilter(raw, 0, (int)rowBytes, height, bpp);
                for (int y = 0; y < height; y++)
                    unpacker.UnpackRow(raw, y * (int)(rowBytes + 1) + 1, width, y, 0, 1);
            }
            else
            {
                long total = 0;
                for (int p = 0; p < 7; p++)
                {
                    int pw = PassSize(width, Adam7XStart[p], Adam7XStep[p]);
                    int ph = PassSize(height, Adam7YStart[p], Adam7YStep[p]);
                    if (pw == 0 || ph == 0) continue;
                    total += (((long)pw * bitsPerPixel + 7) / 8 + 1) * ph;
                }
                raw = Zlib.Decompress(idatBytes, 0, idatBytes.Length, total);
                int off = 0;
                for (int p = 0; p < 7; p++)
                {
                    int pw = PassSize(width, Adam7XStart[p], Adam7XStep[p]);
                    int ph = PassSize(height, Adam7YStart[p], Adam7YStep[p]);
                    if (pw == 0 || ph == 0) continue;
                    int rowBytes = (int)(((long)pw * bitsPerPixel + 7) / 8);
                    Unfilter(raw, off, rowBytes, ph, bpp);
                    for (int r = 0; r < ph; r++)
                        unpacker.UnpackRow(raw, off + r * (rowBytes + 1) + 1, pw,
                            Adam7YStart[p] + r * Adam7YStep[p], Adam7XStart[p], Adam7XStep[p]);
                    off += (rowBytes + 1) * ph;
                }
            }
            return img;
        }

        private static readonly int[] Adam7XStart = { 0, 4, 0, 2, 0, 1, 0 };
        private static readonly int[] Adam7YStart = { 0, 0, 4, 0, 2, 0, 1 };
        private static readonly int[] Adam7XStep = { 8, 8, 4, 4, 2, 2, 1 };
        private static readonly int[] Adam7YStep = { 8, 8, 8, 4, 4, 2, 2 };

        private static int PassSize(int size, int start, int step)
        {
            return size <= start ? 0 : (size - start + step - 1) / step;
        }

        private static int SamplesPerPixel(int colorType)
        {
            switch (colorType)
            {
                case 0: return 1;
                case 2: return 3;
                case 3: return 1;
                case 4: return 2;
                case 6: return 4;
                default: return 0;
            }
        }

        private static bool ValidDepth(int colorType, int d)
        {
            switch (colorType)
            {
                case 0: return d == 1 || d == 2 || d == 4 || d == 8 || d == 16;
                case 3: return d == 1 || d == 2 || d == 4 || d == 8;
                default: return d == 8 || d == 16;
            }
        }

        /// <summary>Reverses PNG row filters in place. Each row is preceded by its filter-type byte.</summary>
        private static void Unfilter(byte[] buf, int offset, int rowBytes, int rows, int bpp)
        {
            int stride = rowBytes + 1;
            for (int y = 0; y < rows; y++)
            {
                int row = offset + y * stride;
                int ft = buf[row];
                int cur = row + 1;
                int prev = y > 0 ? cur - stride : -1;
                switch (ft)
                {
                    case 0:
                        break;
                    case 1:
                        for (int i = bpp; i < rowBytes; i++) buf[cur + i] = (byte)(buf[cur + i] + buf[cur + i - bpp]);
                        break;
                    case 2:
                        if (prev >= 0)
                            for (int i = 0; i < rowBytes; i++) buf[cur + i] = (byte)(buf[cur + i] + buf[prev + i]);
                        break;
                    case 3:
                        for (int i = 0; i < rowBytes; i++)
                        {
                            int left = i >= bpp ? buf[cur + i - bpp] : 0;
                            int up = prev >= 0 ? buf[prev + i] : 0;
                            buf[cur + i] = (byte)(buf[cur + i] + ((left + up) >> 1));
                        }
                        break;
                    case 4:
                        for (int i = 0; i < rowBytes; i++)
                        {
                            int a = i >= bpp ? buf[cur + i - bpp] : 0;
                            int b = prev >= 0 ? buf[prev + i] : 0;
                            int c = (i >= bpp && prev >= 0) ? buf[prev + i - bpp] : 0;
                            buf[cur + i] = (byte)(buf[cur + i] + Paeth(a, b, c));
                        }
                        break;
                    default:
                        throw new InvalidDataException("Invalid PNG filter type " + ft + ".");
                }
            }
        }

        private static int Paeth(int a, int b, int c)
        {
            int p = a + b - c;
            int pa = Math.Abs(p - a), pb = Math.Abs(p - b), pc = Math.Abs(p - c);
            if (pa <= pb && pa <= pc) return a;
            return pb <= pc ? b : c;
        }

        /// <summary>Turns unfiltered PNG rows into interleaved 8/16-bit samples.</summary>
        private sealed class Unpacker
        {
            private readonly PixelImage _img;
            private readonly int _colorType;
            private readonly int _depth;
            private readonly byte[] _palette;
            private readonly byte[] _trns;
            private readonly int _keyGray = -1;
            private readonly int _keyR = -1, _keyG = -1, _keyB = -1;

            public Unpacker(PixelImage img, int colorType, int depth, byte[] palette, byte[] trns)
            {
                _img = img;
                _colorType = colorType;
                _depth = depth;
                _palette = palette;
                _trns = trns;
                if (trns != null && colorType == 0 && trns.Length >= 2) _keyGray = (trns[0] << 8) | trns[1];
                if (trns != null && colorType == 2 && trns.Length >= 6)
                {
                    _keyR = (trns[0] << 8) | trns[1];
                    _keyG = (trns[2] << 8) | trns[3];
                    _keyB = (trns[4] << 8) | trns[5];
                }
            }

            public void UnpackRow(byte[] raw, int off, int count, int y, int x0, int xStep)
            {
                int ch = _img.Channels;
                int width = _img.Width;
                if (_img.BitDepth == 16)
                {
                    var o = _img.Data16;
                    int spp = SamplesPerPixel(_colorType);
                    for (int i = 0; i < count; i++)
                    {
                        int x = x0 + i * xStep;
                        int d = (y * width + x) * ch;
                        int s = off + i * spp * 2;
                        for (int c = 0; c < spp; c++) o[d + c] = (ushort)((raw[s + c * 2] << 8) | raw[s + c * 2 + 1]);
                        if (_colorType == 0 && ch == 2) o[d + 1] = o[d] == _keyGray ? (ushort)0 : (ushort)65535;
                        if (_colorType == 2 && ch == 4)
                            o[d + 3] = (o[d] == _keyR && o[d + 1] == _keyG && o[d + 2] == _keyB) ? (ushort)0 : (ushort)65535;
                    }
                    return;
                }

                var dst = _img.Data8;
                if (_depth == 8)
                {
                    int spp = SamplesPerPixel(_colorType);
                    for (int i = 0; i < count; i++)
                    {
                        int x = x0 + i * xStep;
                        int d = (y * width + x) * ch;
                        int s = off + i * spp;
                        if (_colorType == 3)
                        {
                            WritePalette(dst, d, raw[s]);
                            continue;
                        }
                        for (int c = 0; c < spp; c++) dst[d + c] = raw[s + c];
                        if (_colorType == 0 && ch == 2) dst[d + 1] = raw[s] == _keyGray ? (byte)0 : (byte)255;
                        if (_colorType == 2 && ch == 4)
                            dst[d + 3] = (raw[s] == _keyR && raw[s + 1] == _keyG && raw[s + 2] == _keyB) ? (byte)0 : (byte)255;
                    }
                    return;
                }

                // 1, 2 or 4 bits per sample: gray or palette index, MSB first.
                int mask = (1 << _depth) - 1;
                int scale = _colorType == 0 ? 255 / mask : 1;
                for (int i = 0; i < count; i++)
                {
                    int bit = i * _depth;
                    int v = (raw[off + (bit >> 3)] >> (8 - _depth - (bit & 7))) & mask;
                    int x = x0 + i * xStep;
                    int d = (y * width + x) * ch;
                    if (_colorType == 3)
                    {
                        WritePalette(dst, d, v);
                    }
                    else
                    {
                        dst[d] = (byte)(v * scale);
                        if (ch == 2) dst[d + 1] = v == _keyGray ? (byte)0 : (byte)255;
                    }
                }
            }

            private void WritePalette(byte[] dst, int d, int index)
            {
                int p = index * 3;
                if (p + 2 < _palette.Length)
                {
                    dst[d] = _palette[p];
                    dst[d + 1] = _palette[p + 1];
                    dst[d + 2] = _palette[p + 2];
                }
                else
                {
                    dst[d] = dst[d + 1] = dst[d + 2] = 0;
                }
                if (_img.Channels == 4) dst[d + 3] = index < _trns.Length ? _trns[index] : (byte)255;
            }
        }

        private static byte[] ReadIccp(byte[] data, int start, int len)
        {
            int end = start + len;
            int p = start;
            while (p < end && data[p] != 0) p++; // profile name
            p += 2; // null separator + compression method
            if (p >= end) return null;
            try
            {
                return Zlib.Decompress(data, p, end - p, -1);
            }
            catch (Exception)
            {
                return null; // a broken profile should not stop the image from loading
            }
        }

        // =========================================================================================
        // Encoding
        // =========================================================================================

        public static byte[] Encode(PixelImage img, PngEncodeOptions options = null)
        {
            if (img == null) throw new ArgumentNullException(nameof(img));
            if (options == null) options = new PngEncodeOptions();
            int colorType;
            switch (img.Channels)
            {
                case 1: colorType = 0; break;
                case 2: colorType = 4; break;
                case 3: colorType = 2; break;
                default: colorType = 6; break;
            }
            int bytesPerSample = img.BitDepth / 8;
            int bpp = img.Channels * bytesPerSample;
            int rowBytes = img.Width * bpp;
            int stride = rowBytes + 1;
            long total = (long)stride * img.Height;
            if (total > int.MaxValue) throw new InvalidOperationException("Image is too large to encode as a single PNG buffer.");
            var filtered = new byte[total];

            Parallel.For(0, img.Height, () => new FilterScratch(rowBytes), (y, state, scratch) =>
            {
                WriteRawRow(img, y, scratch.Cur);
                if (y > 0) WriteRawRow(img, y - 1, scratch.Prev);
                FilterRow(scratch, y > 0, bpp, filtered, y * stride);
                return scratch;
            }, s => { });

            byte[] z = Zlib.Compress(filtered, 0, filtered.Length, options.Compression);

            using (var ms = new MemoryStream(z.Length + 1024))
            {
                ms.Write(Signature, 0, 8);
                var ihdr = new byte[13];
                WriteInt32BE(ihdr, 0, img.Width);
                WriteInt32BE(ihdr, 4, img.Height);
                ihdr[8] = (byte)img.BitDepth;
                ihdr[9] = (byte)colorType;
                ihdr[10] = 0; // deflate
                ihdr[11] = 0; // adaptive filtering
                ihdr[12] = 0; // no interlace
                WriteChunk(ms, "IHDR", ihdr);

                var meta = img.Metadata;
                bool hasIcc = meta != null && meta.IccProfile != null && meta.IccProfile.Length > 0;
                if (hasIcc) WriteChunk(ms, "iCCP", BuildIccp(meta.IccProfile));
                if (meta != null && meta.PngChunks != null)
                {
                    foreach (var c in meta.PngChunks)
                    {
                        if (hasIcc && (c.Type == "sRGB" || c.Type == "gAMA" || c.Type == "cHRM")) continue; // iCCP wins
                        if (c.Type == "sBIT") continue; // significant bits no longer describe resampled data
                        WriteChunk(ms, c.Type, c.Data);
                    }
                }
                if (meta != null && meta.Exif != null && meta.Exif.Length > 0) WriteChunk(ms, "eXIf", meta.Exif);

                const int MaxIdat = 1 << 20;
                for (int p = 0; p < z.Length; p += MaxIdat)
                {
                    int n = Math.Min(MaxIdat, z.Length - p);
                    var chunk = new byte[n];
                    Buffer.BlockCopy(z, p, chunk, 0, n);
                    WriteChunk(ms, "IDAT", chunk);
                }
                WriteChunk(ms, "IEND", new byte[0]);
                return ms.ToArray();
            }
        }

        private sealed class FilterScratch
        {
            public readonly byte[] Cur;
            public readonly byte[] Prev;
            public readonly byte[] Trial;

            public FilterScratch(int rowBytes)
            {
                Cur = new byte[rowBytes];
                Prev = new byte[rowBytes];
                Trial = new byte[rowBytes];
            }
        }

        private static void WriteRawRow(PixelImage img, int y, byte[] dst)
        {
            int n = img.Stride;
            int s = y * n;
            if (img.BitDepth == 8)
            {
                Buffer.BlockCopy(img.Data8, s, dst, 0, n);
            }
            else
            {
                var a = img.Data16;
                for (int i = 0; i < n; i++)
                {
                    ushort v = a[s + i];
                    dst[i * 2] = (byte)(v >> 8);
                    dst[i * 2 + 1] = (byte)v;
                }
            }
        }

        /// <summary>Picks the filter with the smallest sum of absolute signed residuals (the libpng heuristic).</summary>
        private static void FilterRow(FilterScratch s, bool hasPrev, int bpp, byte[] output, int outOffset)
        {
            byte[] cur = s.Cur, prev = s.Prev, trial = s.Trial;
            int n = cur.Length;
            int bestType = 0;
            long bestSum = long.MaxValue;

            for (int ft = 0; ft < 5; ft++)
            {
                if (!hasPrev && (ft == 2)) continue; // Up == None on the first row
                long sum = 0;
                for (int i = 0; i < n; i++)
                {
                    int a = i >= bpp ? cur[i - bpp] : 0;
                    int b = hasPrev ? prev[i] : 0;
                    int c = (hasPrev && i >= bpp) ? prev[i - bpp] : 0;
                    int pred;
                    switch (ft)
                    {
                        case 0: pred = 0; break;
                        case 1: pred = a; break;
                        case 2: pred = b; break;
                        case 3: pred = (a + b) >> 1; break;
                        default: pred = Paeth(a, b, c); break;
                    }
                    byte r = (byte)(cur[i] - pred);
                    trial[i] = r;
                    sum += r < 128 ? r : 256 - r;
                    if (sum >= bestSum) break;
                }
                if (sum < bestSum)
                {
                    bestSum = sum;
                    bestType = ft;
                    output[outOffset] = (byte)ft;
                    Buffer.BlockCopy(trial, 0, output, outOffset + 1, n);
                }
            }
            if (bestSum == long.MaxValue)
            {
                output[outOffset] = 0;
                Buffer.BlockCopy(cur, 0, output, outOffset + 1, n);
            }
        }

        private static byte[] BuildIccp(byte[] icc)
        {
            byte[] name = Encoding.ASCII.GetBytes("ICC Profile");
            byte[] z = Zlib.Compress(icc, 0, icc.Length, CompressionLevel.Optimal);
            var d = new byte[name.Length + 2 + z.Length];
            Buffer.BlockCopy(name, 0, d, 0, name.Length);
            d[name.Length] = 0;
            d[name.Length + 1] = 0; // compression method: deflate
            Buffer.BlockCopy(z, 0, d, name.Length + 2, z.Length);
            return d;
        }

        private static void WriteChunk(Stream s, string type, byte[] data)
        {
            var len = new byte[4];
            WriteInt32BE(len, 0, data.Length);
            s.Write(len, 0, 4);
            byte[] t = Encoding.ASCII.GetBytes(type);
            s.Write(t, 0, 4);
            s.Write(data, 0, data.Length);
            uint crc = Crc32.Compute(t, data);
            var c = new byte[4];
            WriteInt32BE(c, 0, (int)crc);
            s.Write(c, 0, 4);
        }

        private static int ReadInt32BE(byte[] b, int p)
        {
            return (b[p] << 24) | (b[p + 1] << 16) | (b[p + 2] << 8) | b[p + 3];
        }

        private static void WriteInt32BE(byte[] b, int p, int v)
        {
            b[p] = (byte)(v >> 24);
            b[p + 1] = (byte)(v >> 16);
            b[p + 2] = (byte)(v >> 8);
            b[p + 3] = (byte)v;
        }

        private static byte[] Slice(byte[] d, int start, int len)
        {
            var r = new byte[len];
            Buffer.BlockCopy(d, start, r, 0, len);
            return r;
        }
    }
}
