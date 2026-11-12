using System;
using System.IO;
using System.Threading.Tasks;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Codecs
{
    public sealed class JpegEncodeOptions
    {
        /// <summary>IJG quality 1-100. 100 makes every quantisation step 1 (the most a JPEG can keep).</summary>
        public int Quality = 100;
    }

    /// <summary>
    /// Baseline JPEG encoder built for maximum quality:
    ///  * quality 100 by default (all quantisation steps = 1),
    ///  * 4:4:4 - no chroma subsampling,
    ///  * colour conversion and DCT in floating point straight from the source samples
    ///    (16-bit input is not rounded to 8-bit before the transform),
    ///  * Huffman tables optimised for each image (smaller files, identical pixels).
    /// Grayscale images are written as single-component JPEGs. Alpha is flattened onto white.
    /// EXIF and ICC profiles from the metadata are written as APP1/APP2 segments.
    /// </summary>
    public static class JpegEncoder
    {
        private static readonly int[] ZigZag =
        {
            0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18, 11, 4, 5,
            12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28,
            35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
            58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
        };

        private static readonly int[] StdLuminanceQ =
        {
            16, 11, 10, 16, 24, 40, 51, 61,
            12, 12, 14, 19, 26, 58, 60, 55,
            14, 13, 16, 24, 40, 57, 69, 56,
            14, 17, 22, 29, 51, 87, 80, 62,
            18, 22, 37, 56, 68, 109, 103, 77,
            24, 35, 55, 64, 81, 104, 113, 92,
            49, 64, 78, 87, 103, 121, 120, 101,
            72, 92, 95, 98, 112, 100, 103, 99,
        };

        private static readonly int[] StdChrominanceQ =
        {
            17, 18, 24, 47, 99, 99, 99, 99,
            18, 21, 26, 66, 99, 99, 99, 99,
            24, 26, 56, 99, 99, 99, 99, 99,
            47, 66, 99, 99, 99, 99, 99, 99,
            99, 99, 99, 99, 99, 99, 99, 99,
            99, 99, 99, 99, 99, 99, 99, 99,
            99, 99, 99, 99, 99, 99, 99, 99,
            99, 99, 99, 99, 99, 99, 99, 99,
        };

        private static readonly float[] AanScale =
        {
            1.0f, 1.387039845f, 1.306562965f, 1.175875602f, 1.0f, 0.785694958f, 0.541196100f, 0.275899379f,
        };

        public static byte[] Encode(PixelImage img, JpegEncodeOptions options = null)
        {
            if (img == null) throw new ArgumentNullException(nameof(img));
            if (options == null) options = new JpegEncodeOptions();
            if (img.Width > 65535 || img.Height > 65535) throw new ArgumentException("JPEG is limited to 65535 x 65535 pixels.");
            int quality = Math.Max(1, Math.Min(100, options.Quality));

            bool gray = img.IsGray;
            int ncomp = gray ? 1 : 3;
            int[] qLuma = ScaleQuant(StdLuminanceQ, quality);
            int[] qChroma = ScaleQuant(StdChrominanceQ, quality);
            float[] divLuma = Divisors(qLuma);
            float[] divChroma = Divisors(qChroma);

            int bw = (img.Width + 7) / 8;
            int bh = (img.Height + 7) / 8;
            int blocks = bw * bh;
            var coef = new short[ncomp][];
            for (int c = 0; c < ncomp; c++) coef[c] = new short[(long)blocks * 64];

            // ---- 1. colour conversion + FDCT + quantisation, in parallel per block row
            Parallel.For(0, bh, () => new float[3 * 64], (by, state, buf) =>
            {
                var block = new float[64];
                for (int bx = 0; bx < bw; bx++)
                {
                    LoadBlock(img, bx * 8, by * 8, buf, gray);
                    int bi = (by * bw + bx) * 64;
                    for (int c = 0; c < ncomp; c++)
                    {
                        Array.Copy(buf, c * 64, block, 0, 64);
                        ForwardDct(block);
                        float[] div = c == 0 ? divLuma : divChroma;
                        short[] dst = coef[c];
                        for (int k = 0; k < 64; k++)
                        {
                            int n = ZigZag[k];
                            float v = block[n] * div[n];
                            int q = (int)Math.Floor(v + 0.5f);
                            int lim = k == 0 ? 2047 : 1023;
                            if (q > lim) q = lim;
                            if (q < -lim) q = -lim;
                            dst[bi + k] = (short)q;
                        }
                    }
                }
                return buf;
            }, b => { });

            // ---- 2. symbol statistics -> optimal Huffman tables
            var dcFreq = new long[2][];
            var acFreq = new long[2][];
            for (int t = 0; t < 2; t++)
            {
                dcFreq[t] = new long[257];
                acFreq[t] = new long[257];
            }
            for (int c = 0; c < ncomp; c++)
            {
                int t = c == 0 ? 0 : 1;
                GatherStats(coef[c], blocks, dcFreq[t], acFreq[t]);
            }
            int tables = ncomp == 1 ? 1 : 2;
            var dc = new HuffTable[tables];
            var ac = new HuffTable[tables];
            for (int t = 0; t < tables; t++)
            {
                dc[t] = HuffTable.Optimal(dcFreq[t]);
                ac[t] = HuffTable.Optimal(acFreq[t]);
            }

            // ---- 3. write the file
            var w = new JpegWriter(EstimateSize(img));
            w.Marker(0xD8); // SOI
            var meta = img.Metadata;
            byte[] exif = meta != null ? meta.Exif : null;
            bool wroteExif = false;
            if (exif != null && exif.Length > 0 && exif.Length <= 65533 - 6)
            {
                var seg = new byte[6 + exif.Length];
                seg[0] = (byte)'E'; seg[1] = (byte)'x'; seg[2] = (byte)'i'; seg[3] = (byte)'f'; seg[4] = 0; seg[5] = 0;
                Buffer.BlockCopy(exif, 0, seg, 6, exif.Length);
                w.Segment(0xE1, seg);
                wroteExif = true;
            }
            if (!wroteExif)
            {
                // JFIF: version 1.01, no units, 1:1 pixel aspect, no thumbnail.
                w.Segment(0xE0, new byte[] { (byte)'J', (byte)'F', (byte)'I', (byte)'F', 0, 1, 1, 0, 0, 1, 0, 1, 0, 0 });
            }
            if (meta != null && meta.IccProfile != null && meta.IccProfile.Length > 0)
                WriteIcc(w, meta.IccProfile);

            // DQT
            var dqt = new byte[tables * 65];
            for (int t = 0; t < tables; t++)
            {
                int[] q = t == 0 ? qLuma : qChroma;
                dqt[t * 65] = (byte)t; // 8-bit precision, table id t
                for (int k = 0; k < 64; k++) dqt[t * 65 + 1 + k] = (byte)q[ZigZag[k]];
            }
            w.Segment(0xDB, dqt);

            // SOF0 (baseline)
            var sof = new byte[6 + ncomp * 3];
            sof[0] = 8;
            sof[1] = (byte)(img.Height >> 8);
            sof[2] = (byte)img.Height;
            sof[3] = (byte)(img.Width >> 8);
            sof[4] = (byte)img.Width;
            sof[5] = (byte)ncomp;
            for (int c = 0; c < ncomp; c++)
            {
                sof[6 + c * 3] = (byte)(c + 1);
                sof[7 + c * 3] = 0x11; // 1x1 sampling: no chroma subsampling
                sof[8 + c * 3] = (byte)(c == 0 ? 0 : 1);
            }
            w.Segment(0xC0, sof);

            // DHT
            for (int t = 0; t < tables; t++)
            {
                w.Segment(0xC4, dc[t].Serialize((byte)(0x00 | t)));
                w.Segment(0xC4, ac[t].Serialize((byte)(0x10 | t)));
            }

            // SOS
            var sos = new byte[4 + ncomp * 2];
            sos[0] = (byte)ncomp;
            for (int c = 0; c < ncomp; c++)
            {
                sos[1 + c * 2] = (byte)(c + 1);
                sos[2 + c * 2] = (byte)(c == 0 ? 0x00 : 0x11);
            }
            sos[1 + ncomp * 2] = 0;   // Ss
            sos[2 + ncomp * 2] = 63;  // Se
            sos[3 + ncomp * 2] = 0;   // Ah/Al
            w.Segment(0xDA, sos);

            // Entropy-coded data, MCU = one 8x8 block per component.
            var lastDc = new int[ncomp];
            for (int b = 0; b < blocks; b++)
            {
                for (int c = 0; c < ncomp; c++)
                {
                    int t = c == 0 ? 0 : 1;
                    EncodeBlock(w, coef[c], b * 64, ref lastDc[c], dc[t], ac[t]);
                }
            }
            w.FlushBits();
            w.Marker(0xD9); // EOI
            return w.ToArray();
        }

        private static int EstimateSize(PixelImage img)
        {
            long est = (long)img.Width * img.Height * (img.IsGray ? 1 : 3) / 2 + 65536;
            return (int)Math.Min(est, int.MaxValue / 2);
        }

        private static int[] ScaleQuant(int[] baseTable, int quality)
        {
            int scale = quality < 50 ? 5000 / quality : 200 - quality * 2;
            var q = new int[64];
            for (int i = 0; i < 64; i++)
            {
                int v = (baseTable[i] * scale + 50) / 100;
                q[i] = v < 1 ? 1 : (v > 255 ? 255 : v);
            }
            return q;
        }

        /// <summary>Reciprocal divisors that also undo the AAN DCT's output scaling.</summary>
        private static float[] Divisors(int[] q)
        {
            var d = new float[64];
            for (int row = 0; row < 8; row++)
                for (int col = 0; col < 8; col++)
                    d[row * 8 + col] = (float)(1.0 / (q[row * 8 + col] * AanScale[row] * AanScale[col] * 8.0));
            return d;
        }

        /// <summary>
        /// Fills buf with level-shifted component samples for the 8x8 block at (x0,y0):
        /// Y at [0..64), Cb at [64..128), Cr at [128..192) (or gray at [0..64)).
        /// Edges are padded by replication. Alpha is flattened onto white.
        /// </summary>
        private static void LoadBlock(PixelImage img, int x0, int y0, float[] buf, bool gray)
        {
            int w = img.Width, h = img.Height, ch = img.Channels;
            bool alpha = img.HasAlpha;
            float scale = 255f / img.MaxValue;
            float aScale = 1f / img.MaxValue;
            for (int yy = 0; yy < 8; yy++)
            {
                int sy = Math.Min(y0 + yy, h - 1);
                for (int xx = 0; xx < 8; xx++)
                {
                    int sx = Math.Min(x0 + xx, w - 1);
                    int s = (sy * w + sx) * ch;
                    float r, g, b;
                    float a = 1f;
                    if (img.BitDepth == 8)
                    {
                        var d = img.Data8;
                        if (gray)
                        {
                            r = g = b = d[s];
                            if (alpha) a = d[s + 1] * aScale;
                        }
                        else
                        {
                            r = d[s];
                            g = d[s + 1];
                            b = d[s + 2];
                            if (alpha) a = d[s + 3] * aScale;
                        }
                    }
                    else
                    {
                        var d = img.Data16;
                        if (gray)
                        {
                            r = g = b = d[s] * scale;
                            if (alpha) a = d[s + 1] * aScale;
                        }
                        else
                        {
                            r = d[s] * scale;
                            g = d[s + 1] * scale;
                            b = d[s + 2] * scale;
                            if (alpha) a = d[s + 3] * aScale;
                        }
                    }
                    if (alpha)
                    {
                        float bg = 255f * (1f - a);
                        r = r * a + bg;
                        g = g * a + bg;
                        b = b * a + bg;
                    }
                    int i = yy * 8 + xx;
                    if (gray)
                    {
                        buf[i] = r - 128f;
                    }
                    else
                    {
                        // JFIF full-range YCbCr, kept in float.
                        float y = 0.299f * r + 0.587f * g + 0.114f * b;
                        float cb = -0.168735892f * r - 0.331264108f * g + 0.5f * b;
                        float cr = 0.5f * r - 0.418687589f * g - 0.081312411f * b;
                        buf[i] = Clamp(y, 0f, 255f) - 128f;
                        buf[64 + i] = Clamp(cb, -128f, 127f);
                        buf[128 + i] = Clamp(cr, -128f, 127f);
                    }
                }
            }
        }

        private static float Clamp(float v, float lo, float hi)
        {
            return v < lo ? lo : (v > hi ? hi : v);
        }

        /// <summary>Arai-Agui-Nakajima floating-point forward DCT (as libjpeg's jfdctflt.c). Output is scaled; see Divisors.</summary>
        private static void ForwardDct(float[] d)
        {
            for (int pass = 0; pass < 2; pass++)
            {
                int step = pass == 0 ? 1 : 8;   // element step within a line
                int next = pass == 0 ? 8 : 1;   // step between lines
                for (int line = 0; line < 8; line++)
                {
                    int b = line * next;
                    float d0 = d[b], d1 = d[b + step], d2 = d[b + 2 * step], d3 = d[b + 3 * step];
                    float d4 = d[b + 4 * step], d5 = d[b + 5 * step], d6 = d[b + 6 * step], d7 = d[b + 7 * step];

                    float tmp0 = d0 + d7, tmp7 = d0 - d7;
                    float tmp1 = d1 + d6, tmp6 = d1 - d6;
                    float tmp2 = d2 + d5, tmp5 = d2 - d5;
                    float tmp3 = d3 + d4, tmp4 = d3 - d4;

                    float tmp10 = tmp0 + tmp3, tmp13 = tmp0 - tmp3;
                    float tmp11 = tmp1 + tmp2, tmp12 = tmp1 - tmp2;
                    d[b] = tmp10 + tmp11;
                    d[b + 4 * step] = tmp10 - tmp11;
                    float z1 = (tmp12 + tmp13) * 0.707106781f;
                    d[b + 2 * step] = tmp13 + z1;
                    d[b + 6 * step] = tmp13 - z1;

                    tmp10 = tmp4 + tmp5;
                    tmp11 = tmp5 + tmp6;
                    tmp12 = tmp6 + tmp7;
                    float z5 = (tmp10 - tmp12) * 0.382683433f;
                    float z2 = 0.541196100f * tmp10 + z5;
                    float z4 = 1.306562965f * tmp12 + z5;
                    float z3 = tmp11 * 0.707106781f;
                    float z11 = tmp7 + z3, z13 = tmp7 - z3;
                    d[b + 5 * step] = z13 + z2;
                    d[b + 3 * step] = z13 - z2;
                    d[b + 1 * step] = z11 + z4;
                    d[b + 7 * step] = z11 - z4;
                }
            }
        }

        private static int BitLength(int v)
        {
            if (v < 0) v = -v;
            int n = 0;
            while (v != 0)
            {
                n++;
                v >>= 1;
            }
            return n;
        }

        private static void GatherStats(short[] coef, int blocks, long[] dcFreq, long[] acFreq)
        {
            int last = 0;
            for (int b = 0; b < blocks; b++)
            {
                int o = b * 64;
                int diff = coef[o] - last;
                last = coef[o];
                dcFreq[BitLength(diff)]++;
                int run = 0;
                for (int k = 1; k < 64; k++)
                {
                    int v = coef[o + k];
                    if (v == 0)
                    {
                        run++;
                        continue;
                    }
                    while (run > 15)
                    {
                        acFreq[0xF0]++;
                        run -= 16;
                    }
                    acFreq[(run << 4) | BitLength(v)]++;
                    run = 0;
                }
                if (run > 0) acFreq[0x00]++;
            }
        }

        private static void EncodeBlock(JpegWriter w, short[] coef, int o, ref int lastDc, HuffTable dc, HuffTable ac)
        {
            int diff = coef[o] - lastDc;
            lastDc = coef[o];
            int nb = BitLength(diff);
            w.WriteBits(dc.Code[nb], dc.Size[nb]);
            if (nb > 0) w.WriteBits((diff < 0 ? diff - 1 : diff) & ((1 << nb) - 1), nb);

            int run = 0;
            for (int k = 1; k < 64; k++)
            {
                int v = coef[o + k];
                if (v == 0)
                {
                    run++;
                    continue;
                }
                while (run > 15)
                {
                    w.WriteBits(ac.Code[0xF0], ac.Size[0xF0]);
                    run -= 16;
                }
                nb = BitLength(v);
                int sym = (run << 4) | nb;
                w.WriteBits(ac.Code[sym], ac.Size[sym]);
                w.WriteBits((v < 0 ? v - 1 : v) & ((1 << nb) - 1), nb);
                run = 0;
            }
            if (run > 0) w.WriteBits(ac.Code[0x00], ac.Size[0x00]);
        }

        private static void WriteIcc(JpegWriter w, byte[] icc)
        {
            const int MaxChunk = 65533 - 14;
            int count = (icc.Length + MaxChunk - 1) / MaxChunk;
            if (count > 255) return; // absurdly large profile; skip rather than write a broken file
            for (int i = 0; i < count; i++)
            {
                int start = i * MaxChunk;
                int n = Math.Min(MaxChunk, icc.Length - start);
                var seg = new byte[14 + n];
                byte[] tag = System.Text.Encoding.ASCII.GetBytes("ICC_PROFILE");
                Buffer.BlockCopy(tag, 0, seg, 0, tag.Length);
                seg[11] = 0;
                seg[12] = (byte)(i + 1);
                seg[13] = (byte)count;
                Buffer.BlockCopy(icc, start, seg, 14, n);
                w.Segment(0xE2, seg);
            }
        }

        // -----------------------------------------------------------------------------------------

        private sealed class HuffTable
        {
            public readonly byte[] Bits = new byte[17]; // Bits[1..16]: number of codes of each length
            public byte[] Values;
            public readonly int[] Code = new int[256];
            public readonly int[] Size = new int[256];

            /// <summary>Optimal length-limited Huffman table (ITU T.81 Annex K.2, as libjpeg's jpeg_gen_optimal_table).</summary>
            public static HuffTable Optimal(long[] freqIn)
            {
                const int MaxCodeLen = 32;
                var freq = new long[257];
                Array.Copy(freqIn, freq, 256);
                bool any = false;
                for (int i = 0; i < 256; i++) any |= freq[i] != 0;
                if (!any) freq[0] = 1;
                freq[256] = 1; // reserved so no real code is all ones

                var codeSize = new int[257];
                var others = new int[257];
                for (int i = 0; i < 257; i++) others[i] = -1;

                while (true)
                {
                    int c1 = -1;
                    long v = long.MaxValue;
                    for (int i = 0; i <= 256; i++)
                        if (freq[i] != 0 && freq[i] <= v) { v = freq[i]; c1 = i; }
                    int c2 = -1;
                    v = long.MaxValue;
                    for (int i = 0; i <= 256; i++)
                        if (freq[i] != 0 && freq[i] <= v && i != c1) { v = freq[i]; c2 = i; }
                    if (c2 < 0) break;

                    freq[c1] += freq[c2];
                    freq[c2] = 0;
                    codeSize[c1]++;
                    while (others[c1] >= 0)
                    {
                        c1 = others[c1];
                        codeSize[c1]++;
                    }
                    others[c1] = c2;
                    codeSize[c2]++;
                    while (others[c2] >= 0)
                    {
                        c2 = others[c2];
                        codeSize[c2]++;
                    }
                }

                var bits = new int[MaxCodeLen + 1];
                for (int i = 0; i <= 256; i++)
                {
                    if (codeSize[i] == 0) continue;
                    if (codeSize[i] > MaxCodeLen) throw new InvalidOperationException("Huffman code length overflow.");
                    bits[codeSize[i]]++;
                }
                // Limit code lengths to 16 bits.
                for (int i = MaxCodeLen; i > 16; i--)
                {
                    while (bits[i] > 0)
                    {
                        int j = i - 2;
                        while (bits[j] == 0) j--;
                        bits[i] -= 2;
                        bits[i - 1]++;
                        bits[j + 1] += 2;
                        bits[j]--;
                    }
                }
                // Remove the reserved code point (it has the longest code).
                int k = 16;
                while (bits[k] == 0) k--;
                bits[k]--;

                var t = new HuffTable();
                int total = 0;
                for (int i = 1; i <= 16; i++)
                {
                    t.Bits[i] = (byte)bits[i];
                    total += bits[i];
                }
                t.Values = new byte[total];
                int p = 0;
                for (int len = 1; len <= MaxCodeLen && p < total; len++)
                    for (int s = 0; s < 256 && p < total; s++)
                        if (codeSize[s] == len) t.Values[p++] = (byte)s;

                // Canonical codes (ITU T.81 Annex C).
                int code = 0;
                p = 0;
                for (int len = 1; len <= 16; len++)
                {
                    for (int i = 0; i < t.Bits[len]; i++)
                    {
                        int sym = t.Values[p++];
                        t.Code[sym] = code;
                        t.Size[sym] = len;
                        code++;
                    }
                    code <<= 1;
                }
                return t;
            }

            public byte[] Serialize(byte classAndId)
            {
                var seg = new byte[17 + Values.Length];
                seg[0] = classAndId;
                for (int i = 1; i <= 16; i++) seg[i] = Bits[i];
                Buffer.BlockCopy(Values, 0, seg, 17, Values.Length);
                return seg;
            }
        }

        private sealed class JpegWriter
        {
            private byte[] _buf;
            private int _len;
            private ulong _bitBuf;
            private int _bitCount;

            public JpegWriter(int capacity)
            {
                _buf = new byte[Math.Max(1024, capacity)];
            }

            private void Ensure(int extra)
            {
                if (_len + extra <= _buf.Length) return;
                long n = Math.Max((long)_buf.Length * 2, (long)_len + extra + 1024);
                if (n > int.MaxValue - 64) n = int.MaxValue - 64;
                Array.Resize(ref _buf, (int)n);
            }

            private void Put(byte b)
            {
                if (_len == _buf.Length) Ensure(1);
                _buf[_len++] = b;
            }

            public void Marker(byte m)
            {
                Ensure(2);
                _buf[_len++] = 0xFF;
                _buf[_len++] = m;
            }

            public void Segment(byte marker, byte[] payload)
            {
                int len = payload.Length + 2;
                if (len > 65535) throw new InvalidDataException("JPEG segment too large.");
                Ensure(4 + payload.Length);
                _buf[_len++] = 0xFF;
                _buf[_len++] = marker;
                _buf[_len++] = (byte)(len >> 8);
                _buf[_len++] = (byte)len;
                Buffer.BlockCopy(payload, 0, _buf, _len, payload.Length);
                _len += payload.Length;
            }

            public void WriteBits(int value, int count)
            {
                if (count == 0) return;
                _bitBuf = (_bitBuf << count) | ((uint)value & ((1u << count) - 1));
                _bitCount += count;
                while (_bitCount >= 8)
                {
                    _bitCount -= 8;
                    byte b = (byte)(_bitBuf >> _bitCount);
                    Put(b);
                    if (b == 0xFF) Put(0x00); // byte stuffing
                }
            }

            public void FlushBits()
            {
                if (_bitCount > 0) WriteBits(0x7F, 8 - _bitCount); // pad with 1-bits
                _bitBuf = 0;
                _bitCount = 0;
            }

            public byte[] ToArray()
            {
                var r = new byte[_len];
                Buffer.BlockCopy(_buf, 0, r, 0, _len);
                return r;
            }
        }
    }
}
