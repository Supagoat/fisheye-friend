using System;
using System.IO;
using System.IO.Compression;

namespace FisheyeFriend.Core.Codecs
{
    /// <summary>CRC-32 (IEEE 802.3), as used by PNG chunks.</summary>
    internal static class Crc32
    {
        private static readonly uint[] Table = BuildTable();

        private static uint[] BuildTable()
        {
            var t = new uint[256];
            for (uint n = 0; n < 256; n++)
            {
                uint c = n;
                for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
                t[n] = c;
            }
            return t;
        }

        public static uint Update(uint crc, byte[] data, int offset, int count)
        {
            uint c = crc;
            for (int i = 0; i < count; i++) c = Table[(c ^ data[offset + i]) & 0xFF] ^ (c >> 8);
            return c;
        }

        public static uint Compute(byte[] a, byte[] b)
        {
            uint c = 0xFFFFFFFFu;
            c = Update(c, a, 0, a.Length);
            if (b != null) c = Update(c, b, 0, b.Length);
            return c ^ 0xFFFFFFFFu;
        }
    }

    /// <summary>zlib (RFC 1950) wrapper around the platform's raw DEFLATE implementation.</summary>
    internal static class Zlib
    {
        public static uint Adler32(byte[] data, int offset, int count)
        {
            const uint Mod = 65521;
            uint a = 1, b = 0;
            int i = offset, end = offset + count;
            while (i < end)
            {
                int n = Math.Min(5552, end - i); // largest n with no 32-bit overflow before the modulo
                for (int k = 0; k < n; k++)
                {
                    a += data[i++];
                    b += a;
                }
                a %= Mod;
                b %= Mod;
            }
            return (b << 16) | a;
        }

        private const int ParallelChunk = 4 << 20;
        private static int s_parallelState; // 0 = untested, 1 = works, -1 = unavailable

        /// <summary>
        /// zlib-compresses data. Large inputs are split into 4 MB pieces compressed on all cores and
        /// joined at sync-flush boundaries (the technique pigz uses), which yields one standard zlib
        /// stream. The platform's DEFLATE is checked once for sync-flush support; if it lacks it,
        /// compression falls back to a single thread.
        /// </summary>
        public static byte[] Compress(byte[] data, int offset, int count, CompressionLevel level)
        {
            byte[][] parts;
            if (count > 2 * ParallelChunk && Environment.ProcessorCount > 1 && ParallelFlushWorks())
                parts = CompressParts(data, offset, count, level, ParallelChunk);
            else
                parts = new[] { RawDeflate(data, offset, count, level, true) };

            long total = 2 + 4;
            foreach (var p in parts) total += p.Length;
            var result = new byte[total];
            result[0] = 0x78;
            result[1] = level == CompressionLevel.Optimal ? (byte)0xDA : (byte)0x01;
            int pos = 2;
            foreach (var p in parts)
            {
                Buffer.BlockCopy(p, 0, result, pos, p.Length);
                pos += p.Length;
            }
            uint ad = Adler32(data, offset, count);
            result[pos] = (byte)(ad >> 24);
            result[pos + 1] = (byte)(ad >> 16);
            result[pos + 2] = (byte)(ad >> 8);
            result[pos + 3] = (byte)ad;
            return result;
        }

        private static byte[][] CompressParts(byte[] data, int offset, int count, CompressionLevel level, int chunk)
        {
            int n = (count + chunk - 1) / chunk;
            var parts = new byte[n][];
            System.Threading.Tasks.Parallel.For(0, n, i =>
            {
                int start = offset + i * chunk;
                int len = Math.Min(chunk, offset + count - start);
                parts[i] = RawDeflate(data, start, len, level, i == n - 1);
            });
            return parts;
        }

        /// <summary>Raw DEFLATE. When final is false the stream ends with a sync flush (byte-aligned, not final).</summary>
        private static byte[] RawDeflate(byte[] data, int offset, int count, CompressionLevel level, bool final)
        {
            using (var ms = new MemoryStream(count / 2 + 64))
            {
                var ds = new DeflateStream(ms, level, true);
                ds.Write(data, offset, count);
                if (final)
                {
                    ds.Dispose();
                    return ms.ToArray();
                }
                ds.Flush();
                byte[] bytes = ms.ToArray();
                ds.Dispose(); // writes a final block after the bytes we already took; ignored
                return bytes;
            }
        }

        private static bool ParallelFlushWorks()
        {
            int state = System.Threading.Volatile.Read(ref s_parallelState);
            if (state != 0) return state > 0;
            bool ok;
            try
            {
                var sample = new byte[3 * 70000];
                var rnd = new Random(12345);
                for (int i = 0; i < sample.Length; i++) sample[i] = (byte)(i % 251 < 120 ? rnd.Next(256) : (i * 7) & 0xFF);
                var parts = CompressParts(sample, 0, sample.Length, CompressionLevel.Fastest, 70000);
                using (var ms = new MemoryStream())
                {
                    foreach (var p in parts) ms.Write(p, 0, p.Length);
                    ms.Position = 0;
                    var back = new byte[sample.Length];
                    using (var ds = new DeflateStream(ms, CompressionMode.Decompress))
                    {
                        int read = 0, nr;
                        while (read < back.Length && (nr = ds.Read(back, read, back.Length - read)) > 0) read += nr;
                        ok = read == back.Length;
                    }
                    for (int i = 0; ok && i < back.Length; i++) ok = back[i] == sample[i];
                }
            }
            catch (Exception)
            {
                ok = false;
            }
            System.Threading.Volatile.Write(ref s_parallelState, ok ? 1 : -1);
            return ok;
        }

        /// <summary>Inflates a zlib stream. If expectedLength is known (>= 0) exactly that many bytes are read.</summary>
        public static byte[] Decompress(byte[] data, int offset, int count, long expectedLength)
        {
            if (count < 2) throw new InvalidDataException("zlib stream is too short.");
            int cmf = data[offset], flg = data[offset + 1];
            if ((cmf & 0x0F) != 8 || ((cmf << 8) | flg) % 31 != 0)
                throw new InvalidDataException("Not a valid zlib stream.");
            if ((flg & 0x20) != 0) throw new InvalidDataException("zlib preset dictionaries are not supported.");

            using (var input = new MemoryStream(data, offset + 2, count - 2, false))
            using (var ds = new DeflateStream(input, CompressionMode.Decompress))
            {
                if (expectedLength >= 0)
                {
                    if (expectedLength > int.MaxValue) throw new InvalidDataException("Image data is too large.");
                    var result = new byte[expectedLength];
                    int read = 0;
                    while (read < result.Length)
                    {
                        int n = ds.Read(result, read, result.Length - read);
                        if (n <= 0) throw new InvalidDataException("Image data is truncated.");
                        read += n;
                    }
                    return result;
                }
                using (var outMs = new MemoryStream())
                {
                    ds.CopyTo(outMs);
                    return outMs.ToArray();
                }
            }
        }
    }
}
