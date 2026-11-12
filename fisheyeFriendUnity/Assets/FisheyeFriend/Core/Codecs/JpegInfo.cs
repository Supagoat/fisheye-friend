using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace FisheyeFriend.Core.Codecs
{
    /// <summary>
    /// Header facts and metadata read from a JPEG's marker segments (no pixel decoding).
    /// JPEG pixels are decoded by the host (Unity's ImageConversion in the app).
    /// </summary>
    public sealed class JpegInfo
    {
        public int Width;
        public int Height;
        public int Components;
        public int Precision = 8;
        public bool Progressive;
        public bool Arithmetic;
        /// <summary>Adobe APP14 colour transform (0 = none/CMYK, 1 = YCbCr, 2 = YCCK), or -1 if absent.</summary>
        public int AdobeTransform = -1;
        /// <summary>EXIF as a TIFF structure (no "Exif\0\0" prefix), or null.</summary>
        public byte[] Exif;
        /// <summary>ICC profile reassembled from APP2 chunks, or null.</summary>
        public byte[] IccProfile;

        public bool IsGrayscale => Components == 1;
        public bool IsCmyk => Components == 4;

        public static bool IsJpeg(byte[] data)
        {
            return data != null && data.Length > 3 && data[0] == 0xFF && data[1] == 0xD8 && data[2] == 0xFF;
        }

        public static JpegInfo Read(byte[] d)
        {
            if (!IsJpeg(d)) throw new InvalidDataException("Not a JPEG file.");
            var info = new JpegInfo();
            var icc = new SortedDictionary<int, byte[]>();
            int pos = 2;
            while (pos + 4 <= d.Length)
            {
                if (d[pos] != 0xFF)
                {
                    pos++;
                    continue;
                }
                byte m = d[pos + 1];
                if (m == 0xFF)
                {
                    pos++; // fill byte
                    continue;
                }
                pos += 2;
                if (m == 0xD8 || m == 0x01 || (m >= 0xD0 && m <= 0xD7)) continue; // no length
                if (m == 0xD9) break;
                if (pos + 2 > d.Length) break;
                int len = (d[pos] << 8) | d[pos + 1];
                if (len < 2 || pos + len > d.Length) break;
                int p = pos + 2;
                int n = len - 2;

                if (m == 0xE1 && n >= 6 && d[p] == 'E' && d[p + 1] == 'x' && d[p + 2] == 'i' && d[p + 3] == 'f' && d[p + 4] == 0 && info.Exif == null)
                {
                    info.Exif = new byte[n - 6];
                    Buffer.BlockCopy(d, p + 6, info.Exif, 0, n - 6);
                }
                else if (m == 0xE2 && n >= 14 && Encoding.ASCII.GetString(d, p, 11) == "ICC_PROFILE" && d[p + 11] == 0)
                {
                    int seq = d[p + 12];
                    var chunk = new byte[n - 14];
                    Buffer.BlockCopy(d, p + 14, chunk, 0, n - 14);
                    icc[seq] = chunk;
                }
                else if (m == 0xEE && n >= 12 && Encoding.ASCII.GetString(d, p, 5) == "Adobe")
                {
                    info.AdobeTransform = d[p + 11];
                }
                else if (m >= 0xC0 && m <= 0xCF && m != 0xC4 && m != 0xC8 && m != 0xCC && n >= 6)
                {
                    info.Precision = d[p];
                    info.Height = (d[p + 1] << 8) | d[p + 2];
                    info.Width = (d[p + 3] << 8) | d[p + 4];
                    info.Components = d[p + 5];
                    info.Progressive = m == 0xC2 || m == 0xC6 || m == 0xCA || m == 0xCE;
                    info.Arithmetic = m >= 0xC9;
                }
                else if (m == 0xDA)
                {
                    break; // start of scan: all header segments have been seen
                }
                pos += len;
            }

            if (icc.Count > 0)
            {
                int total = 0;
                foreach (var c in icc.Values) total += c.Length;
                info.IccProfile = new byte[total];
                int o = 0;
                foreach (var c in icc.Values)
                {
                    Buffer.BlockCopy(c, 0, info.IccProfile, o, c.Length);
                    o += c.Length;
                }
            }
            return info;
        }
    }
}
