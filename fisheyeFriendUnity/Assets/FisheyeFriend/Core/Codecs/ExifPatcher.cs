using System;

namespace FisheyeFriend.Core.Codecs
{
    /// <summary>
    /// Reads and patches an EXIF (TIFF) block: the orientation tag, the pixel-dimension tags, and the link to
    /// the embedded thumbnail. Everything else - camera, lens, exposure, GPS - is left untouched.
    /// Blocks are patched on a copy; if one cannot be parsed it is returned unchanged.
    /// </summary>
    public static class ExifPatcher
    {
        private const int TagImageWidth = 0x0100;
        private const int TagImageLength = 0x0101;
        private const int TagOrientation = 0x0112;
        private const int TagExifIfd = 0x8769;
        private const int TagPixelXDimension = 0xA002;
        private const int TagPixelYDimension = 0xA003;

        /// <summary>
        /// The EXIF orientation (1-8, see <see cref="Imaging.Orientation"/>) stored in IFD0, or 0 if the block
        /// has no valid orientation tag.
        /// </summary>
        public static int ReadOrientation(byte[] exif)
        {
            try
            {
                bool le;
                int entry = FindIfd0Entry(exif, TagOrientation, out le);
                if (entry < 0 || U16(exif, entry + 2, le) != 3 || U32(exif, entry + 4, le) < 1) return 0;
                int v = U16(exif, entry + 8, le);
                return v >= 1 && v <= 8 ? v : 0;
            }
            catch (Exception)
            {
                return 0;
            }
        }

        /// <summary>Returns a copy with the IFD0 orientation tag set to <paramref name="orientation"/> (if the tag exists).</summary>
        public static byte[] WithOrientation(byte[] exif, int orientation)
        {
            try
            {
                bool le;
                int entry = FindIfd0Entry(exif, TagOrientation, out le);
                if (entry < 0 || U16(exif, entry + 2, le) != 3) return exif;
                var t = (byte[])exif.Clone();
                W16(t, entry + 8, le, orientation);
                return t;
            }
            catch (Exception)
            {
                return exif;
            }
        }

        /// <summary>
        /// Returns a copy without the link to the embedded thumbnail (IFD1). The thumbnail is a small copy of the
        /// squeezed capture in its stored orientation, so after a de-squeeze it would show the wrong shape; with the
        /// link removed, viewers build a fresh thumbnail from the image itself.
        /// </summary>
        public static byte[] WithoutThumbnail(byte[] exif)
        {
            try
            {
                bool le;
                int ifd0 = Ifd0(exif, out le);
                if (ifd0 < 0) return exif;
                int next = ifd0 + 2 + U16(exif, ifd0, le) * 12;
                if (next + 4 > exif.Length || U32(exif, next, le) == 0) return exif;
                var t = (byte[])exif.Clone();
                W32(t, next, le, 0);
                return t;
            }
            catch (Exception)
            {
                return exif;
            }
        }

        /// <summary>Offset of IFD0, or -1 if the block is not a TIFF structure.</summary>
        private static int Ifd0(byte[] t, out bool le)
        {
            le = true;
            if (t == null || t.Length < 8) return -1;
            if (t[0] == 'I' && t[1] == 'I') le = true;
            else if (t[0] == 'M' && t[1] == 'M') le = false;
            else return -1;
            if (U16(t, 2, le) != 42) return -1;
            long ifd0 = U32(t, 4, le);
            return ifd0 >= 8 && ifd0 + 2 <= t.Length ? (int)ifd0 : -1;
        }

        /// <summary>Offset of the 12-byte IFD0 entry for <paramref name="tag"/>, or -1.</summary>
        private static int FindIfd0Entry(byte[] t, int tag, out bool le)
        {
            int ifd0 = Ifd0(t, out le);
            if (ifd0 < 0) return -1;
            int count = U16(t, ifd0, le);
            for (int i = 0; i < count; i++)
            {
                int e = ifd0 + 2 + i * 12;
                if (e + 12 > t.Length) return -1;
                if (U16(t, e, le) == tag) return e;
            }
            return -1;
        }

        public static byte[] WithDimensions(byte[] exif, int width, int height)
        {
            if (exif == null || exif.Length < 8) return exif;
            var t = (byte[])exif.Clone();
            try
            {
                bool le;
                if (t[0] == 'I' && t[1] == 'I') le = true;
                else if (t[0] == 'M' && t[1] == 'M') le = false;
                else return exif;
                if (U16(t, 2, le) != 42) return exif;
                int ifd0 = (int)U32(t, 4, le);
                int exifIfd = -1;
                PatchIfd(t, ifd0, le, width, height, TagImageWidth, TagImageLength, ref exifIfd);
                if (exifIfd > 0)
                {
                    int dummy = -1;
                    PatchIfd(t, exifIfd, le, width, height, TagPixelXDimension, TagPixelYDimension, ref dummy);
                }
                return t;
            }
            catch (Exception)
            {
                return exif;
            }
        }

        private static void PatchIfd(byte[] t, int ifd, bool le, int width, int height, int tagW, int tagH, ref int exifIfd)
        {
            if (ifd < 8 || ifd + 2 > t.Length) return;
            int count = U16(t, ifd, le);
            for (int i = 0; i < count; i++)
            {
                int e = ifd + 2 + i * 12;
                if (e + 12 > t.Length) return;
                int tag = U16(t, e, le);
                int type = U16(t, e + 2, le);
                uint n = U32(t, e + 4, le);
                if (tag == TagExifIfd && (type == 4 || type == 13) && n == 1)
                {
                    exifIfd = (int)U32(t, e + 8, le);
                }
                else if ((tag == tagW || tag == tagH) && n == 1)
                {
                    int v = tag == tagW ? width : height;
                    if (type == 3 && v <= 0xFFFF)
                    {
                        W16(t, e + 8, le, v);
                    }
                    else if (type == 4)
                    {
                        W32(t, e + 8, le, (uint)v);
                    }
                }
            }
        }

        private static int U16(byte[] b, int p, bool le)
        {
            return le ? b[p] | (b[p + 1] << 8) : (b[p] << 8) | b[p + 1];
        }

        private static uint U32(byte[] b, int p, bool le)
        {
            return le
                ? (uint)(b[p] | (b[p + 1] << 8) | (b[p + 2] << 16) | (b[p + 3] << 24))
                : (uint)((b[p] << 24) | (b[p + 1] << 16) | (b[p + 2] << 8) | b[p + 3]);
        }

        private static void W16(byte[] b, int p, bool le, int v)
        {
            if (le)
            {
                b[p] = (byte)v;
                b[p + 1] = (byte)(v >> 8);
            }
            else
            {
                b[p] = (byte)(v >> 8);
                b[p + 1] = (byte)v;
            }
        }

        private static void W32(byte[] b, int p, bool le, uint v)
        {
            if (le)
            {
                b[p] = (byte)v;
                b[p + 1] = (byte)(v >> 8);
                b[p + 2] = (byte)(v >> 16);
                b[p + 3] = (byte)(v >> 24);
            }
            else
            {
                b[p] = (byte)(v >> 24);
                b[p + 1] = (byte)(v >> 16);
                b[p + 2] = (byte)(v >> 8);
                b[p + 3] = (byte)v;
            }
        }
    }
}
