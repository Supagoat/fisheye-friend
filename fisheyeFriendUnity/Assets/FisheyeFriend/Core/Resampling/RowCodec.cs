using System;
using FisheyeFriend.Core.Color;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Resampling
{
    /// <summary>
    /// Converts image rows between integer samples and float working buffers.
    /// Colour channels are decoded to linear light (or kept as gamma values), alpha is always linear,
    /// and colour is premultiplied by alpha so transparent pixels don't bleed their colour into neighbours.
    /// </summary>
    internal sealed class RowCodec
    {
        private static float[] s_gamma8;
        private static float[] s_gamma16;

        private readonly PixelImage _img;
        private readonly float[] _colorLut;
        private readonly float _alphaScale;

        public RowCodec(PixelImage img, LightSpace decodeTo)
        {
            _img = img;
            bool linear = decodeTo != LightSpace.Gamma;
            if (img.BitDepth == 8) _colorLut = linear ? Srgb.Decode8 : Gamma8;
            else _colorLut = linear ? Srgb.Decode16 : Gamma16;
            _alphaScale = 1f / img.MaxValue;
        }

        private static float[] Gamma8 => s_gamma8 ?? (s_gamma8 = BuildIdentity(255));
        private static float[] Gamma16 => s_gamma16 ?? (s_gamma16 = BuildIdentity(65535));

        private static float[] BuildIdentity(int max)
        {
            var t = new float[max + 1];
            for (int i = 0; i <= max; i++) t[i] = i / (float)max;
            return t;
        }

        /// <summary>Decodes source row y into dst[offset .. offset + Width*Channels).</summary>
        public void Decode(int y, float[] dst, int offset)
        {
            int ch = _img.Channels;
            int n = _img.Width;
            int s = y * _img.Stride;
            var lut = _colorLut;
            bool alpha = _img.HasAlpha;
            int cc = alpha ? ch - 1 : ch;
            if (_img.BitDepth == 8)
            {
                var a = _img.Data8;
                if (!alpha)
                {
                    int count = n * ch;
                    for (int i = 0; i < count; i++) dst[offset + i] = lut[a[s + i]];
                    return;
                }
                for (int x = 0; x < n; x++, s += ch, offset += ch)
                {
                    float al = a[s + cc] * _alphaScale;
                    for (int c = 0; c < cc; c++) dst[offset + c] = lut[a[s + c]] * al;
                    dst[offset + cc] = al;
                }
            }
            else
            {
                var a = _img.Data16;
                if (!alpha)
                {
                    int count = n * ch;
                    for (int i = 0; i < count; i++) dst[offset + i] = lut[a[s + i]];
                    return;
                }
                for (int x = 0; x < n; x++, s += ch, offset += ch)
                {
                    float al = a[s + cc] * _alphaScale;
                    for (int c = 0; c < cc; c++) dst[offset + c] = lut[a[s + c]] * al;
                    dst[offset + cc] = al;
                }
            }
        }

        /// <summary>
        /// Encodes a working row (in Linear or Gamma domain, premultiplied if the image has alpha)
        /// into row y of dst, with clamping and round-to-nearest quantisation.
        /// </summary>
        public static void Encode(float[] src, int offset, PixelImage dst, int y, LightSpace from)
        {
            int ch = dst.Channels;
            int n = dst.Width;
            int d = y * dst.Stride;
            bool alpha = dst.HasAlpha;
            int cc = alpha ? ch - 1 : ch;
            float max = dst.MaxValue;
            if (from == LightSpace.Sigmoidal) throw new ArgumentException("Convert sigmoidal rows to linear before encoding.");
            bool linear = from == LightSpace.Linear;
            if (dst.BitDepth == 8)
            {
                var o = dst.Data8;
                for (int x = 0; x < n; x++, d += ch, offset += ch)
                {
                    float inv = 1f;
                    if (alpha)
                    {
                        float al = Clamp01(src[offset + cc]);
                        o[d + cc] = (byte)(al * max + 0.5f);
                        inv = al > 1e-7f ? 1f / al : 0f;
                    }
                    for (int c = 0; c < cc; c++)
                    {
                        float v = src[offset + c] * inv;
                        v = linear ? Srgb.Encode(v) : Clamp01(v);
                        o[d + c] = (byte)(v * max + 0.5f);
                    }
                }
            }
            else
            {
                var o = dst.Data16;
                for (int x = 0; x < n; x++, d += ch, offset += ch)
                {
                    float inv = 1f;
                    if (alpha)
                    {
                        float al = Clamp01(src[offset + cc]);
                        o[d + cc] = (ushort)(al * max + 0.5f);
                        inv = al > 1e-7f ? 1f / al : 0f;
                    }
                    for (int c = 0; c < cc; c++)
                    {
                        float v = src[offset + c] * inv;
                        v = linear ? Srgb.Encode(v) : Clamp01(v);
                        o[d + c] = (ushort)(v * max + 0.5f);
                    }
                }
            }
        }

        /// <summary>
        /// Converts a premultiplied row with an alpha channel: colour is un-premultiplied, converted and
        /// premultiplied again (the light-space curves don't commute with alpha).
        /// </summary>
        public static void ConvertWithAlpha(float[] buf, int offset, int count, int channels, LightSpace from, LightSpace to)
        {
            if (from == to) return;
            int cc = channels - 1;
            var tmp = new float[cc];
            for (int p = offset; p < offset + count; p += channels)
            {
                float a = buf[p + cc];
                if (!(a > 1e-7f))
                {
                    for (int c = 0; c < cc; c++) buf[p + c] = 0f;
                    continue;
                }
                float inv = 1f / a;
                for (int c = 0; c < cc; c++) tmp[c] = buf[p + c] * inv;
                Convert(tmp, 0, cc, from, to);
                for (int c = 0; c < cc; c++) buf[p + c] = tmp[c] * a;
            }
        }

        /// <summary>Converts count values in place from one light space to another (via linear light).</summary>
        public static void Convert(float[] buf, int offset, int count, LightSpace from, LightSpace to)
        {
            if (from == to) return;
            var sig = SigmoidCurve.Default;
            int end = offset + count;
            if (from == LightSpace.Gamma)
            {
                if (to == LightSpace.Linear)
                    for (int i = offset; i < end; i++) buf[i] = Srgb.DecodeFloat(buf[i]);
                else
                    for (int i = offset; i < end; i++) buf[i] = sig.ToSigmoid(Srgb.DecodeFloat(buf[i]));
            }
            else if (from == LightSpace.Linear)
            {
                if (to == LightSpace.Gamma)
                    for (int i = offset; i < end; i++) buf[i] = Srgb.Encode(buf[i]);
                else
                    for (int i = offset; i < end; i++) buf[i] = sig.ToSigmoid(buf[i]);
            }
            else
            {
                if (to == LightSpace.Linear)
                    for (int i = offset; i < end; i++) buf[i] = sig.FromSigmoid(buf[i]);
                else
                    for (int i = offset; i < end; i++) buf[i] = Srgb.Encode(sig.FromSigmoid(buf[i]));
            }
        }

        private static float Clamp01(float v)
        {
            if (!(v > 0f)) return 0f;
            return v > 1f ? 1f : v;
        }
    }
}
