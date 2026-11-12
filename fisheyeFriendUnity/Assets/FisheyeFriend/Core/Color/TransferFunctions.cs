using System;

namespace FisheyeFriend.Core.Color
{
    /// <summary>
    /// sRGB transfer function with lookup tables accurate enough for 16-bit output.
    /// Decoding integer samples is exact (one table per bit depth); encoding uses a 64K-entry
    /// table with linear interpolation (max error about 1e-7, far below one 16-bit code).
    /// </summary>
    public static class Srgb
    {
        private const int EncodeLutSize = 65536;
        private static readonly float[] s_decode8 = BuildDecode(255);
        private static float[] s_decode16;
        private static readonly float[] s_encode = BuildEncode();
        private static readonly object s_lock = new object();

        public static double ToLinear(double v)
        {
            return v <= 0.04045 ? v / 12.92 : Math.Pow((v + 0.055) / 1.055, 2.4);
        }

        public static double FromLinear(double l)
        {
            return l <= 0.0031308 ? l * 12.92 : 1.055 * Math.Pow(l, 1.0 / 2.4) - 0.055;
        }

        /// <summary>Table mapping 8-bit sRGB codes to linear light in [0,1].</summary>
        public static float[] Decode8 => s_decode8;

        /// <summary>Table mapping 16-bit sRGB codes to linear light in [0,1] (built on first use).</summary>
        public static float[] Decode16
        {
            get
            {
                if (s_decode16 == null)
                    lock (s_lock)
                        if (s_decode16 == null) s_decode16 = BuildDecode(65535);
                return s_decode16;
            }
        }

        /// <summary>Linear light -> sRGB-encoded value in [0,1]. Input is clamped to [0,1].</summary>
        public static float Encode(float linear)
        {
            if (!(linear > 0f)) return 0f; // also catches NaN
            if (linear >= 1f) return 1f;
            float p = linear * (EncodeLutSize - 1);
            int i = (int)p;
            float f = p - i;
            float a = s_encode[i];
            return a + (s_encode[i + 1] - a) * f;
        }

        /// <summary>sRGB-encoded value in [0,1] -> linear light (input clamped), via interpolated table.</summary>
        public static float DecodeFloat(float encoded)
        {
            if (!(encoded > 0f)) return 0f;
            if (encoded >= 1f) return 1f;
            var t = Decode16;
            float p = encoded * 65535f;
            int i = (int)p;
            if (i >= 65535) return 1f;
            float f = p - i;
            float a = t[i];
            return a + (t[i + 1] - a) * f;
        }

        private static float[] BuildDecode(int max)
        {
            var t = new float[max + 1];
            for (int i = 0; i <= max; i++) t[i] = (float)ToLinear(i / (double)max);
            return t;
        }

        private static float[] BuildEncode()
        {
            var t = new float[EncodeLutSize + 1];
            for (int i = 0; i <= EncodeLutSize; i++)
                t[i] = (float)FromLinear(Math.Min(1.0, i / (double)(EncodeLutSize - 1)));
            return t;
        }
    }

    /// <summary>
    /// Sigmoidal light space used while enlarging (as in mpv's sigmoid-upscaling and ImageMagick's
    /// sigmoidal-contrast resizing). Linear light is pushed through an inverse logistic curve before
    /// resampling and back afterwards, so overshoot from sharp filters is compressed near black and
    /// white instead of appearing as halos. Defaults follow mpv: centre 0.75, slope 6.5.
    /// </summary>
    public sealed class SigmoidCurve
    {
        private const int LutSize = 65536;
        public readonly double Center;
        public readonly double Slope;
        private readonly float[] _toSig;
        private readonly float[] _fromSig;

        private static SigmoidCurve s_default;
        public static SigmoidCurve Default => s_default ?? (s_default = new SigmoidCurve(0.75, 6.5));

        public SigmoidCurve(double center, double slope)
        {
            Center = center;
            Slope = slope;
            double offset = Logistic(0.0);
            double scale = Logistic(1.0) - offset;
            _toSig = new float[LutSize + 1];
            _fromSig = new float[LutSize + 1];
            for (int i = 0; i <= LutSize; i++)
            {
                double x = Math.Min(1.0, i / (double)(LutSize - 1));
                // linear -> sigmoid space: inverse logistic of the rescaled value
                double l = x * scale + offset;
                double s = Center - Math.Log(1.0 / l - 1.0) / Slope;
                _toSig[i] = (float)Math.Max(0.0, Math.Min(1.0, s));
                // sigmoid space -> linear
                double back = (Logistic(x) - offset) / scale;
                _fromSig[i] = (float)Math.Max(0.0, Math.Min(1.0, back));
            }
        }

        private double Logistic(double x)
        {
            return 1.0 / (1.0 + Math.Exp(Slope * (Center - x)));
        }

        public float ToSigmoid(float linear)
        {
            return Lookup(_toSig, linear);
        }

        public float FromSigmoid(float sig)
        {
            return Lookup(_fromSig, sig);
        }

        private static float Lookup(float[] t, float v)
        {
            if (!(v > 0f)) return t[0];
            if (v >= 1f) return t[LutSize - 1];
            float p = v * (LutSize - 1);
            int i = (int)p;
            float f = p - i;
            float a = t[i];
            return a + (t[i + 1] - a) * f;
        }
    }
}
