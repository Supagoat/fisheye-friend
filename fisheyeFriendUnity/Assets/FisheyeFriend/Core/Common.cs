using System;
using System.Threading;

namespace FisheyeFriend.Core
{
    /// <summary>Resampling algorithms available for the de-squeeze.</summary>
    public enum ResampleAlgorithm
    {
        /// <summary>Elliptical weighted average with a Jinc-windowed Jinc ("EWA Lanczos Sharp").</summary>
        EwaLanczosSharp = 0,
        /// <summary>Separable 3-lobe Lanczos (windowed sinc).</summary>
        Lanczos3 = 1,
        /// <summary>Separable 4-lobe Lanczos.</summary>
        Lanczos4 = 2,
        /// <summary>Separable Panorama Tools 36-tap spline (3-lobe).</summary>
        Spline36 = 3,
        /// <summary>Separable Panorama Tools 64-tap spline (4-lobe).</summary>
        Spline64 = 4,
        /// <summary>John Costella's Magic Kernel Sharp 2021.</summary>
        MagicKernelSharp2021 = 5,
        /// <summary>Cubic, B=0 C=0.5.</summary>
        CatmullRom = 6,
        /// <summary>Cubic, B=C=1/3.</summary>
        Mitchell = 7,
        /// <summary>Triangle filter.</summary>
        Bilinear = 8,
        /// <summary>Nearest neighbour when enlarging, exact area average when reducing.</summary>
        Box = 9,
    }

    /// <summary>The light space a resampling pass is computed in.</summary>
    public enum LightSpace
    {
        /// <summary>Physically linear light (sRGB decoded). Keeps the brightness of fine detail correct.</summary>
        Linear = 0,
        /// <summary>The file's own gamma-encoded (sRGB) values. Gentlest halos on dark edges.</summary>
        Gamma = 1,
        /// <summary>
        /// Linear light through an inverse sigmoid (mpv's "sigmoid-upscaling", ImageMagick's sigmoidal-contrast resize):
        /// compresses overshoot near black and white so sharp filters don't produce halos.
        /// Images with alpha use Linear instead.
        /// </summary>
        Sigmoidal = 2,
    }

    /// <summary>
    /// Thread-safe progress reporting that the UI can poll. Work is reported as a fraction of the
    /// current stage, which is mapped into an overall range set with <see cref="BeginStage"/>.
    /// </summary>
    public sealed class ProgressTracker
    {
        private float _value;
        private string _stage = "";
        private float _rangeStart;
        private float _rangeSpan = 1f;
        private readonly ProgressTracker _parent;
        private readonly float _pa, _pb;

        public ProgressTracker() { }

        private ProgressTracker(ProgressTracker parent, float a, float b)
        {
            _parent = parent;
            _pa = a;
            _pb = b;
        }

        /// <summary>A tracker whose 0..1 maps onto [a, b] of this tracker's current stage (null-safe helper).</summary>
        public static ProgressTracker Slice(ProgressTracker parent, float a, float b)
        {
            return parent == null ? null : new ProgressTracker(parent, a, b);
        }

        /// <summary>Overall progress, 0..1.</summary>
        public float Value => Volatile.Read(ref _value);

        /// <summary>Name of the current stage, for display.</summary>
        public string Stage => Volatile.Read(ref _stage);

        public void BeginStage(string stage, float start, float span)
        {
            Volatile.Write(ref _stage, stage);
            _rangeStart = start;
            _rangeSpan = span;
            Volatile.Write(ref _value, start);
        }

        /// <summary>Reports progress within the current stage (0..1).</summary>
        public void Report(float fractionOfStage)
        {
            if (fractionOfStage < 0) fractionOfStage = 0;
            if (fractionOfStage > 1) fractionOfStage = 1;
            if (_parent != null)
            {
                _parent.Report(_pa + (_pb - _pa) * fractionOfStage);
                return;
            }
            Volatile.Write(ref _value, _rangeStart + _rangeSpan * fractionOfStage);
        }

        public void Reset()
        {
            Volatile.Write(ref _stage, "");
            _rangeStart = 0;
            _rangeSpan = 1;
            Volatile.Write(ref _value, 0f);
        }
    }
}
