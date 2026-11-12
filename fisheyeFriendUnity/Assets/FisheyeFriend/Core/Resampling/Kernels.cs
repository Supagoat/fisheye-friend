using System;

namespace FisheyeFriend.Core.Resampling
{
    /// <summary>A 1-D reconstruction filter: weight as a function of distance (in source pixels at 1:1).</summary>
    public sealed class Kernel1D
    {
        public readonly string Name;
        /// <summary>Support radius at 1:1 scale; weights are zero beyond it.</summary>
        public readonly double Radius;
        private readonly Func<double, double> _f;

        public Kernel1D(string name, double radius, Func<double, double> f)
        {
            Name = name;
            Radius = radius;
            _f = f;
        }

        public double Evaluate(double x)
        {
            x = Math.Abs(x);
            return x >= Radius ? 0.0 : _f(x);
        }
    }

    /// <summary>The separable kernels, plus the radial Jinc kernel used by the EWA resampler.</summary>
    public static class Kernels
    {
        public static Kernel1D For(ResampleAlgorithm algorithm)
        {
            switch (algorithm)
            {
                case ResampleAlgorithm.Lanczos3: return Lanczos(3);
                case ResampleAlgorithm.Lanczos4: return Lanczos(4);
                case ResampleAlgorithm.Spline36: return Spline36;
                case ResampleAlgorithm.Spline64: return Spline64;
                case ResampleAlgorithm.MagicKernelSharp2021: return MagicKernelSharp2021;
                case ResampleAlgorithm.CatmullRom: return CatmullRom;
                case ResampleAlgorithm.Mitchell: return Mitchell;
                case ResampleAlgorithm.Bilinear: return Triangle;
                case ResampleAlgorithm.Box: return Box;
                default: throw new ArgumentException(algorithm + " is not a separable kernel.");
            }
        }

        public static double Sinc(double x)
        {
            if (Math.Abs(x) < 1e-12) return 1.0;
            double px = Math.PI * x;
            return Math.Sin(px) / px;
        }

        public static Kernel1D Lanczos(int lobes)
        {
            return new Kernel1D("Lanczos" + lobes, lobes, x => Sinc(x) * Sinc(x / lobes));
        }

        public static readonly Kernel1D Box = new Kernel1D("Box", 0.5, x => 1.0);

        public static readonly Kernel1D Triangle = new Kernel1D("Triangle", 1.0, x => 1.0 - x);

        public static readonly Kernel1D CatmullRom = new Kernel1D("CatmullRom", 2.0, x => Cubic(x, 0.0, 0.5));

        public static readonly Kernel1D Mitchell = new Kernel1D("Mitchell", 2.0, x => Cubic(x, 1.0 / 3.0, 1.0 / 3.0));

        /// <summary>Mitchell-Netravali family of cubics with parameters B and C.</summary>
        public static double Cubic(double x, double b, double c)
        {
            x = Math.Abs(x);
            if (x < 1.0)
                return ((12 - 9 * b - 6 * c) * x * x * x + (-18 + 12 * b + 6 * c) * x * x + (6 - 2 * b)) / 6.0;
            if (x < 2.0)
                return ((-b - 6 * c) * x * x * x + (6 * b + 30 * c) * x * x + (-12 * b - 48 * c) * x + (8 * b + 24 * c)) / 6.0;
            return 0.0;
        }

        /// <summary>Helmut Dersch's Panorama Tools spline, 36 taps for 2-D (3 lobes).</summary>
        public static readonly Kernel1D Spline36 = new Kernel1D("Spline36", 3.0, x =>
        {
            if (x < 1.0) return ((13.0 / 11.0 * x - 453.0 / 209.0) * x - 3.0 / 209.0) * x + 1.0;
            if (x < 2.0) { x -= 1.0; return ((-6.0 / 11.0 * x + 270.0 / 209.0) * x - 156.0 / 209.0) * x; }
            x -= 2.0; return ((1.0 / 11.0 * x - 45.0 / 209.0) * x + 26.0 / 209.0) * x;
        });

        /// <summary>Helmut Dersch's Panorama Tools spline, 64 taps for 2-D (4 lobes).</summary>
        public static readonly Kernel1D Spline64 = new Kernel1D("Spline64", 4.0, x =>
        {
            if (x < 1.0) return ((49.0 / 41.0 * x - 6387.0 / 2911.0) * x - 3.0 / 2911.0) * x + 1.0;
            if (x < 2.0) { x -= 1.0; return ((-24.0 / 41.0 * x + 4032.0 / 2911.0) * x - 2328.0 / 2911.0) * x; }
            if (x < 3.0) { x -= 2.0; return ((6.0 / 41.0 * x - 1008.0 / 2911.0) * x + 582.0 / 2911.0) * x; }
            x -= 3.0; return ((-1.0 / 41.0 * x + 168.0 / 2911.0) * x - 97.0 / 2911.0) * x;
        });

        /// <summary>
        /// Magic Kernel Sharp 2021 (John Costella), as implemented in libvips (VIPS_KERNEL_MKS2021).
        /// </summary>
        public static readonly Kernel1D MagicKernelSharp2021 = new Kernel1D("MKS2021", 4.5, x =>
        {
            if (x >= 3.5) return (4.0 * x * x - 36.0 * x + 81.0) / -1152.0;
            if (x >= 2.5) return (4.0 * x * x - 27.0 * x + 45.0) / 144.0;
            if (x >= 1.5) return (24.0 * x * x - 113.0 * x + 130.0) / -144.0;
            if (x >= 0.5) return (140.0 * x * x - 379.0 * x + 239.0) / 144.0;
            return 577.0 / 576.0 - 239.0 * x * x / 144.0;
        });

        // ---------------------------------------------------------------------------------------
        // Radial (EWA) kernel
        // ---------------------------------------------------------------------------------------

        /// <summary>First zero of Jinc(x) = 2*J1(pi x)/(pi x).</summary>
        public const double JincZero1 = 1.2196698912665045;
        /// <summary>Third zero of Jinc - the support of a 3-lobe Jinc-windowed Jinc ("EWA Lanczos").</summary>
        public const double JincZero3 = 3.2383154841662362;
        /// <summary>
        /// Nicolas Robidoux's "LanczosSharp" blur: scales EWA Lanczos so that it reproduces
        /// constant and linear signals as closely as possible on a regular grid (same value mpv and ImageMagick use).
        /// </summary>
        public const double LanczosSharpBlur = 0.9812505644269356;

        /// <summary>Jinc(x) = 2*J1(pi x)/(pi x), normalised so Jinc(0) = 1.</summary>
        public static double Jinc(double x)
        {
            x = Math.Abs(x);
            if (x < 1e-9) return 1.0;
            double z = Math.PI * x;
            return 2.0 * BesselJ1(z) / z;
        }

        /// <summary>Jinc-windowed Jinc with 3 lobes (EWA Lanczos).</summary>
        public static double JincLanczos3(double r)
        {
            r = Math.Abs(r);
            if (r >= JincZero3) return 0.0;
            return Jinc(r) * Jinc(r * JincZero1 / JincZero3);
        }

        /// <summary>
        /// Bessel function of the first kind, order 1, by its power series. Used only to build lookup
        /// tables; accurate to ~1e-12 over the range we need (|z| &lt; 12).
        /// </summary>
        public static double BesselJ1(double z)
        {
            double half = z * 0.5;
            double term = half; // m = 0: (z/2)^1 / (0! 1!)
            double sum = term;
            double h2 = half * half;
            for (int m = 1; m < 60; m++)
            {
                term *= -h2 / (m * (m + 1.0));
                sum += term;
                if (Math.Abs(term) < 1e-17 * Math.Abs(sum)) break;
            }
            return sum;
        }
    }
}
