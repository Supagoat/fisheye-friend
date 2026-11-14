using System;
using System.Runtime.InteropServices;
using FisheyeFriend.Core.Color;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;

namespace FisheyeFriend.Core.Gpu
{
    /// <summary>
    /// Everything the GPU defish (FisheyeFriendDefish.compute) needs, prepared on the CPU without Unity: the model's
    /// constants, the grid's nodes in source pixels, the per-column and per-row cell tables, the resampler's own tables,
    /// and the photo's samples packed into 32-bit words. The Unity side (FisheyeFriend.Gpu.GpuDefisher) only uploads these
    /// and dispatches; the bench's emulator runs the same kernels on them in C#.
    /// </summary>
    public static class GpuPlan
    {
        public const int MaxThetaTerms = 8, MaxHorizontalTerms = 4, MaxVerticalRows = 4, MaxVerticalColumns = 4;

        /// <summary>Whether the GPU kernels can hold this model (the built-in one fits with room to spare).</summary>
        public static bool ModelFits(DefishModel m, out string why)
        {
            why = null;
            if (m == null) why = "there is no model (the projection is point pairs)";
            else if (m.Theta.Length > MaxThetaTerms) why = "the lens has more than " + MaxThetaTerms + " terms";
            else if (m.Horizontal.Length > MaxHorizontalTerms) why = "the projection has more than " + MaxHorizontalTerms + " horizontal terms";
            else if (m.Vertical.GetLength(0) > MaxVerticalRows || m.Vertical.GetLength(1) > MaxVerticalColumns)
                why = "the projection has more than " + MaxVerticalRows + " x " + MaxVerticalColumns + " vertical terms";
            return why == null;
        }

        /// <summary>Whether the GPU resampler does what these options ask (the separable kernels; not the radial one).</summary>
        public static bool OptionsFit(WarpOptions o, out string why)
        {
            why = o.Kernel == WarpKernel.RadialLanczos ? "the radial Lanczos filter runs on the CPU only" : null;
            return why == null;
        }

        // ================================================================== the model

        /// <summary>The model's constants as the MapModel kernel takes them (unused terms 0, which adds nothing).</summary>
        public sealed class ModelParams
        {
            public float[] Theta = new float[MaxThetaTerms];
            public float[] Horizontal = new float[MaxHorizontalTerms];
            /// <summary>Row by row, 4 x 4.</summary>
            public float[] Vertical = new float[MaxVerticalRows * MaxVerticalColumns];
            public float Panini, PaniniVertical, K, Strength, AmountX, AmountY, Aspect, RhoNorm;
            public float[] ShapeX = new float[4], ShapeY = new float[4];
            public bool Flat;

            public static ModelParams For(DefishModel m, DefishAdjustment adjustment, double aspect)
            {
                string why;
                if (!ModelFits(m, out why)) throw new NotSupportedException("The GPU can't evaluate this model: " + why + ".");
                var a = (adjustment ?? new DefishAdjustment()).Clone();
                a.Validate();
                var p = new ModelParams();
                for (int i = 0; i < m.Theta.Length; i++) p.Theta[i] = (float)m.Theta[i];
                for (int i = 0; i < m.Horizontal.Length; i++) p.Horizontal[i] = (float)m.Horizontal[i];
                for (int i = 0; i < m.Vertical.GetLength(0); i++)
                    for (int j = 0; j < m.Vertical.GetLength(1); j++) p.Vertical[i * MaxVerticalColumns + j] = (float)m.Vertical[i, j];
                p.Panini = (float)m.Panini;
                p.PaniniVertical = (float)m.PaniniVertical;
                // as DefishModelMap: the scale that keeps the photo's edge midpoints at the frame's edges
                p.K = (float)(m.Scale * m.ProjectX(m.LensAngle(0.5)) / m.ProjectX(m.LensAngle(0.5 * a.Strength)));
                p.Strength = (float)a.Strength;
                p.AmountX = (float)a.AmountX;
                p.AmountY = (float)a.AmountY;
                p.Aspect = (float)aspect;
                p.RhoNorm = (float)(1 / Math.Sqrt(0.25 + 0.25 * aspect * aspect));
                for (int i = 0; i < 4; i++)
                {
                    p.ShapeX[i] = (float)a.ShapeX[i];
                    p.ShapeY[i] = (float)a.ShapeY[i];
                }
                p.Flat = a.ShapesFlat;
                return p;
            }
        }

        /// <summary>The (u, v) of every node of the geometry's grid for this frame shape, row by row, interleaved.</summary>
        public static float[] GridPoints(double aspect, int nodesAcross)
        {
            int nu, nv;
            double u0, v0, hu, hv;
            DefishGeometry.GridLayout(aspect, nodesAcross, out nu, out nv, out u0, out v0, out hu, out hv);
            var p = new float[2 * nu * nv];
            for (int j = 0, k = 0; j < nv; j++)
            {
                double v = v0 + j * hv;
                for (int i = 0; i < nu; i++, k += 2)
                {
                    p[k] = (float)(u0 + i * hu);
                    p[k + 1] = (float)v;
                }
            }
            return p;
        }

        /// <summary>A geometry from the mapped grid points (s, t interleaved, as the MapModel kernel returns them).</summary>
        public static DefishGeometry GeometryFromMapped(DefishModel m, DefishAdjustment adjustment, double aspect, int nodesAcross, float[] mapped)
        {
            int n = mapped.Length / 2;
            var s = new float[n];
            var t = new float[n];
            for (int i = 0; i < n; i++)
            {
                s[i] = mapped[2 * i];
                t[i] = mapped[2 * i + 1];
            }
            return DefishGeometry.FromNodes(m, adjustment, aspect, nodesAcross, s, t);
        }

        // ================================================================== the resampler

        [StructLayout(LayoutKind.Sequential)]
        public struct Node
        {
            public int IX, IY;
            public float FX, FY;
        }

        [StructLayout(LayoutKind.Sequential)]
        public struct Cell
        {
            public int I;
            public float A;
        }

        /// <summary>The geometry as the Warp kernel reads it, for one photo size and output.</summary>
        public sealed class WarpGrid
        {
            public Node[] Nodes;
            public int NodesU;
            public Cell[] Cols, Rows;
            public float ColStep, RowStep;

            /// <summary>
            /// Each node in source pixels (pixel centres at whole numbers), split into a whole pixel and a fraction in double
            /// precision; and for each output column and row, the cell its centre falls in (as DefishWarpMap and
            /// DefishGeometry.Map find it).
            /// </summary>
            public static WarpGrid For(DefishGeometry g, int srcW, int srcH, int outW, int outH, DefishCanvas c)
            {
                var wg = new WarpGrid { NodesU = g.NodesU, Nodes = new Node[g.GridNodes], Cols = new Cell[outW], Rows = new Cell[outH] };
                for (int k = 0; k < wg.Nodes.Length; k++)
                {
                    float s, t;
                    g.Node(k, out s, out t);
                    double sx = s * (double)srcW - 0.5, sy = t * (double)srcH - 0.5;
                    double ix = Math.Floor(sx), iy = Math.Floor(sy);
                    wg.Nodes[k] = new Node { IX = (int)ix, IY = (int)iy, FX = (float)(sx - ix), FY = (float)(sy - iy) };
                }
                double du = (c.U1 - c.U0) / outW, dv = (c.V1 - c.V0) / outH;
                for (int x = 0; x < outW; x++)
                {
                    int i;
                    double a;
                    g.CellAtU(c.U0 + (x + 0.5) * du, out i, out a);
                    wg.Cols[x] = new Cell { I = i, A = (float)a };
                }
                for (int y = 0; y < outH; y++)
                {
                    int j;
                    double b;
                    g.CellAtV(c.V0 + (y + 0.5) * dv, out j, out b);
                    wg.Rows[y] = new Cell { I = j, A = (float)b };
                }
                wg.ColStep = (float)(du / g.CellU);
                wg.RowStep = (float)(dv / g.CellV);
                return wg;
            }
        }

        /// <summary>The resampler's options and tables as the Warp kernel takes them.</summary>
        public sealed class WarpParams
        {
            public int Channels, Depth, OutDepth, Light, WeightLutSize, ApCount;
            public bool Alpha, Antialias, SharpenOnly, Extend;
            public float AlphaScale, R, LutScale, AaStrength, Aperture, ApLogMin, ApLogMax, OutMax;
            public float[] Bg = new float[4];
            public float[] ColorLut, WeightLut, ApTable, FromSig, Encode;

            public static WarpParams For(PixelImage src, WarpOptions o)
            {
                string why;
                if (!OptionsFit(o, out why)) throw new NotSupportedException("The GPU can't resample with these options: " + why + ".");
                var p = new WarpParams
                {
                    Channels = src.Channels,
                    Depth = src.BitDepth,
                    OutDepth = o.OutputBitDepth == 0 ? src.BitDepth : o.OutputBitDepth,
                    Light = (int)o.Light,
                    Alpha = src.HasAlpha,
                    AlphaScale = 1f / src.MaxValue,
                    R = (float)(o.Kernel == WarpKernel.Bicubic ? 2 : o.Radius),
                    Antialias = o.Antialias,
                    AaStrength = (float)o.AntialiasStrength,
                    Aperture = (float)o.Aperture,
                    SharpenOnly = o.ApertureSharpenOnly,
                    Extend = o.ExtendEdges,
                    ColorLut = WarpResampler.DecodeTable(src.BitDepth, o.Light),
                    WeightLutSize = WarpResampler.WeightTableSize,
                    ApLogMin = (float)WarpResampler.ApertureLogMin,
                    ApLogMax = (float)WarpResampler.ApertureLogMax,
                    FromSig = SigmoidCurve.Default.FromTable,
                    Encode = Srgb.EncodeTable,
                };
                double R = o.Kernel == WarpKernel.Bicubic ? 2 : o.Radius;
                p.WeightLut = WarpResampler.WeightTable(o.Kernel, R);
                p.LutScale = (float)((WarpResampler.WeightTableSize - 1) / R);
                var ap = WarpResampler.ApertureTableValues();
                p.ApTable = Array.ConvertAll(ap, v => (float)v);
                p.ApCount = ap.Length;
                p.OutMax = p.OutDepth == 8 ? 255f : 65535f;
                int cc = p.Alpha ? p.Channels - 1 : p.Channels;
                for (int c = 0; c < cc; c++)
                {
                    float g = o.Background != null && c < o.Background.Length ? o.Background[c] : 0f;
                    p.Bg[c] = WarpResampler.WorkingValue(g, o.Light);
                }
                return p;
            }

            /// <summary>32-bit words per output pixel.</summary>
            public int WordsPerPixel => OutDepth == 8 ? 1 : 2;
        }

        // ================================================================== samples

        /// <summary>The photo's samples packed into 32-bit words (4 per word at 8 bits, 2 at 16; little-endian, as stored).</summary>
        public static uint[] Pack(PixelImage img)
        {
            if (img.BitDepth == 8)
            {
                var w = new uint[(img.Data8.Length + 3) / 4];
                Buffer.BlockCopy(img.Data8, 0, w, 0, img.Data8.Length);
                return w;
            }
            var w16 = new uint[(img.Data16.Length + 1) / 2];
            Buffer.BlockCopy(img.Data16, 0, w16, 0, img.Data16.Length * 2);
            return w16;
        }

        /// <summary>Puts <paramref name="rows"/> output rows from <paramref name="y0"/> (as the Warp kernel wrote them) into <paramref name="dst"/>.</summary>
        public static void Unpack(uint[] words, PixelImage dst, int y0, int rows)
        {
            int w = dst.Width, ch = dst.Channels;
            if (dst.BitDepth == 8)
            {
                var d = dst.Data8;
                for (int y = 0; y < rows; y++)
                {
                    int o = (y0 + y) * w * ch, k = y * w;
                    for (int x = 0; x < w; x++, k++, o += ch)
                    {
                        uint v = words[k];
                        d[o] = (byte)v;
                        if (ch > 1) d[o + 1] = (byte)(v >> 8);
                        if (ch > 2) d[o + 2] = (byte)(v >> 16);
                        if (ch > 3) d[o + 3] = (byte)(v >> 24);
                    }
                }
            }
            else
            {
                var d = dst.Data16;
                for (int y = 0; y < rows; y++)
                {
                    int o = (y0 + y) * w * ch, k = 2 * y * w;
                    for (int x = 0; x < w; x++, k += 2, o += ch)
                    {
                        uint a = words[k], b = words[k + 1];
                        d[o] = (ushort)a;
                        if (ch > 1) d[o + 1] = (ushort)(a >> 16);
                        if (ch > 2) d[o + 2] = (ushort)b;
                        if (ch > 3) d[o + 3] = (ushort)(b >> 16);
                    }
                }
            }
        }

        // ================================================================== a whole defish

        /// <summary>
        /// A defish prepared for the GPU, as <see cref="Defisher.Process"/> does it on the CPU: the photo turned into the
        /// lens' frame, the output's size and the part of the frame it shows, the grid, the options and the packed samples.
        /// </summary>
        public sealed class Job
        {
            public PixelImage Upright, Sensor;
            public int Turn, OutW, OutH;
            public DefishCanvas Canvas;
            public WarpOptions Options;
            public WarpGrid Grid;
            public WarpParams Params;
            public uint[] Samples;

            public static Job Prepare(PixelImage upright, DefishGeometry g, DefishSettings settings)
            {
                if (upright == null) throw new ArgumentNullException(nameof(upright));
                if (g == null) throw new ArgumentNullException(nameof(g));
                if (g.Folds)
                    throw new InvalidOperationException("These defish settings fold the picture over itself near the edges. Lower the horizontal or vertical amount, or the curve.");
                var s = settings.Clone();
                s.Validate();
                var j = new Job { Upright = upright, Turn = Defisher.SensorTurn(upright) };
                j.Sensor = Orientation.Apply(upright, j.Turn);
                Defisher.Plan(g, j.Sensor.Width, j.Sensor.Height, s, out j.Canvas, out j.OutW, out j.OutH);
                j.Options = Defisher.WarpOptionsFor(j.Sensor, s);
                j.Grid = WarpGrid.For(g, j.Sensor.Width, j.Sensor.Height, j.OutW, j.OutH, j.Canvas);
                j.Params = WarpParams.For(j.Sensor, j.Options);
                j.Samples = Pack(j.Sensor);
                return j;
            }

            /// <summary>The empty result, in the lens' frame, for <see cref="Unpack"/> to fill.</summary>
            public PixelImage NewOutput()
            {
                return new PixelImage(OutW, OutH, Sensor.Channels, Params.OutDepth)
                {
                    Metadata = Sensor.Metadata != null ? Sensor.Metadata.Clone() : new ImageMetadata(),
                    SourceFormat = Sensor.SourceFormat,
                };
            }

            /// <summary>The filled result turned back upright, with the photo's metadata (as Defisher.Process returns it).</summary>
            public PixelImage Finish(PixelImage r)
            {
                if (Turn != Orientation.Normal) r = Orientation.Apply(r, Orientation.Inverse(Turn));
                r.Metadata = Upright.Metadata != null ? Upright.Metadata.Clone() : new ImageMetadata();
                r.SourceFormat = Upright.SourceFormat;
                return r;
            }

            /// <summary>Rows per dispatch: about 4 million pixels, so each stays well inside the GPU's time limits.</summary>
            public int StripRows => Math.Max(1, Math.Min(OutH, (1 << 22) / Math.Max(1, OutW)));
        }
    }
}
