using System;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Gpu;
using FisheyeFriend.Core.Imaging;
using UnityEngine;
using UnityEngine.Rendering;

namespace FisheyeFriend.Gpu
{
    /// <summary>
    /// The defish on the GPU (FisheyeFriendDefish.compute): the same model and resampler as the CPU code, computed in 32-bit
    /// floats and arranged to keep their precision (see the shader's notes). The CPU-side preparation is
    /// <see cref="GpuPlan"/>; this class uploads it, dispatches the kernels and reads the results back.
    ///
    /// Call it on Unity's main thread. <see cref="ProcessAsync"/> does the CPU work (turning, packing, unpacking) on worker
    /// threads and reads the GPU's results back without blocking the frame.
    /// </summary>
    public sealed class GpuDefisher
    {
        public const string ShaderName = "FisheyeFriendDefish";

        private readonly ComputeShader _cs;
        private readonly int _kMap, _kWarp;

        /// <summary>Rows per Warp dispatch; 0 = the job's own choice (about 4 million pixels). The tests make it small to check the strips join.</summary>
        public int StripRowsOverride { get; set; }

        private GpuDefisher(ComputeShader cs, int kMap, int kWarp)
        {
            _cs = cs;
            _kMap = kMap;
            _kWarp = kWarp;
        }

        /// <summary>The GPU defisher, or null (with the reason) where this computer or build can't run it.</summary>
        public static GpuDefisher Create(out string why)
        {
            why = null;
            if (!SystemInfo.supportsComputeShaders)
            {
                why = "this computer's graphics don't support compute shaders";
                return null;
            }
            var cs = Resources.Load<ComputeShader>(ShaderName);
            if (cs == null)
            {
                why = "the shader " + ShaderName + " is missing from the build";
                return null;
            }
            int kMap, kWarp;
            try
            {
                kMap = cs.FindKernel("MapModel");
                kWarp = cs.FindKernel("Warp");
            }
            catch (ArgumentException)
            {
                why = "the shader " + ShaderName + " has no defish kernels";
                return null;
            }
            if (!cs.IsSupported(kMap) || !cs.IsSupported(kWarp))
            {
                why = "the shader didn't compile for this graphics device (see the console)";
                return null;
            }
            return new GpuDefisher(cs, kMap, kWarp);
        }

        // ================================================================== the model

        /// <summary>The model at the given points: for each (u, v) of the defished frame (interleaved), the photo's (s, t).</summary>
        public float[] MapPoints(DefishModel model, DefishAdjustment adjustment, double aspect, float[] points)
        {
            var p = GpuPlan.ModelParams.For(model, adjustment, aspect);
            SetModel(p);
            int n = points.Length / 2;
            var result = new float[points.Length];
            const int chunk = 64 * 65535;                 // the most points one dispatch can cover
            for (int start = 0; start < n; start += chunk)
            {
                int count = Math.Min(chunk, n - start);
                var part = new float[2 * count];
                Array.Copy(points, 2 * start, part, 0, 2 * count);
                var pts = new ComputeBuffer(count, 8);
                var outB = new ComputeBuffer(count, 8);
                try
                {
                    pts.SetData(part);
                    _cs.SetInt("_PointCount", count);
                    _cs.SetBuffer(_kMap, "_Points", pts);
                    _cs.SetBuffer(_kMap, "_Mapped", outB);
                    _cs.Dispatch(_kMap, (count + 63) / 64, 1, 1);
                    outB.GetData(part);
                    Array.Copy(part, 0, result, 2 * start, 2 * count);
                }
                finally
                {
                    pts.Release();
                    outB.Release();
                }
            }
            return result;
        }

        /// <summary>The geometry with the model evaluated on the GPU at every node (32-bit floats: within about 0.003 px of the CPU's).</summary>
        public DefishGeometry BuildGeometry(DefishModel model, DefishAdjustment adjustment, double aspect, int nodesAcross = 1025)
        {
            var mapped = MapPoints(model, adjustment, aspect, GpuPlan.GridPoints(aspect, nodesAcross));
            return GpuPlan.GeometryFromMapped(model, adjustment, aspect, nodesAcross, mapped);
        }

        private void SetModel(GpuPlan.ModelParams p)
        {
            _cs.SetVectorArray("_Theta", new[] { V(p.Theta, 0), V(p.Theta, 4) });
            _cs.SetVector("_Horizontal", V(p.Horizontal, 0));
            _cs.SetVectorArray("_Vertical", new[] { V(p.Vertical, 0), V(p.Vertical, 4), V(p.Vertical, 8), V(p.Vertical, 12) });
            _cs.SetFloat("_Panini", p.Panini);
            _cs.SetFloat("_PaniniVertical", p.PaniniVertical);
            _cs.SetFloat("_K", p.K);
            _cs.SetFloat("_Strength", p.Strength);
            _cs.SetFloat("_AmountX", p.AmountX);
            _cs.SetFloat("_AmountY", p.AmountY);
            _cs.SetFloat("_Aspect", p.Aspect);
            _cs.SetFloat("_RhoNorm", p.RhoNorm);
            _cs.SetVector("_ShapeX", V(p.ShapeX, 0));
            _cs.SetVector("_ShapeY", V(p.ShapeY, 0));
            _cs.SetInt("_Flat", p.Flat ? 1 : 0);
        }

        private static Vector4 V(float[] a, int o)
        {
            return new Vector4(a[o], a[o + 1], a[o + 2], a[o + 3]);
        }

        // ================================================================== the defish

        /// <summary>Whether this photo and these settings can be defished on the GPU (else the reason).</summary>
        public static bool CanProcess(PixelImage upright, DefishSettings settings, out string why)
        {
            var s = settings.Clone();
            s.Validate();
            if (!GpuPlan.OptionsFit(Defisher.WarpOptionsFor(upright, s), out why)) return false;
            long bytes = upright.BitDepth == 8 ? (long)upright.Data8.Length : 2L * upright.Data16.Length;
            long limit = SystemInfo.maxGraphicsBufferSize;
            if (limit > 0 && bytes > limit)
            {
                why = "the photo (" + (bytes >> 20) + " MB) is larger than the graphics card's largest buffer (" + (limit >> 20) + " MB)";
                return false;
            }
            return true;
        }

        /// <summary>Defishes an upright photo, as <see cref="Defisher.Process"/> does, waiting for the GPU (for tests and tools).</summary>
        public PixelImage Process(PixelImage upright, DefishGeometry g, DefishSettings settings)
        {
            string why;
            if (!CanProcess(upright, settings, out why)) throw new NotSupportedException("Can't defish on the GPU: " + why + ".");
            var job = GpuPlan.Job.Prepare(upright, g, settings);
            var output = job.NewOutput();
            using (var run = new WarpRun(this, job))
            {
                int strip = Strip(job);
                for (int y0 = 0; y0 < job.OutH; y0 += strip)
                {
                    int rows = Math.Min(strip, job.OutH - y0);
                    run.Dispatch(y0, rows);
                    GpuPlan.Unpack(run.Read(rows), output, y0, rows);
                }
            }
            return job.Finish(output);
        }

        /// <summary>
        /// Defishes an upright photo, as <see cref="Defisher.Process"/> does, keeping the window alive: preparation and
        /// unpacking on worker threads; each strip is read back as soon as the GPU has it, and a frame is let through about
        /// five times a second (for the progress bar). It never waits a frame per strip, so the app's frame rate (lowered
        /// to save battery, see FramePacer) doesn't slow it down. Start it on the main thread.
        /// </summary>
        public async Task<PixelImage> ProcessAsync(PixelImage upright, DefishGeometry g, DefishSettings settings,
            ProgressTracker progress = null, CancellationToken ct = default)
        {
            string why;
            if (!CanProcess(upright, settings, out why)) throw new NotSupportedException("Can't defish on the GPU: " + why + ".");
            var job = await Task.Run(() => GpuPlan.Job.Prepare(upright, g, settings), ct);
            ct.ThrowIfCancellationRequested();
            var output = job.NewOutput();
            var unpacking = new System.Collections.Generic.List<Task>();
            using (var run = new WarpRun(this, job))
            {
                int strip = Strip(job);
                var sinceFrame = System.Diagnostics.Stopwatch.StartNew();
                for (int y0 = 0; y0 < job.OutH; y0 += strip)
                {
                    ct.ThrowIfCancellationRequested();
                    int rows = Math.Min(strip, job.OutH - y0), start = y0;
                    run.Dispatch(y0, rows);
                    var words = run.Read(rows);                   // waits for this strip on the GPU (tens of milliseconds)
                    unpacking.Add(Task.Run(() => GpuPlan.Unpack(words, output, start, rows), ct));
                    if (progress != null) progress.Report((float)(y0 + rows) / job.OutH);
                    if (sinceFrame.ElapsedMilliseconds > 200 && y0 + rows < job.OutH)
                    {
                        await Task.Yield();                       // one frame, so the window and the progress bar update
                        sinceFrame.Restart();
                    }
                }
            }
            await Task.WhenAll(unpacking);
            ct.ThrowIfCancellationRequested();
            return await Task.Run(() => job.Finish(output), ct);
        }

        private int Strip(GpuPlan.Job job)
        {
            return StripRowsOverride > 0 ? Math.Min(StripRowsOverride, job.OutH) : job.StripRows;
        }

        /// <summary>The buffers of one defish on the GPU.</summary>
        private sealed class WarpRun : IDisposable
        {
            private readonly GpuDefisher _d;
            private readonly GpuPlan.Job _job;
            private readonly ComputeBuffer _src, _nodes, _cols, _rows, _colorLut, _weightLut, _apTable, _fromSig, _encode, _out;
            private readonly int _wpp;

            public WarpRun(GpuDefisher d, GpuPlan.Job job)
            {
                _d = d;
                _job = job;
                var p = job.Params;
                var gr = job.Grid;
                _wpp = p.WordsPerPixel;
                try
                {
                    _src = Upload(job.Samples, 4);
                    _nodes = Upload(gr.Nodes, Marshal.SizeOf(typeof(GpuPlan.Node)));
                    _cols = Upload(gr.Cols, Marshal.SizeOf(typeof(GpuPlan.Cell)));
                    _rows = Upload(gr.Rows, Marshal.SizeOf(typeof(GpuPlan.Cell)));
                    _colorLut = Upload(p.ColorLut, 4);
                    _weightLut = Upload(p.WeightLut, 4);
                    _apTable = Upload(p.ApTable, 4);
                    _fromSig = Upload(p.FromSig, 4);
                    _encode = Upload(p.Encode, 4);
                    int strip = d.Strip(job);
                    _out = new ComputeBuffer(Math.Max(1, strip * job.OutW * _wpp), 4);
                }
                catch
                {
                    Dispose();
                    throw;
                }
            }

            private static ComputeBuffer Upload(Array data, int stride)
            {
                var b = new ComputeBuffer(Math.Max(1, data.Length), stride);
                if (data.Length > 0) b.SetData(data);
                return b;
            }

            public void Dispatch(int y0, int rows)
            {
                var cs = _d._cs;
                int k = _d._kWarp;
                var p = _job.Params;
                var gr = _job.Grid;
                cs.SetBuffer(k, "_Src", _src);
                cs.SetBuffer(k, "_Nodes", _nodes);
                cs.SetBuffer(k, "_Cols", _cols);
                cs.SetBuffer(k, "_Rows", _rows);
                cs.SetBuffer(k, "_ColorLut", _colorLut);
                cs.SetBuffer(k, "_WeightLut", _weightLut);
                cs.SetBuffer(k, "_ApTable", _apTable);
                cs.SetBuffer(k, "_FromSig", _fromSig);
                cs.SetBuffer(k, "_Encode", _encode);
                cs.SetBuffer(k, "_Out", _out);
                cs.SetInt("_NodesU", gr.NodesU);
                cs.SetFloat("_ColStep", gr.ColStep);
                cs.SetFloat("_RowStep", gr.RowStep);
                cs.SetInt("_SrcW", _job.Sensor.Width);
                cs.SetInt("_SrcH", _job.Sensor.Height);
                cs.SetInt("_Ch", p.Channels);
                cs.SetInt("_Depth", p.Depth);
                cs.SetInt("_Alpha", p.Alpha ? 1 : 0);
                cs.SetFloat("_AlphaScale", p.AlphaScale);
                cs.SetInt("_WeightLutSize", p.WeightLutSize);
                cs.SetFloat("_R", p.R);
                cs.SetFloat("_LutScale", p.LutScale);
                cs.SetInt("_Antialias", p.Antialias ? 1 : 0);
                cs.SetFloat("_AaStrength", p.AaStrength);
                cs.SetFloat("_Aperture", p.Aperture);
                cs.SetInt("_SharpenOnly", p.SharpenOnly ? 1 : 0);
                cs.SetInt("_ApCount", p.ApCount);
                cs.SetFloat("_ApLogMin", p.ApLogMin);
                cs.SetFloat("_ApLogMax", p.ApLogMax);
                cs.SetInt("_Extend", p.Extend ? 1 : 0);
                cs.SetVector("_Bg", new Vector4(p.Bg[0], p.Bg[1], p.Bg[2], p.Bg[3]));
                cs.SetInt("_Light", p.Light);
                cs.SetInt("_OutW", _job.OutW);
                cs.SetInt("_OutH", _job.OutH);
                cs.SetInt("_StripY0", y0);
                cs.SetInt("_StripRows", rows);
                cs.SetInt("_OutDepth", p.OutDepth);
                cs.SetFloat("_OutMax", p.OutMax);
                cs.Dispatch(k, (_job.OutW + 7) / 8, (rows + 7) / 8, 1);
            }

            public uint[] Read(int rows)
            {
                var r = new uint[rows * _job.OutW * _wpp];
                _out.GetData(r, 0, 0, r.Length);
                return r;
            }

            public void Dispose()
            {
                foreach (var b in new[] { _src, _nodes, _cols, _rows, _colorLut, _weightLut, _apTable, _fromSig, _encode, _out })
                    if (b != null) b.Release();
            }
        }
    }
}
