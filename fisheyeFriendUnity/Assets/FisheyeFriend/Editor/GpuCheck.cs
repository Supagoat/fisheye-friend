using System;
using System.Diagnostics;
using System.IO;
using System.Text;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Testing;
using FisheyeFriend.Gpu;
using UnityEditor;
using UnityEngine;

namespace FisheyeFriend.EditorTools
{
    /// <summary>
    /// Fisheye Friend > Check GPU Against CPU: the GPU conformance checks of ConformanceTests in one go, plus a timing on a
    /// full-size photo, written to Tools/TestData/gpu-report.txt (and summarised in a dialog). Handy where running the
    /// Test Runner and copying its output is a chore.
    /// </summary>
    internal static class GpuCheck
    {
        private static StringBuilder _log;
        private static int _fails;

        private static void Line(string s)
        {
            _log.AppendLine(s);
        }

        private static void Check(bool ok, string what)
        {
            Line((ok ? "PASS  " : "FAIL  ") + what);
            if (!ok) _fails++;
        }

        [MenuItem("Fisheye Friend/Check GPU Against CPU")]
        private static void Run()
        {
            _log = new StringBuilder();
            _fails = 0;
            string dataDir = Path.GetFullPath(Path.Combine(Application.dataPath, "..", "Tools", "TestData"));
            string reportPath = Path.Combine(dataDir, "gpu-report.txt");
            Line("Fisheye Friend GPU check, " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"));
            Line("Unity " + Application.unityVersion + ", " + SystemInfo.operatingSystem);
            Line("GPU: " + SystemInfo.graphicsDeviceName + " (" + SystemInfo.graphicsDeviceType + ", " + SystemInfo.graphicsMemorySize + " MB); " +
                 "largest buffer " + (SystemInfo.maxGraphicsBufferSize >> 20) + " MB; CPU: " + SystemInfo.processorType + " x" + SystemInfo.processorCount);
            try
            {
                var cs = Resources.Load<ComputeShader>(GpuDefisher.ShaderName);
                if (cs != null)
                {
                    var messages = ShaderUtil.GetComputeShaderMessages(cs);
                    Line("Shader messages: " + messages.Length);
                    foreach (var m in messages) Line("  " + m.severity + " line " + m.line + ": " + m.message + " (" + m.platform + ")");
                }
                RunChecks(dataDir);
            }
            catch (Exception e)
            {
                Check(false, "the check stopped: " + e);
            }
            finally
            {
                EditorUtility.ClearProgressBar();
            }
            Line(_fails == 0 ? "ALL PASSED" : _fails + " FAILED");
            File.WriteAllText(reportPath, _log.ToString());
            UnityEngine.Debug.Log("[Fisheye Friend] GPU check: " + (_fails == 0 ? "all passed" : _fails + " failed") + ". Report: " + reportPath + "\n" + _log);
            EditorUtility.DisplayDialog("GPU check", (_fails == 0 ? "All checks passed." : _fails + " checks failed.") + "\n\nThe report is in Tools/TestData/gpu-report.txt.", "OK");
        }

        private static void RunChecks(string dataDir)
        {
            string why;
            var gpu = GpuDefisher.Create(out why);
            if (gpu == null)
            {
                Check(false, "no GPU defish: " + why);
                return;
            }
            var model = DefishModel.Parse(File.ReadAllText(Path.Combine(Application.streamingAssetsPath, "defish-hemi-model.json")));
            var golden = Conformance.LoadGolden(File.ReadAllText(Path.Combine(dataDir, "defish-golden.json")));

            // the model kernel: golden points, and a whole grid against the CPU's
            double worstPts = 0, worstNodes = 0;
            bool foldsAgree = true;
            for (int ci = 0; ci < golden.Count; ci++)
            {
                var c = golden[ci];
                EditorUtility.DisplayProgressBar("GPU check", "Model: " + c.Name, 0.3f * ci / golden.Count);
                var pts = new float[2 * c.Points.Length];
                for (int i = 0; i < c.Points.Length; i++)
                {
                    pts[2 * i] = (float)c.Points[i][0];
                    pts[2 * i + 1] = (float)c.Points[i][1];
                }
                var mapped = gpu.MapPoints(model, c.Adjustment, c.Aspect, pts);
                var byPoint = new System.Collections.Generic.Dictionary<(double, double), double[]>();
                for (int i = 0; i < c.Points.Length; i++) byPoint[(c.Points[i][0], c.Points[i][1])] = new double[] { mapped[2 * i], mapped[2 * i + 1] };
                var e = Conformance.CompareMapping(c, (u, v) => byPoint[(u, v)]);
                var gg = gpu.BuildGeometry(model, c.Adjustment, c.Aspect);
                var gc = DefishGeometry.FromModel(model, c.Adjustment, c.Aspect);
                double nd = 0;
                for (int i = 0; i < gc.GridNodes; i++)
                {
                    float s1, t1, s2, t2;
                    gg.Node(i, out s1, out t1);
                    gc.Node(i, out s2, out t2);
                    nd = Math.Max(nd, Math.Max(Math.Abs(s1 - s2), Math.Abs(t1 - t2) * c.Aspect) * Conformance.ReferenceWidth);
                }
                foldsAgree &= gg.Folds == gc.Folds;
                Line("  " + c.Name.PadRight(32) + " golden points: " + e + "   grid nodes vs CPU: max " + nd.ToString("0.00000") + " px");
                worstPts = Math.Max(worstPts, e.MaxPx);
                worstNodes = Math.Max(worstNodes, nd);
            }
            Check(worstPts < 0.01, "GPU model vs golden points: max " + worstPts.ToString("0.00000") + " px on an " + Conformance.ReferenceWidth + " px frame (need < 0.01)");
            Check(worstNodes < 0.01, "GPU grid vs CPU grid: max " + worstNodes.ToString("0.00000") + " px (need < 0.01)");
            Check(foldsAgree, "GPU and CPU grids agree on folding");

            // the resampler: every render case against the CPU (same geometry) and the golden render
            var cases = Conformance.RenderCases;
            for (int k = 0; k < cases.Length; k++)
            {
                var rc = cases[k];
                EditorUtility.DisplayProgressBar("GPU check", "Render: " + rc.Id, 0.3f + 0.4f * k / cases.Length);
                var img = rc.MakeImage();
                var g = DefishGeometry.FromModel(model, rc.Settings.Adjust, rc.SensorAspect);
                var cpu = Defisher.Process(img, g, rc.Settings);
                var r = gpu.Process(img, g, rc.Settings);
                var d = Conformance.Compare(r, cpu);
                var dg = Conformance.Compare(r, PngCodec.Decode(File.ReadAllBytes(Path.Combine(dataDir, rc.GoldenFile))));
                var full = gpu.Process(img, gpu.BuildGeometry(model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
                var df = Conformance.Compare(full, cpu);
                Line("  " + rc.Id);
                Line("      GPU vs CPU (same geometry):   " + d);
                Line("      GPU vs golden render:         " + dg);
                Line("      all on the GPU vs CPU:        " + df);
                bool exact = rc.Id == "rgb8-none";
                double need = rc.Depth == 8 ? 0.999 : 0.99;
                Check(d.SameLayout && d.Max <= (exact ? 0 : 1) && d.IdenticalFraction >= (exact ? 1 : need),
                    rc.Id + ": GPU matches the CPU (max " + d.Max + ", " + (d.IdenticalFraction * 100).ToString("0.000") + "% identical)");
                Check(dg.SameLayout && dg.Max <= (exact ? 0 : 1) && dg.IdenticalFraction >= (exact ? 1 : need),
                    rc.Id + ": GPU matches the golden render (max " + dg.Max + ", " + (dg.IdenticalFraction * 100).ToString("0.000") + "% identical)");
            }

            // strips join seamlessly
            {
                var rc = Array.Find(cases, x => x.Id == "rgb8-adjusted-crop");
                var img = rc.MakeImage();
                var g = DefishGeometry.FromModel(model, rc.Settings.Adjust, rc.SensorAspect);
                var whole = gpu.Process(img, g, rc.Settings);
                gpu.StripRowsOverride = 7;
                var strips = gpu.Process(img, g, rc.Settings);
                gpu.StripRowsOverride = 0;
                var d = Conformance.Compare(strips, whole);
                Check(d.Max == 0, "strips of 7 rows equal one dispatch: " + d);
            }

            // speed: a full-size photo, CPU against GPU
            {
                EditorUtility.DisplayProgressBar("GPU check", "Timing a 24 MP photo (CPU, then GPU)", 0.75f);
                var img = Conformance.Scene(6000, 4000, 3, 8, 3);
                var s = Array.Find(cases, x => x.Id == "rgb8-default").Settings;
                var g = DefishGeometry.FromModel(model, s.Adjust, 4000.0 / 6000.0);
                var sw = Stopwatch.StartNew();
                var cpu = Defisher.Process(img, g, s);
                double tCpu = sw.Elapsed.TotalSeconds;
                gpu.Process(img, g, s);                              // the first run sets the shader up
                sw.Restart();
                var r = gpu.Process(img, g, s);
                double tGpu = sw.Elapsed.TotalSeconds;
                var d = Conformance.Compare(r, cpu);
                Line("  24 MP (6000 x 4000) RGB, defaults: CPU " + tCpu.ToString("0.00") + " s, GPU " + tGpu.ToString("0.00") + " s (" +
                     (tCpu / Math.Max(1e-6, tGpu)).ToString("0.0") + "x), output " + r.Width + " x " + r.Height);
                Check(d.SameLayout && d.Max <= 1 && d.IdenticalFraction >= 0.999, "24 MP photo: GPU matches the CPU: " + d);
            }
        }
    }
}
