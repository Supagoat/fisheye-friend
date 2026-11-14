using System.Collections.Generic;
using System.IO;
using System.Linq;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Gpu;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Testing;
using FisheyeFriend.Gpu;
using NUnit.Framework;
using UnityEngine;

namespace FisheyeFriend.Tests
{
    /// <summary>
    /// The defish conformance suite (see FisheyeFriend.Core.Testing.Conformance), run against every implementation.
    ///
    /// <b>CPU</b>: the model against golden points computed independently (Tools/TestData/golden.py, double precision),
    /// and the full defish against golden renders of test images, at many adjustments and every option.
    ///
    /// <b>GPU</b>: the same golden points through the shader's model kernel, and the shader's resampler against the CPU
    /// defish and the golden renders. The GPU computes in 32-bit floats, so it is held to: positions within 1/100 px,
    /// and pixels equal except where a value sits on a rounding boundary, then by one level at most.
    ///
    /// Each test prints what it measured (Test Runner: select the test, see the output below).
    /// </summary>
    public class ConformanceTests
    {
        private static string DataDir => Path.GetFullPath(Path.Combine(Application.dataPath, "..", "Tools", "TestData"));

        private static DefishModel Model => DefishModel.Parse(File.ReadAllText(Path.Combine(Application.streamingAssetsPath, "defish-hemi-model.json")));

        private static List<Conformance.GoldenCase> _golden;

        private static List<Conformance.GoldenCase> Golden
        {
            get
            {
                if (_golden == null) _golden = Conformance.LoadGolden(File.ReadAllText(Path.Combine(DataDir, "defish-golden.json")));
                return _golden;
            }
        }

        public static IEnumerable<string> RenderIds()
        {
            return Conformance.RenderCases.Select(r => r.Id);
        }

        private static Conformance.RenderCase Case(string id)
        {
            return Conformance.RenderCases.First(r => r.Id == id);
        }

        private static PixelImage GoldenRender(Conformance.RenderCase rc)
        {
            return PngCodec.Decode(File.ReadAllBytes(Path.Combine(DataDir, rc.GoldenFile)));
        }

        private static void Report(string line)
        {
            TestContext.WriteLine(line);
        }

        // ============================================================================ step 1: the CPU defish

        [Test]
        public void Cpu_ModelMatchesGoldenPoints()
        {
            var model = Model;
            double worst = 0;
            foreach (var c in Golden)
            {
                var m = new DefishModelMap(model, c.Adjustment, c.Aspect);
                var e = Conformance.CompareMapping(c, (u, v) =>
                {
                    double s, t;
                    m.Reset();
                    m.Map(u, v, out s, out t);
                    return new[] { s, t };
                });
                Report(c.Name + ": " + e);
                worst = System.Math.Max(worst, e.MaxPx);
            }
            Assert.Less(worst, 1e-4, "px on an " + Conformance.ReferenceWidth + " px frame");
        }

        [Test]
        public void Cpu_GeometryMatchesGoldenPoints()
        {
            var model = Model;
            double worst = 0;
            foreach (var c in Golden)
            {
                var g = DefishGeometry.FromModel(model, c.Adjustment, c.Aspect);
                var e = Conformance.CompareMapping(c, (u, v) =>
                {
                    double s, t, a, b, cc, d;
                    g.Map(u, v, out s, out t, out a, out b, out cc, out d);
                    return new[] { s, t };
                });
                Report(c.Name + ": " + e);
                worst = System.Math.Max(worst, e.MaxPx);
            }
            Assert.Less(worst, 0.002, "px on an " + Conformance.ReferenceWidth + " px frame");
        }

        [TestCaseSource(nameof(RenderIds))]
        public void Cpu_RenderMatchesGolden(string id)
        {
            var rc = Case(id);
            var r = Defisher.Process(rc.MakeImage(), DefishGeometry.FromModel(Model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
            var d = Conformance.Compare(r, GoldenRender(rc));
            Report(rc.Description + "\n" + d);
            Assert.IsTrue(d.SameLayout, d.ToString());
            // The golden renders were made on another computer: its maths library may round the odd value differently.
            Assert.LessOrEqual(d.Max, 1, d.ToString());
            Assert.GreaterOrEqual(d.IdenticalFraction, 0.9999, d.ToString());
        }

        [Test]
        public void Cpu_NoCorrectionGivesThePhotoBack()
        {
            var rc = Case("rgb8-none");
            var img = rc.MakeImage();
            var r = Defisher.Process(img, DefishGeometry.FromModel(Model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
            var d = Conformance.Compare(r, img);
            Report(d.ToString());
            Assert.AreEqual(0, d.Max, d.ToString());
        }

        // ============================================================================ step 2: the GPU defish

        private static GpuDefisher _gpu;
        private static string _gpuWhy;

        private static GpuDefisher Gpu
        {
            get
            {
                if (_gpu == null && _gpuWhy == null) _gpu = GpuDefisher.Create(out _gpuWhy);
                if (_gpu == null) Assert.Ignore("No GPU defish here: " + _gpuWhy);
                _gpu.StripRowsOverride = 0;
                return _gpu;
            }
        }

        [Test]
        public void Gpu_ModelMatchesGoldenPoints()
        {
            var gpu = Gpu;
            var model = Model;
            double worst = 0;
            foreach (var c in Golden)
            {
                var pts = new float[2 * c.Points.Length];
                for (int i = 0; i < c.Points.Length; i++)
                {
                    pts[2 * i] = (float)c.Points[i][0];
                    pts[2 * i + 1] = (float)c.Points[i][1];
                }
                var mapped = gpu.MapPoints(model, c.Adjustment, c.Aspect, pts);
                var byPoint = new Dictionary<(double, double), double[]>();
                for (int i = 0; i < c.Points.Length; i++)
                    byPoint[(c.Points[i][0], c.Points[i][1])] = new double[] { mapped[2 * i], mapped[2 * i + 1] };
                var e = Conformance.CompareMapping(c, (u, v) => byPoint[(u, v)]);
                Report(c.Name + ": " + e);
                worst = System.Math.Max(worst, e.MaxPx);
            }
            Assert.Less(worst, 0.01, "px on an " + Conformance.ReferenceWidth + " px frame");
        }

        [Test]
        public void Gpu_GeometryMatchesCpuGeometry()
        {
            var gpu = Gpu;
            var model = Model;
            double worst = 0;
            foreach (var c in Golden)
            {
                var gg = gpu.BuildGeometry(model, c.Adjustment, c.Aspect);
                var gc = DefishGeometry.FromModel(model, c.Adjustment, c.Aspect);
                double d = 0;
                for (int i = 0; i < gc.GridNodes; i++)
                {
                    float s1, t1, s2, t2;
                    gg.Node(i, out s1, out t1);
                    gc.Node(i, out s2, out t2);
                    d = System.Math.Max(d, System.Math.Max(System.Math.Abs(s1 - s2), System.Math.Abs(t1 - t2) * c.Aspect) * Conformance.ReferenceWidth);
                }
                Report(c.Name + ": nodes within " + d.ToString("0.00000") + " px of the CPU's (" + gc.GridNodes + " nodes)");
                worst = System.Math.Max(worst, d);
                Assert.AreEqual(gc.Folds, gg.Folds, c.Name + ": fold check");
            }
            Assert.Less(worst, 0.01, "px on an " + Conformance.ReferenceWidth + " px frame");
        }

        /// <summary>The GPU's resampling against the CPU's, with the same geometry: equal but for rounding boundaries.</summary>
        [TestCaseSource(nameof(RenderIds))]
        public void Gpu_RenderMatchesCpu(string id)
        {
            var gpu = Gpu;
            var rc = Case(id);
            var img = rc.MakeImage();
            var g = DefishGeometry.FromModel(Model, rc.Settings.Adjust, rc.SensorAspect);
            var cpu = Defisher.Process(img, g, rc.Settings);
            var r = gpu.Process(img, g, rc.Settings);
            var d = Conformance.Compare(r, cpu);
            Report(rc.Description + "\nGPU vs CPU: " + d);
            Assert.IsTrue(d.SameLayout, d.ToString());
            if (id == "rgb8-none")
            {
                Assert.AreEqual(0, d.Max, "with no correction the GPU must give the photo back exactly: " + d);
                return;
            }
            Assert.LessOrEqual(d.Max, 1, d.ToString());
            Assert.GreaterOrEqual(d.IdenticalFraction, rc.Depth == 8 ? 0.999 : 0.99, d.ToString());
        }

        [TestCaseSource(nameof(RenderIds))]
        public void Gpu_RenderMatchesGolden(string id)
        {
            var gpu = Gpu;
            var rc = Case(id);
            var r = gpu.Process(rc.MakeImage(), DefishGeometry.FromModel(Model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
            var d = Conformance.Compare(r, GoldenRender(rc));
            Report(rc.Description + "\nGPU vs golden: " + d);
            Assert.IsTrue(d.SameLayout, d.ToString());
            Assert.LessOrEqual(d.Max, 1, d.ToString());
            Assert.GreaterOrEqual(d.IdenticalFraction, rc.Depth == 8 ? 0.999 : 0.99, d.ToString());
        }

        /// <summary>The result doesn't depend on how the output is split into dispatches.</summary>
        [Test]
        public void Gpu_StripsJoinSeamlessly()
        {
            var gpu = Gpu;
            var rc = Case("rgb8-adjusted-crop");
            var img = rc.MakeImage();
            var g = DefishGeometry.FromModel(Model, rc.Settings.Adjust, rc.SensorAspect);
            var whole = gpu.Process(img, g, rc.Settings);
            try
            {
                gpu.StripRowsOverride = 7;
                var strips = gpu.Process(img, g, rc.Settings);
                var d = Conformance.Compare(strips, whole);
                Report("in strips of 7 rows vs one dispatch: " + d);
                Assert.AreEqual(0, d.Max, d.ToString());
            }
            finally
            {
                gpu.StripRowsOverride = 0;
            }
        }

        /// <summary>Everything on the GPU (its own geometry too), for the record: within a level or so of the CPU.</summary>
        [TestCaseSource(nameof(RenderIds))]
        public void Gpu_GeometryAndRenderMatchCpu(string id)
        {
            var gpu = Gpu;
            var rc = Case(id);
            var img = rc.MakeImage();
            var model = Model;
            var cpu = Defisher.Process(img, DefishGeometry.FromModel(model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
            var r = gpu.Process(img, gpu.BuildGeometry(model, rc.Settings.Adjust, rc.SensorAspect), rc.Settings);
            var d = Conformance.Compare(r, cpu);
            Report(rc.Description + "\nGPU geometry and resampling vs CPU: " + d);
            Assert.IsTrue(d.SameLayout, d.ToString());
            Assert.LessOrEqual(d.Max, rc.Depth == 8 ? 1 : 8, d.ToString());
        }
    }
}
