using System;
using System.IO;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Core.Resampling;
using NUnit.Framework;
using UnityEngine;

namespace FisheyeFriend.Tests
{
    /// <summary>EditMode tests for the core library (Window > General > Test Runner).</summary>
    public class CoreTests
    {
        private static DefishGeometry _geometry;

        private static DefishGeometry Geometry
        {
            get
            {
                if (_geometry == null)
                {
                    string path = Path.Combine(Application.streamingAssetsPath, "defish-hemi-596.json");
                    _geometry = new DefishGeometry(HemiMapping.Parse(File.ReadAllText(path)));
                }
                return _geometry;
            }
        }

        private sealed class IdentityMap : WarpMap
        {
            public override void Map(double x, double y, out double sx, out double sy, out double j00, out double j01, out double j10, out double j11)
            {
                sx = x; sy = y; j00 = 1; j01 = 0; j10 = 0; j11 = 1;
            }
        }

        private static PixelImage Noise(int w, int h, int ch, int depth, int seed)
        {
            var img = new PixelImage(w, h, ch, depth);
            var rnd = new System.Random(seed);
            if (depth == 8) rnd.NextBytes(img.Data8);
            else for (int i = 0; i < img.Data16.Length; i++) img.Data16[i] = (ushort)rnd.Next(65536);
            if (ch == 4)
                for (int i = 3; i < w * h * 4; i += 4)
                    if (depth == 8) img.Data8[i] = 255; else img.Data16[i] = 65535;
            return img;
        }

        private static int MaxDiff(PixelImage a, PixelImage b, int channels)
        {
            int m = 0;
            for (int i = 0; i < a.Width * a.Height; i++)
                for (int c = 0; c < channels; c++)
                {
                    int va = a.BitDepth == 8 ? a.Data8[i * a.Channels + c] : a.Data16[i * a.Channels + c];
                    int vb = b.BitDepth == 8 ? b.Data8[i * b.Channels + c] : b.Data16[i * b.Channels + c];
                    m = Math.Max(m, Math.Abs(va - vb));
                }
            return m;
        }

        [Test]
        public void Mapping_LoadsTheHemiPairs()
        {
            var m = Geometry.Mapping;
            Assert.AreEqual(596, m.Pairs.Count);
            Assert.AreEqual(6000, m.SourceWidth);
            Assert.AreEqual(4000, m.SourceHeight);
        }

        [Test]
        public void Mapping_ParsesPixelUnitsAndReportsErrors()
        {
            var alt = HemiMapping.Parse("{\"source\":{\"width\":100,\"height\":50},\"target\":{\"width\":100,\"height\":50},\"pairs\":[" +
                                        "{\"src\":{\"x\":10,\"y\":5},\"dst\":{\"xn\":0.1,\"yn\":0.1}},{\"src\":{\"x\":90,\"y\":5},\"dst\":{\"x\":90,\"y\":5}}," +
                                        "{\"src\":{\"x\":50,\"y\":45},\"dst\":{\"x\":50,\"y\":45}}]}");
            Assert.AreEqual(3, alt.Pairs.Count);
            Assert.AreEqual(0.1, alt.Pairs[0].SrcX, 1e-12);
            Assert.AreEqual(0.9, alt.Pairs[2].DstY, 1e-12);
            Assert.Throws<FormatException>(() => HemiMapping.Parse("{\"pairs\": [1, 2"));
        }

        [Test]
        public void Geometry_GridMatchesTheExactSplines()
        {
            Assert.Less(Geometry.MaxGridError(7952, 5304, 5000, 3), 0.01);
        }

        [Test]
        public void Geometry_UnmapInvertsMap()
        {
            var rnd = new System.Random(5);
            for (int k = 0; k < 500; k++)
            {
                double u = 0.02 + 0.96 * rnd.NextDouble(), v = 0.02 + 0.96 * rnd.NextDouble(), s, t, a, b, c, d, uu, vv;
                Geometry.Map(u, v, out s, out t, out a, out b, out c, out d);
                Assert.IsTrue(Geometry.Unmap(s, t, out uu, out vv));
                Assert.AreEqual(u, uu, 1e-8);
                Assert.AreEqual(v, vv, 1e-8);
            }
        }

        [Test]
        public void Warp_IdentityIsExact()
        {
            foreach (var depth in new[] { 8, 16 })
                foreach (var ch in new[] { 1, 3, 4 })
                {
                    var img = Noise(37, 23, ch, depth, ch * depth);
                    foreach (var k in new[] { WarpKernel.TensorLanczos, WarpKernel.Bicubic })
                        foreach (var l in new[] { LightSpace.Gamma, LightSpace.Sigmoidal, LightSpace.Linear })
                        {
                            var r = WarpResampler.Warp(img, 37, 23, new IdentityMap(), new WarpOptions
                            {
                                Kernel = k, Radius = k == WarpKernel.Bicubic ? 2 : 4, Light = l, Aperture = 1, ApertureSharpenOnly = true,
                                Antialias = false, ExtendEdges = true,
                            });
                            Assert.LessOrEqual(MaxDiff(img, r, ch), depth == 8 ? 0 : 1, ch + " ch " + depth + "-bit " + k + " " + l);
                        }
                }
        }

        [Test]
        public void Warp_ApertureCoefficientFollowsTheScale()
        {
            Assert.AreEqual(0, WarpResampler.ApertureCoefficient(1), 1e-6);
            Assert.Greater(WarpResampler.ApertureCoefficient(0.55), 0.05);
            Assert.Less(WarpResampler.ApertureCoefficient(1.19), 0);
        }

        [Test]
        public void Defish_FlatFieldStaysFlat()
        {
            foreach (var l in new[] { LightSpace.Gamma, LightSpace.Sigmoidal, LightSpace.Linear })
            {
                var flat = new PixelImage(150, 100, 3, 8);
                for (int i = 0; i < flat.Data8.Length; i += 3) { flat.Data8[i] = 200; flat.Data8[i + 1] = 90; flat.Data8[i + 2] = 17; }
                var r = Defisher.Process(flat, Geometry, new DefishSettings { Light = l });
                for (int i = 0; i < r.Data8.Length; i += 3)
                {
                    Assert.LessOrEqual(Math.Abs(r.Data8[i] - 200), 1);
                    Assert.LessOrEqual(Math.Abs(r.Data8[i + 1] - 90), 1);
                    Assert.LessOrEqual(Math.Abs(r.Data8[i + 2] - 17), 1);
                }
            }
        }

        [Test]
        public void Defish_RgbAndRgbaAgree()
        {
            var rgb = Noise(160, 107, 3, 8, 11);
            var rgba = new PixelImage(160, 107, 4, 8);
            for (int i = 0; i < 160 * 107; i++)
            {
                for (int c = 0; c < 3; c++) rgba.Data8[i * 4 + c] = rgb.Data8[i * 3 + c];
                rgba.Data8[i * 4 + 3] = 255;
            }
            var s = new DefishSettings { Light = LightSpace.Gamma };
            Assert.LessOrEqual(MaxDiff(Defisher.Process(rgb, Geometry, s), Defisher.Process(rgba, Geometry, s), 3), 1);
        }

        [Test]
        public void Defish_OutputSizes()
        {
            var land = Noise(300, 200, 3, 8, 21);
            int w, h;
            Defisher.OutputSize(Geometry, land, new DefishSettings(), out w, out h);
            Assert.AreEqual(300, w);
            Assert.AreEqual(200, h);
            double k = Defisher.CentreScale(Geometry, 300, 200);
            Assert.AreEqual(1.19, k, 0.01);
            Defisher.OutputSize(Geometry, land, new DefishSettings { OutputSize = DefishOutputSize.FullDetail }, out w, out h);
            Assert.AreEqual((int)Math.Round(300 * k), w);
            var c = Defisher.CropCanvas(Geometry, 7952, 5304);
            Assert.Greater(c.Width, 0.995);
            Assert.Less(c.Width, 1);
        }

        [Test]
        public void Defish_PortraitsUseTheSensorFrame()
        {
            var land = Noise(240, 160, 3, 8, 31);
            var s = new DefishSettings();
            var rl = Defisher.Process(land, Geometry, s);
            foreach (var o in new[] { Orientation.Rotate90Clockwise, Orientation.Rotate90CounterClockwise })
            {
                var portrait = Orientation.Apply(land, o);
                portrait.Metadata.SourceOrientation = o;
                var rp = Defisher.Process(portrait, Geometry, s);
                var expect = Orientation.Apply(rl, o);
                Assert.AreEqual(expect.Width, rp.Width);
                Assert.AreEqual(expect.Height, rp.Height);
                Assert.AreEqual(0, MaxDiff(rp, expect, 3));
            }
        }

        [Test]
        public void Defish_KeepsBitDepthAndMetadata()
        {
            var img = Noise(90, 60, 3, 16, 41);
            img.Metadata.IccProfile = new byte[] { 1, 2, 3 };
            var r = Defisher.Process(img, Geometry, new DefishSettings());
            Assert.AreEqual(16, r.BitDepth);
            CollectionAssert.AreEqual(img.Metadata.IccProfile, r.Metadata.IccProfile);
        }

        [Test]
        public void Defish_ViewMappingRoundTrips()
        {
            var rnd = new System.Random(9);
            foreach (var o in new[] { 0, Orientation.Rotate90Clockwise })
                foreach (var os in new[] { DefishOutputSize.SameAsPhoto, DefishOutputSize.FullDetail })
                {
                    var img = o == 0 ? new PixelImage(600, 400, 3, 8) : new PixelImage(400, 600, 3, 8);
                    img.Metadata.SourceOrientation = o;
                    var s = new DefishSettings { OutputSize = os, Edges = DefishEdges.Crop };
                    for (int k = 0; k < 50; k++)
                    {
                        double x = 20 + rnd.NextDouble() * (img.Width - 40), y = 20 + rnd.NextDouble() * (img.Height - 40), rx, ry, bx, by;
                        if (!Defisher.PhotoToResult(Geometry, img, s, x, y, out rx, out ry)) continue;
                        Defisher.ResultToPhoto(Geometry, img, s, rx, ry, out bx, out by);
                        Assert.AreEqual(x, bx, 1e-4);
                        Assert.AreEqual(y, by, 1e-4);
                    }
                }
        }

        [Test]
        public void Png_RoundTripsAllLayouts()
        {
            var rnd = new System.Random(3);
            foreach (var ch in new[] { 1, 2, 3, 4 })
                foreach (var depth in new[] { 8, 16 })
                {
                    var img = new PixelImage(31, 17, ch, depth);
                    if (depth == 8) rnd.NextBytes(img.Data8);
                    else for (int i = 0; i < img.Data16.Length; i++) img.Data16[i] = (ushort)rnd.Next(65536);
                    var back = PngCodec.Decode(PngCodec.Encode(img));
                    Assert.AreEqual(ch, back.Channels);
                    Assert.AreEqual(depth, back.BitDepth);
                    if (depth == 8) CollectionAssert.AreEqual(img.Data8, back.Data8);
                    else CollectionAssert.AreEqual(img.Data16, back.Data16);
                }
        }
    }
}
