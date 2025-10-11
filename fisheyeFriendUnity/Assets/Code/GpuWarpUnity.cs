using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Drawing;
using System.Drawing.Imaging;

namespace FisheyeFriend {
    // Example integration in your existing flow (Unity-only):
    public class GpuWarpUnity : MonoBehaviour {
        public ComputeShader tpsWarpCS;
        public Texture2D sourceTexture;
        public RenderTexture outputRT;

        FisheyeFriend.TPS2D tps;
        GpuWarper warper;

        void Start () {

            Debug.Log("Starting GPU VERSION");
            // After you've built TPS in C# from your mapping JSON (unchanged from your code)
            //_kernel = _cs.FindKernel("Warp");

            WarpConfig config = new WarpConfig();
            config.mapPath = "D:\\tmp\\hemi_pairs_596.json";
            config.srcPath = "D:\\tmp\\wallFish.jpg";
            config.outPath = "D:\\tmp\\wallGPUUnityOut.png";
            Mapping map = Mapping.Read(config.mapPath);
            tps = TPS2D.FitInverseNormalized(map, config.lambda);

            using Bitmap src = new Bitmap(config.srcPath);
            sourceTexture = new Texture2D(src.Width, src.Height);

            // copy the src image to sourceTexture

            MemoryStream ms = new MemoryStream();
            src.Save(ms, ImageFormat.Png);
            var buffer = new byte[ms.Length];
            ms.Position = 0;
            ms.Read(buffer, 0, buffer.Length);
            sourceTexture.LoadImage(buffer);
            sourceTexture.wrapMode = TextureWrapMode.Clamp;
            sourceTexture.filterMode = FilterMode.Bilinear;

            if (outputRT == null) {
                outputRT = new RenderTexture(sourceTexture.width, sourceTexture.height, 0, RenderTextureFormat.ARGB32);
                outputRT.enableRandomWrite = true;
                outputRT.Create();
            }

            warper = new GpuWarper(tpsWarpCS);
            bool clampEdge = /* your flag */ true;

            warper.RenderWholeGPU(sourceTexture, outputRT, tps, clampEdge);

            // If you still need a Bitmap on disk:
            Texture2D readback = new Texture2D(outputRT.width, outputRT.height, TextureFormat.RGBA32, false);
            RenderTexture.active = outputRT;
            readback.ReadPixels(new Rect(0, 0, outputRT.width, outputRT.height), 0, 0);
            readback.Apply();
            File.WriteAllBytes(config.outPath, readback.EncodeToPNG());
            Debug.Log("ALL DONE WITH GPU VERSION");
        }

        void OnDestroy () { warper?.Dispose(); }
    }

    public class GpuWarper : System.IDisposable {
        ComputeShader _cs;
        int _kernel;

        ComputeBuffer _txX, _txY, _txW, _tyX, _tyY, _tyW;
        bool _inited;

        public GpuWarper (ComputeShader cs) {
            _cs = cs;
            _kernel = _cs.FindKernel("Warp");
        }

        void EnsureBuffers ((float[] x, float[] y, float[] w) tx, (float[] x, float[] y, float[] w) ty) {
            ReleaseBuffers();
            _txX = new ComputeBuffer(tx.x.Length, sizeof(float));
            _txY = new ComputeBuffer(tx.y.Length, sizeof(float));
            _txW = new ComputeBuffer(tx.w.Length, sizeof(float));
            _tyX = new ComputeBuffer(ty.x.Length, sizeof(float));
            _tyY = new ComputeBuffer(ty.y.Length, sizeof(float));
            _tyW = new ComputeBuffer(ty.w.Length, sizeof(float));

            _txX.SetData(tx.x); _txY.SetData(tx.y); _txW.SetData(tx.w);
            _tyX.SetData(ty.x); _tyY.SetData(ty.y); _tyW.SetData(ty.w);
            _inited = true;
        }

        public void RenderWholeGPU (Texture src, RenderTexture dst, FisheyeFriend.TPS2D tps, bool clampEdge) {
            int sw = src.width, sh = src.height;
            int dw = dst.width, dh = dst.height;

            var tx = tps.ExportTxArrays();
            var ty = tps.ExportTyArrays();
            if (!_inited) EnsureBuffers(tx, ty);

            // Coeffs
            var txSpline = (FisheyeFriend.ThinPlateSpline)typeof(FisheyeFriend.TPS2D)
                .GetField("tx", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .GetValue(tps);
            var tySpline = (FisheyeFriend.ThinPlateSpline)typeof(FisheyeFriend.TPS2D)
                .GetField("ty", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .GetValue(tps);

            // Set shader resources
            _cs.SetTexture(_kernel, "_Src", src);
            _cs.SetTexture(_kernel, "_Dst", dst);

            _cs.SetBuffer(_kernel, "_TxX", _txX);
            _cs.SetBuffer(_kernel, "_TxY", _txY);
            _cs.SetBuffer(_kernel, "_TxW", _txW);
            _cs.SetBuffer(_kernel, "_TyX", _tyX);
            _cs.SetBuffer(_kernel, "_TyY", _tyY);
            _cs.SetBuffer(_kernel, "_TyW", _tyW);

            _cs.SetFloat("a0x", (float)txSpline.A0);
            _cs.SetFloat("axx", (float)txSpline.AX);
            _cs.SetFloat("ayx", (float)txSpline.AY);
            _cs.SetFloat("a0y", (float)tySpline.A0);
            _cs.SetFloat("axy", (float)tySpline.AX);
            _cs.SetFloat("ayy", (float)tySpline.AY);
            _cs.SetVector("dims", new Vector4(sw, sh, dw, dh));
            _cs.SetInt("nCtrl", tx.x.Length);
            _cs.SetInt("clampEdge", clampEdge ? 1 : 0);

            // Ensure RW format
            if (!dst.IsCreated()) {
                dst.enableRandomWrite = true;
                dst.Create();
            }

            // Dispatch (8x8 thread group)
            int gx = (dw + 7) / 8;
            int gy = (dh + 7) / 8;
            _cs.Dispatch(_kernel, gx, gy, 1);
        }

        public void Dispose () => ReleaseBuffers();
        void ReleaseBuffers () {
            _txX?.Dispose(); _txY?.Dispose(); _txW?.Dispose();
            _tyX?.Dispose(); _tyY?.Dispose(); _tyW?.Dispose();
            _txX = _txY = _txW = _tyX = _tyY = _tyW = null;
            _inited = false;
        }
    }
}