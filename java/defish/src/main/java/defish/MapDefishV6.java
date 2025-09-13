package defish;


//MapDefishV6.java
//From-scratch, data-driven defish using Thin-Plate Splines (inverse warp: target->source).
//- Converts input to TYPE_INT_ARGB (fixes DataBufferByte/DataBufferInt cast issue).
//- Fits two TPS maps X(dst)->src_u and Y(dst)->src_v in normalized coords (centered, half-height scaling).
//- Mirrors control points left/right (and top/bottom by default) to enforce Hemi-like symmetry.
//- Catmull–Rom bicubic sampling; outside-source samples become transparent unless --clampEdge.
//
//Build:
//javac MapDefishV6.java
//Run examples:
//java MapDefishV6 cloister_before.jpg hemi_pairs_384.json hemi_like.png
//java MapDefishV6 cloister_before.jpg hemi_pairs_384.json hemi_like_1769x1197.png 1769 1197
//Options:
//--lambda=1e-5   (TPS smoothing; lower=tighter, higher=smoother)
//--noMirrorY     (don’t mirror top/bottom; only left/right)
//--clampEdge     (clamp sampling outside src instead of transparent)

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MapDefishV6 {

 // ---------- CLI ----------
 public static void main(String[] args) throws Exception {
     if (args.length < 3) {
         System.out.println("Usage: java MapDefishV6 <input.jpg> <mapping.json> <output.png> [outW outH] [--lambda=1e-5] [--noMirrorY] [--clampEdge]");
         return;
     }
     String inPath   = args[0];
     String mapPath  = args[1];
     String outPath  = args[2];

     Integer outWArg = null, outHArg = null;
     double lambda = 1e-5;
     boolean mirrorY = true;
     boolean clampEdge = false;

     for (int i = 3; i < args.length; i++) {
         String s = args[i];
         if (s.startsWith("--lambda=")) {
             lambda = Double.parseDouble(s.substring("--lambda=".length()));
         } else if (s.equals("--noMirrorY")) {
             mirrorY = false;
         } else if (s.equals("--clampEdge")) {
             clampEdge = true;
         } else if (outWArg == null && s.matches("\\d+")) {
             outWArg = Integer.parseInt(s);
         } else if (outHArg == null && s.matches("\\d+")) {
             outHArg = Integer.parseInt(s);
         }
     }

     BufferedImage srcIn = ImageIO.read(new File(inPath));
     if (srcIn == null) throw new IOException("Cannot read input image: " + inPath);

     // Convert to TYPE_INT_ARGB to ensure an int[] backing buffer (fixes your ClassCastException).
     BufferedImage src = toIntARGB(srcIn);

     Mapping m = Mapping.read(new File(mapPath));
     // Default output size: mapping target canvas
     int outW = (outWArg != null) ? outWArg : m.tW;
     int outH = (outHArg != null) ? outHArg : m.tH;

     // Build training data (inverse: target->source), normalized by half-height
     TrainingData td = buildInverseTraining(m);

     // Symmetry augmentation: mirror L/R always, T/B optionally
     td = augmentSymmetry(td, true, mirrorY);

     // Fit TPS for u(x',y') and v(x',y')
     ThinPlateSpline tpsU = ThinPlateSpline.fit(td.xp, td.yp, td.u, lambda);
     ThinPlateSpline tpsV = ThinPlateSpline.fit(td.xp, td.yp, td.v, lambda);

     // Render with bicubic sampling
     BufferedImage out = renderBicubic(src, outW, outH, tpsU, tpsV, clampEdge);

     ImageIO.write(out, outPath.toLowerCase().endsWith(".jpg") ? "jpg" : "png", new File(outPath));
     System.out.println("Wrote " + outPath + "  (" + outW + "x" + outH + ")");
 }

 // Convert any BufferedImage to TYPE_INT_ARGB (so DataBufferInt is guaranteed).
 private static BufferedImage toIntARGB(BufferedImage in) {
     if (in.getType() == BufferedImage.TYPE_INT_ARGB) return in;
     BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
     Graphics2D g = out.createGraphics();
     g.setComposite(AlphaComposite.Src);
     g.drawImage(in, 0, 0, null);
     g.dispose();
     return out;
 }

 // ---------- Mapping parser & normalization ----------
 static class Mapping {
     int sW, sH, tW, tH;
     java.util.List<Pair> pairs = new ArrayList<>();
     static class Pair { double sx, sy, dx, dy; } // pixel coords

     static Mapping read(File f) throws IOException {
         String text = readAll(f);
         Mapping m = new Mapping();
         m.sW = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.sH = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");
         m.tW = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.tH = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");

         Pattern p = Pattern.compile(
             "\\{\\s*\"id\"\\s*:\\s*\\d+[^}]*?\"src\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)[^}]*?\\}[^}]*?\"dst\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)",
             Pattern.DOTALL);
         Matcher mtr = p.matcher(text);
         int count = 0;
         while (mtr.find()) {
             Pair pr = new Pair();
             pr.sx = Double.parseDouble(mtr.group(1));
             pr.sy = Double.parseDouble(mtr.group(2));
             pr.dx = Double.parseDouble(mtr.group(3));
             pr.dy = Double.parseDouble(mtr.group(4));
             m.pairs.add(pr);
             count++;
         }
         if (count == 0) throw new IOException("No pairs found in mapping JSON.");
         System.out.println("Loaded mapping: " + count + " pairs. Source " + m.sW + "x" + m.sH + "  Target " + m.tW + "x" + m.tH);
         return m;
     }
     static String readAll(File f) throws IOException {
         try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(f)))) {
             StringBuilder sb = new StringBuilder(1<<20);
             String line; while ((line = br.readLine()) != null) sb.append(line).append('\n');
             return sb.toString();
         }
     }
     static int readInt(String text, String regex) throws IOException {
         Matcher m = Pattern.compile(regex).matcher(text);
         if (!m.find()) throw new IOException("Missing field for regex: " + regex);
         return Integer.parseInt(m.group(1));
     }
 }

 static class TrainingData {
     double[] xp, yp;  // domain: target-normalized coords (x',y')
     double[] u, v;    // range:  source-normalized coords (u,v)
 }

 static TrainingData buildInverseTraining(Mapping m) {
     double tcx = m.tW * 0.5, tcy = m.tH * 0.5, tR = m.tH * 0.5;
     double scx = m.sW * 0.5, scy = m.sH * 0.5, sR = m.sH * 0.5;

     int N = m.pairs.size();
     TrainingData td = new TrainingData();
     td.xp = new double[N]; td.yp = new double[N];
     td.u  = new double[N]; td.v  = new double[N];

     int k = 0;
     for (Mapping.Pair p : m.pairs) {
         double xn = (p.dx - tcx) / tR;
         double yn = (p.dy - tcy) / tR;
         double un = (p.sx - scx) / sR;
         double vn = (p.sy - scy) / sR;
         td.xp[k] = xn; td.yp[k] = yn; td.u[k] = un; td.v[k] = vn; k++;
     }
     return td;
 }

 static TrainingData augmentSymmetry(TrainingData td, boolean mirrorX, boolean mirrorY) {
     ArrayList<double[]> rows = new ArrayList<>();
     for (int i = 0; i < td.xp.length; i++) {
         rows.add(new double[]{ td.xp[i], td.yp[i], td.u[i], td.v[i] });
         if (mirrorX) rows.add(new double[]{ -td.xp[i], td.yp[i], -td.u[i], td.v[i] });
         if (mirrorY) rows.add(new double[]{ td.xp[i], -td.yp[i], td.u[i], -td.v[i] });
         if (mirrorX && mirrorY) rows.add(new double[]{ -td.xp[i], -td.yp[i], -td.u[i], -td.v[i] });
     }
     TrainingData out = new TrainingData();
     int M = rows.size();
     out.xp = new double[M]; out.yp = new double[M]; out.u = new double[M]; out.v = new double[M];
     for (int i = 0; i < M; i++) {
         double[] r = rows.get(i);
         out.xp[i] = r[0]; out.yp[i] = r[1]; out.u[i] = r[2]; out.v[i] = r[3];
     }
     System.out.println("Symmetry-augmented control points: " + M + " (from " + td.xp.length + ")");
     return out;
 }

 // ---------- Thin-Plate Spline ----------
 static class ThinPlateSpline {
     final double[] x, y;  // control (domain)
     final double[] w;     // radial weights
     final double a0, ax, ay; // affine part

     private ThinPlateSpline(double[] x, double[] y, double[] w, double a0, double ax, double ay) {
         this.x = x; this.y = y; this.w = w; this.a0 = a0; this.ax = ax; this.ay = ay;
     }

     // Fit TPS: given control points (xi, yi) and target values zi
     static ThinPlateSpline fit(double[] xi, double[] yi, double[] zi, double lambda) {
         int n = xi.length;
         int dim = n + 3;
         double[][] A = new double[dim][dim];
         double[]   b = new double[dim];

         // K + lambda*I
         for (int i = 0; i < n; i++) {
             for (int j = i; j < n; j++) {
                 double r = Math.hypot(xi[i]-xi[j], yi[i]-yi[j]);
                 double v = U(r);
                 if (i == j) v += lambda;
                 A[i][j] = v;
                 A[j][i] = v;
             }
         }

         // P (affine) blocks
         for (int i = 0; i < n; i++) {
             A[i][n+0] = 1.0;       A[n+0][i] = 1.0;
             A[i][n+1] = xi[i];     A[n+1][i] = xi[i];
             A[i][n+2] = yi[i];     A[n+2][i] = yi[i];
         }
         // lower-right 3x3 is zeros by construction

         // RHS
         for (int i = 0; i < n; i++) b[i] = zi[i];

         // Solve A * coeff = b
         double[] coeff = solveSymmetric(A, b);

         double[] w = new double[n];
         System.arraycopy(coeff, 0, w, 0, n);
         double a0 = coeff[n+0], ax = coeff[n+1], ay = coeff[n+2];
         return new ThinPlateSpline(xi.clone(), yi.clone(), w, a0, ax, ay);
     }

     double eval(double px, double py) {
         double s = a0 + ax*px + ay*py;
         for (int i = 0; i < x.length; i++) {
             double r = Math.hypot(px - x[i], py - y[i]);
             s += w[i] * U(r);
         }
         return s;
     }

     static double U(double r) {
         if (r <= 1e-12) return 0.0;
         double r2 = r*r;
         return r2 * Math.log(r);
     }

     // Simple symmetric solver with partial pivoting (dense, small)
     static double[] solveSymmetric(double[][] A, double[] b) {
         int n = b.length;
         double[][] M = new double[n][n];
         for (int i = 0; i < n; i++) System.arraycopy(A[i], 0, M[i], 0, n);
         double[] x = new double[n];
         double[] rhs = b.clone();

         for (int k = 0; k < n; k++) {
             int piv = k;
             double max = Math.abs(M[k][k]);
             for (int i = k+1; i < n; i++) {
                 double v = Math.abs(M[i][k]);
                 if (v > max) { max = v; piv = i; }
             }
             if (max < 1e-14) throw new RuntimeException("Singular system in TPS solve at k=" + k);
             if (piv != k) {
                 double[] tmp = M[k]; M[k] = M[piv]; M[piv] = tmp;
                 double t = rhs[k]; rhs[k] = rhs[piv]; rhs[piv] = t;
             }
             double akk = M[k][k];
             for (int i = k+1; i < n; i++) {
                 double f = M[i][k] / akk;
                 rhs[i] -= f * rhs[k];
                 for (int j = k; j < n; j++) M[i][j] -= f * M[k][j];
             }
         }
         for (int i = n-1; i >= 0; i--) {
             double s = rhs[i];
             for (int j = i+1; j < n; j++) s -= M[i][j] * x[j];
             x[i] = s / M[i][i];
         }
         return x;
     }
 }

 // ---------- Rendering with bicubic sampling ----------
 static BufferedImage renderBicubic(BufferedImage src, int outW, int outH,
                                    ThinPlateSpline tpsU, ThinPlateSpline tpsV,
                                    boolean clampEdge) {

     final int W = src.getWidth(), H = src.getHeight();
     final double scx = W * 0.5, scy = H * 0.5, sR = H * 0.5;
     final double ocx = outW * 0.5, ocy = outH * 0.5, oR = outH * 0.5;

     BufferedImage out = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     int[] dst = ((DataBufferInt) out.getRaster().getDataBuffer()).getData();
     int[] srcData = ((DataBufferInt) src.getRaster().getDataBuffer()).getData();
     final int sw = W, sh = H, sstride = W;

     for (int y = 0; y < outH; y++) {
         double yp = (y - ocy) / oR;          // normalized target coords
         int row = y * outW;
         for (int x = 0; x < outW; x++) {
             double xp = (x - ocx) / oR;

             // inverse mapping: (x',y') -> (u,v)
             double u = tpsU.eval(xp, yp);
             double v = tpsV.eval(xp, yp);

             // back to source pixels
             double sx = scx + u * sR;
             double sy = scy + v * sR;

             int argb;
             if (!clampEdge && (sx < -0.5 || sy < -0.5 || sx > W-0.5 || sy > H-0.5)) {
                 argb = 0x00000000; // transparent outside
             } else {
                 argb = sampleBicubic(srcData, sw, sh, sstride, sx, sy, clampEdge);
             }
             dst[row + x] = argb;
         }
     }
     return out;
 }

 // Catmull–Rom bicubic sampling with optional clamp
 static int sampleBicubic(int[] data, int w, int h, int stride, double x, double y, boolean clamp) {
     int x1 = (int)Math.floor(x);
     int y1 = (int)Math.floor(y);
     double tx = x - x1;
     double ty = y - y1;

     int[] xs = new int[4];
     int[] ys = new int[4];
     for (int i = -1; i <= 2; i++) {
         xs[i+1] = x1 + i;
         ys[i+1] = y1 + i;
     }
     if (clamp) {
         for (int i=0;i<4;i++){ xs[i] = clamp(xs[i],0,w-1); ys[i]=clamp(ys[i],0,h-1); }
     }
     double[] wx = catmullRomWeights(tx);
     double[] wy = catmullRomWeights(ty);

     double a=0,r=0,g=0,b=0, Wsum=0;
     for (int j=0;j<4;j++){
         int yy = ys[j];
         if (yy < 0 || yy >= h) continue;
         for (int i=0;i<4;i++){
             int xx = xs[i];
             if (xx < 0 || xx >= w) continue;
             double wgt = wx[i]*wy[j];
             int c = data[yy*stride + xx];
             double ca=(c>>>24)&255, cr=(c>>>16)&255, cg=(c>>>8)&255, cb=c&255;
             a += wgt*ca; r+=wgt*cr; g+=wgt*cg; b+=wgt*cb; Wsum+=wgt;
         }
     }
     if (Wsum <= 1e-9) return 0x00000000;
     int ia = clamp((int)Math.round(a/Wsum),0,255);
     int ir = clamp((int)Math.round(r/Wsum),0,255);
     int ig = clamp((int)Math.round(g/Wsum),0,255);
     int ib = clamp((int)Math.round(b/Wsum),0,255);
     return (ia<<24)|(ir<<16)|(ig<<8)|ib;
 }

 static int clamp(int v,int lo,int hi){ return (v<lo)?lo:((v>hi)?hi:v); }

 static double[] catmullRomWeights(double t) {
     // centered Catmull–Rom (a = -0.5)
     double a = -0.5;
     double t2 = t*t, t3 = t2*t;
     double w0 = a*(-t3 + 2*t2 - t);
     double w1 = (a+2)*t3 + (-a-3)*t2 + 1;
     double w2 = (a+2)*(-t3) + (2*a+3)*t2 + (-a)*t;
     double w3 = a*(t3 - t2);
     return new double[]{w0,w1,w2,w3};
 }
}
