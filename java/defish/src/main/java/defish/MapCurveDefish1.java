package defish;

//FisheyeHemiWarp.java
//Apply a Hemi-like warp using precomputed curves CSV: angle_deg, r_src, r_dst
//- Inverse mapping: for each output pixel (angle theta, radius r_dst) -> r_src via curve inversion
//- Angle interpolation (wrap-aware) + radial interpolation for smoothness
//- Radii normalized by half-height (H/2) about the image center
//- High quality Catmull–Rom bicubic (default) or bilinear sampling
//
//Usage:
//javac FisheyeHemiWarp.java
//java FisheyeHemiWarp <curves.csv> <input.jpg> <output.png> [outW outH] [--bilinear] [--clampEdge] [--center=a,b]
//
//CSV format (headers required): angle_deg,r_src,r_dst
//angle_deg: curve angle in degrees (can be any set; not required to be uniform)
//r_src:     source normalized radius (half-height normalization)
//r_dst:     destination normalized radius (half-height normalization)
//
//Notes:
//- Default output size equals input size; use [outW outH] to set a canvas (e.g., 1769 1197).
//- “center” defaults to image center; override with --center=0.5,0.5 (fractions of width,height).

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Pattern;

public class MapCurveDefish1 {

 public static void main(String[] args) throws Exception {
     if (args.length < 3) {
         System.out.println("Usage: java FisheyeHemiWarp <curves.csv> <input.jpg> <output.png> [outW outH] [--bilinear] [--clampEdge] [--center=a,b]");
         return;
     }
     String curvesPath = args[0];
     String inPath     = args[1];
     String outPath    = args[2];

     Integer outWArg = null, outHArg = null;
     boolean bilinear = false;
     boolean clampEdge = false;
     double centerFracX = 0.5, centerFracY = 0.5;

     for (int i = 3; i < args.length; i++) {
         String s = args[i];
         if (s.equals("--bilinear")) bilinear = true;
         else if (s.equals("--clampEdge")) clampEdge = true;
         else if (s.startsWith("--center=")) {
             String t = s.substring("--center=".length());
             String[] parts = t.split(",");
             if (parts.length == 2) {
                 centerFracX = Double.parseDouble(parts[0]);
                 centerFracY = Double.parseDouble(parts[1]);
             }
         } else if (outWArg == null && s.matches("\\d+")) {
             outWArg = Integer.parseInt(s);
         } else if (outHArg == null && s.matches("\\d+")) {
             outHArg = Integer.parseInt(s);
         }
     }

     // Load curves
     CurveSet curves = CurveSet.load(new File(curvesPath));
     System.out.println("Loaded " + curves.anglesDeg.length + " curves; samples vary per angle.");

     // Read input and convert to ARGB
     BufferedImage srcIn = ImageIO.read(new File(inPath));
     if (srcIn == null) throw new IOException("Cannot read input image: " + inPath);
     BufferedImage src = toIntARGB(srcIn);
     final int sw = src.getWidth(), sh = src.getHeight();

     // Output size
     final int outW = (outWArg != null) ? outWArg : sw;
     final int outH = (outHArg != null) ? outHArg : sh;

     // Centers and radii (half-height normalization)
     final double scx = sw * 0.5, scy = sh * 0.5, sR = sh * 0.5;
     final double ocx = outW * centerFracX, ocy = outH * centerFracY, oR = outH * 0.5;

     // Prepare output
     BufferedImage out = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     int[] outData = ((DataBufferInt) out.getRaster().getDataBuffer()).getData();
     int[] srcData = ((DataBufferInt) src.getRaster().getDataBuffer()).getData();

     final Sampler sampler = bilinear ? new BilinearSampler() : new BicubicSampler();
     final int sStride = sw;

     // Render
     for (int y = 0; y < outH; y++) {
         double vy = (y - ocy) / oR;
         int row = y * outW;
         for (int x = 0; x < outW; x++) {
             double vx = (x - ocx) / oR;

             // Output polar (normalized)
             double phi = Math.atan2(vy, vx);  // [-pi, pi)
             double rDst = Math.hypot(vx, vy);

             // Invert curve(s): rDst -> rSrc at angle phi (with angle interpolation)
             double rSrc = curves.invertRadius(rDst, phi);

             // Map to source pixels
             double sx = scx + rSrc * sR * Math.cos(phi);
             double sy = scy + rSrc * sR * Math.sin(phi);

             int argb;
             if (!clampEdge && (sx < -0.5 || sy < -0.5 || sx > sw - 0.5 || sy > sh - 0.5)) {
                 argb = 0x00000000; // transparent outside
             } else {
                 argb = sampler.sample(srcData, sw, sh, sStride, sx, sy, clampEdge);
             }
             outData[row + x] = argb;
         }
     }

     ImageIO.write(out, outPath.toLowerCase().endsWith(".jpg") ? "jpg" : "png", new File(outPath));
     System.out.println("Wrote " + outPath + "  (" + outW + "x" + outH + "), " + (bilinear ? "bilinear" : "bicubic")
             + (clampEdge ? ", clampEdge" : ", alpha outside"));
 }

 // ---- Curves container & inversion ----
 static class Curve {
     final double angleDeg;          // angle of this curve
     final double[] rSrc;            // source radii (monotone increasing)
     final double[] rDst;            // destination radii (monotone increasing)

     Curve(double angleDeg, double[] rSrc, double[] rDst) {
         this.angleDeg = angleDeg;
         this.rSrc = rSrc;
         this.rDst = rDst;
     }

     // Invert monotone rDst -> rSrc (piecewise linear; clamps at ends)
     double invert(double rdst) {
         // Clamp within sampled range
         if (rdst <= rDst[0]) {
             if (rDst.length >= 2) {
                 int i = firstIncreasingIndex();
                 if (i > 0) i = 0;
                 double t = safeT(rdst, rDst[0], rDst[1]);
                 return rSrc[0] + t * (rSrc[1] - rSrc[0]);
             }
             return rSrc[0];
         }
         int n = rDst.length;
         if (rdst >= rDst[n-1]) {
             if (n >= 2) {
                 double t = safeT(rdst, rDst[n-2], rDst[n-1]);
                 return rSrc[n-2] + t * (rSrc[n-1] - rSrc[n-2]);
             }
             return rSrc[n-1];
         }
         // Binary search
         int lo = 0, hi = n - 1;
         while (hi - lo > 1) {
             int mid = (lo + hi) >>> 1;
             if (rDst[mid] <= rdst) lo = mid; else hi = mid;
         }
         double denom = (rDst[hi] - rDst[lo]);
         double t = denom == 0 ? 0 : (rdst - rDst[lo]) / denom;
         return rSrc[lo] + t * (rSrc[hi] - rSrc[lo]);
     }

     // Find first index where rDst[i+1] > rDst[i], to avoid 0-length segment
     private int firstIncreasingIndex() {
         for (int i = 0; i < rDst.length - 1; i++) {
             if (rDst[i+1] > rDst[i]) return i;
         }
         return 0;
     }

     private static double safeT(double x, double a, double b) {
         double d = b - a;
         if (Math.abs(d) < 1e-12) return 0.0;
         return (x - a) / d;
     }
 }

 static class CurveSet {
     final Curve[] curves;
     final double[] anglesDeg; // sorted

     CurveSet(Curve[] curves) {
         this.curves = curves;
         Arrays.sort(this.curves, Comparator.comparingDouble(c -> c.angleDeg));
         this.anglesDeg = new double[this.curves.length];
         for (int i = 0; i < curves.length; i++) anglesDeg[i] = this.curves[i].angleDeg;
     }

     static CurveSet load(File csv) throws IOException {
         ArrayList<String> lines = new ArrayList<>();
         try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(csv), StandardCharsets.UTF_8))) {
             String L; while ((L = br.readLine()) != null) lines.add(L.trim());
         }
         if (lines.isEmpty()) throw new IOException("Empty curves file.");

         // Header detection
         String header = lines.get(0);
         int idxAngle = -1, idxRSrc = -1, idxRDst = -1;
         String[] head = header.split("\\s*,\\s*");
         for (int i = 0; i < head.length; i++) {
             String h = head[i].toLowerCase(Locale.US);
             if (h.contains("angle")) idxAngle = i;
             else if (h.contains("r_src") || h.equals("rsrc")) idxRSrc = i;
             else if (h.contains("r_dst") || h.equals("rdst")) idxRDst = i;
         }
         if (idxAngle < 0 || idxRSrc < 0 || idxRDst < 0)
             throw new IOException("CSV header must include angle_deg,r_src,r_dst (found: " + header + ")");

         // Group by angle
         TreeMap<Double, ArrayList<double[]>> groups = new TreeMap<>();
         Pattern num = Pattern.compile("[+-]?(?:\\d+\\.?\\d*|\\d*\\.\\d+)(?:[eE][+-]?\\d+)?");

         for (int li = 1; li < lines.size(); li++) {
             String s = lines.get(li);
             if (s.isEmpty()) continue;
             String[] cols = s.split("\\s*,\\s*");
             if (cols.length < Math.max(idxRDst, Math.max(idxAngle, idxRSrc)) + 1) continue;

             try {
                 double ang = Double.parseDouble(cols[idxAngle]);
                 double rs = Double.parseDouble(cols[idxRSrc]);
                 double rd = Double.parseDouble(cols[idxRDst]);
                 groups.computeIfAbsent(ang, k -> new ArrayList<>()).add(new double[]{rs, rd});
             } catch (NumberFormatException ex) {
                 // skip malformed line
             }
         }

         ArrayList<Curve> list = new ArrayList<>();
         for (Map.Entry<Double, ArrayList<double[]>> e : groups.entrySet()) {
             ArrayList<double[]> pts = e.getValue();
             pts.sort(Comparator.comparingDouble(a -> a[0])); // sort by r_src
             // Build monotone arrays in r_dst order for robust inversion
             // Ensure both r_src and r_dst are non-decreasing with r_src; then sort by r_dst for invert()
             ArrayList<double[]> filt = new ArrayList<>();
             double lastRS = -Double.MAX_VALUE, lastRD = -Double.MAX_VALUE;
             for (double[] p : pts) {
                 double rs = p[0], rd = p[1];
                 if (rs < lastRS) continue;
                 if (rd < lastRD) rd = lastRD; // enforce monotone non-decreasing rd
                 filt.add(new double[]{rs, rd});
                 lastRS = rs; lastRD = rd;
             }
             // Sort by r_dst for inversion
             filt.sort(Comparator.comparingDouble(a -> a[1]));
             double[] rsrc = new double[filt.size()];
             double[] rdst = new double[filt.size()];
             for (int i = 0; i < filt.size(); i++) {
                 rsrc[i] = filt.get(i)[0];
                 rdst[i] = filt.get(i)[1];
             }
             if (rsrc.length >= 2) {
                 list.add(new Curve(e.getKey(), rsrc, rdst));
             }
         }
         if (list.isEmpty()) throw new IOException("No usable curves found.");
         return new CurveSet(list.toArray(new Curve[0]));
     }

     // Invert with angle interpolation (wrap-aware).
     double invertRadius(double rDst, double phiRad) {
         double angDeg = Math.toDegrees(phiRad);
         // Map to [0,360)
         angDeg = ((angDeg % 360.0) + 360.0) % 360.0;

         // Find bracketing curves in angle, with wrap-around
         int n = anglesDeg.length;
         // binary search insert point
         int idx = Arrays.binarySearch(anglesDeg, angDeg);
         if (idx >= 0) {
             // exact angle curve exists
             return curves[idx].invert(rDst);
         } else {
             int ip = -idx - 1;               // first angle greater than angDeg
             int i1 = (ip % n + n) % n;       // wrap
             int i0 = (i1 - 1 + n) % n;

             double a0 = anglesDeg[i0];
             double a1 = anglesDeg[i1];

             // Handle wrap gap (e.g., 359° to 0°)
             double span = a1 - a0;
             if (span <= 0) span += 360.0;

             double t;
             if (span <= 0) {
                 t = 0.0;
             } else {
                 double d = angDeg - a0;
                 if (d < 0) d += 360.0; // wrap distance from a0 up to angDeg
                 t = d / span;
             }

             double r0 = curves[i0].invert(rDst);
             double r1 = curves[i1].invert(rDst);
             return r0 + t * (r1 - r0);
         }
     }
 }

 // ---- Sampling helpers ----

 interface Sampler {
     int sample(int[] data, int w, int h, int stride, double x, double y, boolean clampEdge);
 }

 static class BilinearSampler implements Sampler {
     @Override
     public int sample(int[] data, int w, int h, int stride, double x, double y, boolean clamp) {
         if (clamp) {
             if (x < 0) x = 0; if (y < 0) y = 0;
             if (x > w-1) x = w-1; if (y > h-1) y = h-1;
         }
         int x0 = (int)Math.floor(x), y0 = (int)Math.floor(y);
         int x1 = x0 + 1, y1 = y0 + 1;
         double fx = x - x0, fy = y - y0;

         int c00 = getPixel(data, w, h, stride, x0, y0, clamp);
         int c10 = getPixel(data, w, h, stride, x1, y0, clamp);
         int c01 = getPixel(data, w, h, stride, x0, y1, clamp);
         int c11 = getPixel(data, w, h, stride, x1, y1, clamp);

         int a = lerp2((c00>>>24)&255, (c10>>>24)&255, (c01>>>24)&255, (c11>>>24)&255, fx, fy);
         int r = lerp2((c00>>>16)&255, (c10>>>16)&255, (c01>>>16)&255, (c11>>>16)&255, fx, fy);
         int g = lerp2((c00>>>8)&255,  (c10>>>8)&255,  (c01>>>8)&255,  (c11>>>8)&255,  fx, fy);
         int b = lerp2( c00&255,        c10&255,        c01&255,        c11&255,       fx, fy);
         return (a<<24)|(r<<16)|(g<<8)|b;
     }
 }

 static class BicubicSampler implements Sampler {
     @Override
     public int sample(int[] data, int w, int h, int stride, double x, double y, boolean clamp) {
         int x1 = (int)Math.floor(x);
         int y1 = (int)Math.floor(y);
         double tx = x - x1;
         double ty = y - y1;

         int[] xs = new int[4];
         int[] ys = new int[4];
         for (int i = -1; i <= 2; i++) { xs[i+1] = x1 + i; ys[i+1] = y1 + i; }
         double[] wx = catmullRomWeights(tx);
         double[] wy = catmullRomWeights(ty);

         double a=0,r=0,g=0,b=0, Wsum=0;
         for (int j=0;j<4;j++){
             int yy = ys[j];
             for (int i=0;i<4;i++){
                 int xx = xs[i];
                 int c = getPixel(data, w, h, stride, xx, yy, clamp);
                 if (c == 0 && !clamp) continue; // transparent contributes nothing
                 double wgt = wx[i]*wy[j];
                 a += wgt * ((c>>>24)&255);
                 r += wgt * ((c>>>16)&255);
                 g += wgt * ((c>>>8) &255);
                 b += wgt * ( c      &255);
                 Wsum += wgt;
             }
         }
         if (Wsum <= 1e-12) return 0;
         int ia = clamp255((int)Math.round(a/Wsum));
         int ir = clamp255((int)Math.round(r/Wsum));
         int ig = clamp255((int)Math.round(g/Wsum));
         int ib = clamp255((int)Math.round(b/Wsum));
         return (ia<<24)|(ir<<16)|(ig<<8)|ib;
     }
 }

 static int getPixel(int[] data, int w, int h, int stride, int x, int y, boolean clamp) {
     if (clamp) {
         if (x < 0) x = 0; if (y < 0) y = 0;
         if (x > w-1) x = w-1; if (y > h-1) y = h-1;
     } else {
         if (x < 0 || y < 0 || x >= w || y >= h) return 0x00000000;
     }
     return data[y*stride + x];
 }

 static int lerp2(int c00,int c10,int c01,int c11,double fx,double fy){
     double i0 = c00 + (c10 - c00)*fx;
     double i1 = c01 + (c11 - c01)*fx;
     return clamp255((int)Math.round(i0 + (i1 - i0)*fy));
 }

 static int clamp255(int v){ return (v<0)?0:((v>255)?255:v); }

 static double[] catmullRomWeights(double t) {
     double a = -0.5;
     double t2 = t*t, t3 = t2*t;
     double w0 = a*(-t3 + 2*t2 - t);
     double w1 = (a+2)*t3 + (-a-3)*t2 + 1;
     double w2 = (a+2)*(-t3) + (2*a+3)*t2 + (-a)*t;
     double w3 = a*(t3 - t2);
     return new double[]{w0,w1,w2,w3};
 }

 static BufferedImage toIntARGB(BufferedImage in) {
     if (in.getType() == BufferedImage.TYPE_INT_ARGB) return in;
     BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
     Graphics2D g = out.createGraphics();
     g.setComposite(AlphaComposite.Src);
     g.drawImage(in, 0, 0, null);
     g.dispose();
     return out;
 }
}
