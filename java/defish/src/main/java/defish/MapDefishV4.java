package defish;


//HemiPairsDefishV4.java
//Constrained, symmetric data-driven defish (fitted to your hemi_pairs_384.json).
//- Horizontal x’ ~ u * S(|u|, r^2) + tiny even-in-y row shift (keeps verticals visually straight)
//- Vertical   y’ ~ v * G(|u|, v^2, r^2) + tiny symmetric lift
//- Normalization by half-height preserves the overall aspect ratio
//- EdgeGuard compresses radii close to the source rim to avoid ceiling/floor smearing

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;

public class MapDefishV4 {

 // ====== FITTED COEFFICIENTS (from hemi_pairs_384.json) ======
 // x' = u * (aX · phiX) + (hX0 + hX2 v^2 + hX4 v^4)
 //    where phiX = [1, r2, r4, |u|, |u|^2, |u|^3, v2, v4, r2|u|]
 private static final double[] AX = {
     0.82379661,  0.18767394,  0.03028218,  0.02865178, -0.13106995,
     0.01878403,  0.00637488,  0.06073611, -0.06097809
 };
 private static final double[] HX = {  // tiny row shift only
     -0.00024687, -0.00260790, 0.00463361
 };

 // y' = v * ( 1 + CY · phiYrest ) + (DY · psiY)
 // phiYrest = [r2, r4, |u|, |u|^2, |u|^3, |u|^4, v2, v4, r2|u|, |u|^2 v2, |u| v2, |u|^2 v2?]  (see buildPhiY)
 // psiY     = [1, v2, v4, |u|v2, |u|^2 v2]
 // (Center scale is anchored at 1 by construction.)
 private static final double[] CY = {
     -0.12950130,  0.19146071, -0.45776540,  0.68479890,
     -0.32960349, -0.12574818, -0.16286004,  0.30993409,
     -0.23996736,  0.03535305, -0.61623078,  0.36495654
 };
 private static final double[] DY = {
      0.01698363, -0.00850356, -0.00895734, -0.01032530,  0.01070970
 };

 // set true if you want strict top/bottom mirroring (zero the tiny global bias DY[0])
 private static final boolean STRICT_TOP_BOTTOM = false;

 // ====== Edge guard (prevents top/bottom smear). Tweak if needed. ======
 private static final double EDGE_START  = 0.96;  // start compressing radii beyond this (normalized)
 private static final double EDGE_ALPHA  = 0.35;  // how hard we compress (0..1)
 private static final double EDGE_POWER  = 2.0;   // how fast it ramps up

 public static void main(String[] args) throws Exception {
     if (args.length < 2) {
         System.out.println("Usage: java HemiPairsDefishV4 <input.jpg> <output.png> [outW outH]");
         System.out.println("Tip: use 1769 1197 to match your Hemi canvas.");
         return;
     }
     String inPath  = args[0];
     String outPath = args[1];
     BufferedImage src = ImageIO.read(new File(inPath));
     if (src == null) throw new RuntimeException("Cannot read " + inPath);

     int outW = (args.length >= 4) ? Integer.parseInt(args[2]) : src.getWidth();
     int outH = (args.length >= 4) ? Integer.parseInt(args[3]) : src.getHeight();

     BufferedImage dst = remap(src, outW, outH);
     ImageIO.write(dst, outPath.toLowerCase().endsWith(".jpg") ? "jpg" : "png", new File(outPath));
     System.out.println("Wrote " + outPath + " (" + outW + "x" + outH + ")");
 }

 public static BufferedImage remap(BufferedImage src, int outW, int outH) {
     final int W = src.getWidth(), H = src.getHeight();
     final double cx = W * 0.5, cy = H * 0.5;
     final double R  = H * 0.5;     // normalize with half-height (keeps aspect)

     final double ocx = outW * 0.5, ocy = outH * 0.5;
     final double Rout = outH * 0.5;

     BufferedImage dst = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     final int NEWTON_STEPS = 2;
     final double EPS = 1e-3;
     int[] line = new int[outW];

     for (int y = 0; y < outH; y++) {
         double yn = (y - ocy) / Rout;
         for (int x = 0; x < outW; x++) {
             double xn = (x - ocx) / Rout;

             // initial guess: identity
             double u = xn, v = yn;

             // invert forward(u,v) wavyEquals (xn,yn)
             for (int it = 0; it < NEWTON_STEPS; it++) {
                 double[] f0 = forward(u, v);
                 double ex = f0[0] - xn, ey = f0[1] - yn;

                 // Jacobian by small finite differences
                 double[] fx = forward(u + EPS, v);
                 double[] fy = forward(u, v + EPS);
                 double Jxx = (fx[0] - f0[0]) / EPS, Jxy = (fy[0] - f0[0]) / EPS;
                 double Jyx = (fx[1] - f0[1]) / EPS, Jyy = (fy[1] - f0[1]) / EPS;
                 double det = Jxx * Jyy - Jxy * Jyx;
                 if (Math.abs(det) < 1e-7) break;

                 double du = (-ex * Jyy + ey * Jxy) / det;
                 double dv = (-ey * Jxx + ex * Jyx) / det;
                 u += du; v += dv;
             }

             // Edge guard in (u,v) domain to avoid sampling beyond the fisheye circle
             double r = Math.hypot(u, v);
             if (r > EDGE_START) {
                 double t = (r - EDGE_START) / (1.0 - EDGE_START);
                 double compress = 1.0 - EDGE_ALPHA * Math.pow(t, EDGE_POWER);
                 if (compress < 0.0) compress = 0.0;
                 double rNew = EDGE_START + (1.0 - EDGE_START) * compress;
                 double s = (r > 0.0) ? (rNew / r) : 1.0;
                 u *= s; v *= s;
             }

             double sx = cx + u * R;
             double sy = cy + v * R;
             line[x] = sampleBilinearClamp(src, sx, sy);
         }
         dst.setRGB(0, y, outW, 1, line, 0, outW);
     }
     return dst;
 }

 // Forward mapping from normalized (u,v) -> (x', y')
 private static double[] forward(double u, double v) {
     double ax = Math.abs(u);
     double v2 = v*v;
     double r2 = u*u + v*v;

     // ---- Horizontal: x' ----
     double[] phiX = new double[]{
         1.0, r2, r2*r2, ax, ax*ax, ax*ax*ax, v2, v2*v2, r2*ax
     };
     double scaleX = dot(AX, phiX);
     double xp = u * scaleX + (HX[0] + HX[1]*v2 + HX[2]*v2*v2);

     // ---- Vertical: y' (center scale anchored at 1) ----
     double[] phiY = buildPhiY(ax, v2, r2);
     double[] psiY = new double[]{
         STRICT_TOP_BOTTOM ? 0.0 : 1.0,  // optional tiny global bias
         v2, v2*v2, ax*v2, (ax*ax)*v2
     };
     double gy = 1.0 + dot(CY, phiY);      // multiplicative part
     double yp = v * gy + dot(DY, psiY);   // plus small symmetric lift

     return new double[]{xp, yp};
 }

 private static double[] buildPhiY(double ax, double v2, double r2){
     // Rest of vertical basis (no constant term)
     return new double[]{
         r2, r2*r2, ax, ax*ax, ax*ax*ax, Math.pow(ax,4),
         v2, v2*v2, r2*ax, (ax*ax)*v2, ax*v2, (ax*ax)*v2
     };
 }

 private static double dot(double[] a, double[] b){
     double s = 0.0; for (int i=0;i<a.length;i++) s += a[i]*b[i]; return s;
 }

 // Bilinear with clamp
 private static int sampleBilinearClamp(BufferedImage img, double x, double y) {
     int w = img.getWidth(), h = img.getHeight();
     if (x < 0) x = 0; if (y < 0) y = 0;
     if (x > w-1) x = w-1; if (y > h-1) y = h-1;
     int x0 = (int)Math.floor(x), y0 = (int)Math.floor(y);
     int x1 = Math.min(x0+1, w-1), y1 = Math.min(y0+1, h-1);
     double fx = x - x0, fy = y - y0;

     int c00 = img.getRGB(x0,y0), c10 = img.getRGB(x1,y0), c01 = img.getRGB(x0,y1), c11 = img.getRGB(x1,y1);
     int a = lerp2((c00>>>24)&255, (c10>>>24)&255, (c01>>>24)&255, (c11>>>24)&255, fx, fy);
     int r = lerp2((c00>>>16)&255, (c10>>>16)&255, (c01>>>16)&255, (c11>>>16)&255, fx, fy);
     int g = lerp2((c00>>>8)&255,  (c10>>>8)&255,  (c01>>>8)&255,  (c11>>>8)&255,  fx, fy);
     int b = lerp2( c00&255,        c10&255,        c01&255,        c11&255,       fx, fy);
     return (a<<24)|(r<<16)|(g<<8)|b;
 }
 private static int lerp2(int c00,int c10,int c01,int c11,double fx,double fy){
     double i0 = c00 + (c10 - c00)*fx;
     double i1 = c01 + (c11 - c01)*fx;
     return (int)Math.round(i0 + (i1 - i0)*fy);
 }
}
