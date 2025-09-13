package defish;

//HemiPairsDefishV3.java
//Constrained, symmetric mapping fitted to hemi_pairs_384.json
//- Horizontal depends primarily on u (x in fisheye) to keep verticals straight.
//- Small, even-in-y row term is kept but heavily damped (<= ~2 px over full height).
//- Normalization by half-height preserves aspect ratio.

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;

public class MapDefishV3 {

 // ---- BASIS ORDER (keep if you later refit) ----
 // phiX (12): [1, r2, r2^2, |u|, |u|^2, |u|^3, |u|^4, v2, v4, |u|v2, |u|^2 v2, r2|u|]
 // xp = u * dot(ax, phiX) + (hx0 + hx2 v2 + hx4 v4)
 // psiY (10): [1, r2, r2^2, |u|, |u|^2, |u|^3, v2, v4, |u|v2, |u|^2 v2]
 // yp = v * dot(ay, psiY) + (hy0 + hy2 v2 + hy4 v4 + hy_uv |u| v2 + hy_u2v |u|^2 v2)

 // ---- NEW COEFFICIENTS (fitted with corner weighting + ridge) ----
 private static final double[] aX = {
     0.80687108,  0.07778051,  0.03261357,  0.13065566, -0.22614705,
     0.14048651, -0.04932147,  0.11547168,  0.05925338, -0.03968910,
    -0.00920978, -0.01627468
 };
 private static final double[] hX = {
     0.00025784, -0.00512527, 0.00773303   // hx0, hx2, hx4  (tiny, keeps horizon tidy)
 };

 private static final double[] aY = {
     0.84189622,  0.02749204, -0.04290115, -0.00519630, -0.07521969,
     0.03005192, -0.01025863,  0.01647407, -0.01256360,  0.00395744
 };
 private static final double[] hY = {
     0.01736163, -0.01214193, -0.00017223, -0.01332646, 0.01040754  // hy0, hy2, hy4, hy_uv, hy_u2v
 };

 public static void main(String[] args) throws Exception {
     if (args.length < 2) {
         System.out.println("Usage: java HemiPairsDefishV3 <input.jpg> <output.png> [outW outH]");
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
     final double R  = H * 0.5;     // normalize with half-height (preserves aspect ratio)

     final double ocx = outW * 0.5, ocy = outH * 0.5;
     final double Rout = outH * 0.5;

     BufferedImage dst = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);

     // One or two Newton steps are plenty (the field is smooth & monotone).
     final int NEWTON_STEPS = 2;
     final double EPS = 1e-3;

     int[] line = new int[outW];

     for (int y = 0; y < outH; y++) {
         double yn = (y - ocy) / Rout;

         for (int x = 0; x < outW; x++) {
             double xn = (x - ocx) / Rout;

             // initial guess: identity
             double u = xn, v = yn;

             for (int it = 0; it < NEWTON_STEPS; it++) {
                 double[] f0 = forward(u, v);           // (x', y')
                 double ex = f0[0] - xn, ey = f0[1] - yn;

                 // Jacobian by finite differences (small)
                 double[] fx = forward(u + EPS, v);
                 double[] fy = forward(u, v + EPS);
                 double Jxx = (fx[0] - f0[0]) / EPS, Jxy = (fy[0] - f0[0]) / EPS;
                 double Jyx = (fx[1] - f0[1]) / EPS, Jyy = (fy[1] - f0[1]) / EPS;

                 double det = Jxx * Jyy - Jxy * Jyx;
                 if (Math.abs(det) < 1e-6) break;       // very rare near extremes

                 // Newton step to reduce error
                 double du = (-ex * Jyy + ey * Jxy) / det;
                 double dv = (-ey * Jxx + ex * Jyx) / det;
                 u += du; v += dv;
             }

             double sx = cx + u * R;
             double sy = cy + v * R;
             line[x] = sampleBilinearClamp(src, sx, sy);
         }
         dst.setRGB(0, y, outW, 1, line, 0, outW);
     }
     return dst;
 }

 // Forward mapping from normalized fisheye coords (u,v) to normalized output (x',y').
 private static double[] forward(double u, double v) {
     double axu = Math.abs(u);
     double v2  = v * v;
     double r2  = u * u + v * v;

     // phiX
     double[] phiX = new double[] {
         1.0, r2, r2*r2, axu, axu*axu, axu*axu*axu, Math.pow(axu,4),
         v2, v2*v2, axu*v2, (axu*axu)*v2, r2*axu
     };
     double scaleX = dot(aX, phiX);
     double xp = u * scaleX + (hX[0] + hX[1]*v2 + hX[2]*v2*v2);

     // psiY
     double[] psiY = new double[] {
         1.0, r2, r2*r2, axu, axu*axu, axu*axu*axu, v2, v2*v2, axu*v2, (axu*axu)*v2
     };
     double scaleY = dot(aY, psiY);
     double yp = v * scaleY
             + (hY[0] + hY[1]*v2 + hY[2]*v2*v2 + hY[3]*axu*v2 + hY[4]*(axu*axu)*v2);

     return new double[]{xp, yp};
 }

 private static double dot(double[] a, double[] b){
     double s = 0; for (int i=0;i<a.length;i++) s += a[i]*b[i]; return s;
 }

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
