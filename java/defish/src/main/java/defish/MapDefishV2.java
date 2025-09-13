package defish;

//HemiPairsDefishV2.java
//Data-driven defish from paired correspondences (mirrored in X & Y about center).
//Normalization by half-height preserves aspect ratio; verticals remain straight by construction.
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;

public class MapDefishV2 {

 // ====== MIRROR-SYMMETRIC FIELD (expanded basis; corner-weighted fit) ======
 // Basis theta(|u|, v^2, r^2) = {
 //   1, r2, r2^2, r2^3, r2^4, r2^5,
 //   |u|, |u|^2, |u|^3, |u|^4,
 //   v2, v2^2, v2^3, v2^4,
 //   |u|v2, |u|v2^2, |u|v2^3, |u|v2^4,
 //   |u|^2 v2, |u|^2 v2^2, |u|^2 v2^3,
 //   |u|^3 v2, |u|^3 v2^2, |u|^3 v2^3,
 //   |u|^4 v2, |u|^4 v2^2,
 //   r2|u|, r2 v2, r2 |u| v2
 // }
 // Forward map (normalized coords, half-height scaling):
 //   x' = u * (aX · theta) + (hX0 + hX2 v2 + hX4 v2^2 + hX6 v2^3 + hX8 v2^4)
 //   y' = v * (aY · theta) + (hY0 + hY1|u| + hY2|u|^2 + hY3|u|^3 + hY2v v2 + hY4 v2^2 + hY6 v2^3 + hY8 v2^4
 //                + hYuv |u|v2 + hYu2v |u|^2 v2 + hYu3v |u|^3 v2)
 // Note: hY0 is a small global vertical bias. Set to 0 for strict top/bottom symmetry.

 private static final double[] aX = {
     0.8285546592469589, 2.942606590749195, 4.218593544573172, -4.255758794867123, 1.171519099423085,
     -0.1414968454527168, -1.024028760843289, 4.02039262321832, -7.904037550647732, 8.56324392620437,
     -1.077786032418816, -4.166542084893585, 7.541347600593055, -1.581292925528559, 0.3838950651129238,
     7.312218214176077, -11.99439727697023, 1.60268591069505, -0.08905414836212874, -5.352036440993072,
     6.513145196535659, 2.914924683632383, -1.418822691259738, -1.561103054654899, -3.004947407561674,
     0.8477692260063553, -7.520142485533939, -4.255596233261209, 10.2271428978197
 };

 private static final double[] hX = {
     -0.0002755205040804885, 0.01881571573964771, -0.09902842552693064, 0.1469748442529657, -0.06218298047753369
 };

 private static final double[] aY = {
     0.8272326162411245, 0.5968864762157713, -0.1288746923770851, -2.271726882221418, 0.743604675149687,
     -0.09910509641786652, -0.1374640710302285, 0.9730120926482431, -3.435174580536965, 5.149270281310891,
     -0.3761256161444659, 1.616294650220777, 2.726567640481252, -0.7990318390648553, 2.06794426026485,
     -4.374896407964218, -1.287931267782078, 1.02974711375137, -3.447219811982588, 13.86961850303312,
     -4.687098688326272, 8.299837238062532, -11.88109412900631, 2.873037845279539, -3.256573328789558,
     2.568941332056267, -1.367230320280919, -1.830925161746019, 3.924940830091271
 };

 private static final double[] hY = {
     0.01862110983316903, 0.001406758974484824, -0.002455949188385805, 0.0003292349642026195, -0.02654219570716532,
     0.01177767866999832, -0.006947307503569709, 0.003640527843974444, 0.04101879198684013, -0.08615568808356162,
     0.04804523387643669
 };

 // Set this to true to force strict top/bottom symmetry (zero vertical bias).
 private static final boolean ZERO_VERTICAL_BIAS = false;

 public static void main(String[] args) throws Exception {
     if (args.length < 2) {
         System.out.println("Usage: java HemiPairsDefishV2 <input.jpg> <output.png> [outW outH]");
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
     int[] line = new int[outW];

     for (int y = 0; y < outH; y++) {
         double yn = (y - ocy) / Rout;
         for (int x = 0; x < outW; x++) {
             double xn = (x - ocx) / Rout;

             // Invert forward(u,v) wavy equals (xn,yn) with one damped Newton step
             double u = xn, v = yn;
             double[] f0 = forward(u, v);
             double ex = f0[0] - xn, ey = f0[1] - yn;
             final double eps = 1e-3;
             double[] fx = forward(u + eps, v), fy = forward(u, v + eps);
             double Jxx = (fx[0] - f0[0]) / eps, Jxy = (fy[0] - f0[0]) / eps;
             double Jyx = (fx[1] - f0[1]) / eps, Jyy = (fy[1] - f0[1]) / eps;
             double det = Jxx * Jyy - Jxy * Jyx;
             if (Math.abs(det) < 1e-6) det = (det >= 0 ? 1e-6 : -1e-6);
             double du = (-ex * Jyy + ey * Jxy) / det;
             double dv = (-ey * Jxx + ex * Jyx) / det;
             u += du; v += dv;

             double sx = cx + u * R;
             double sy = cy + v * R;
             line[x] = sampleBilinearClamp(src, sx, sy);
         }
         dst.setRGB(0, y, outW, 1, line, 0, outW);
     }
     return dst;
 }

 private static double[] forward(double u, double v) {
     double ax = Math.abs(u);
     double v2 = v*v;
     double r2 = u*u + v*v;

     double[] phi = new double[]{
         1.0, r2, r2*r2, r2*r2*r2, Math.pow(r2,4), Math.pow(r2,5),
         ax, ax*ax, ax*ax*ax, Math.pow(ax,4),
         v2, v2*v2, v2*v2*v2, v2*v2*v2*v2,
         ax*v2, ax*v2*v2, ax*v2*v2*v2, ax*v2*v2*v2*v2,
         (ax*ax)*v2, (ax*ax)*v2*v2, (ax*ax)*v2*v2*v2,
         (ax*ax*ax)*v2, (ax*ax*ax)*v2*v2, (ax*ax*ax)*v2*v2*v2,
         Math.pow(ax,4)*v2, Math.pow(ax,4)*v2*v2,
         r2*ax, r2*v2, r2*ax*v2
     };

     double scaleX = dot(aX, phi);
     double xp = u * scaleX + (hX[0] + hX[1]*v2 + hX[2]*v2*v2 + hX[3]*v2*v2*v2 + hX[4]*v2*v2*v2*v2);

     double scaleY = dot(aY, phi);
     double biasY0 = ZERO_VERTICAL_BIAS ? 0.0 : hY[0];
     double yp = v * scaleY
             + ( biasY0
                 + hY[1]*ax + hY[2]*ax*ax + hY[3]*ax*ax*ax
                 + hY[4]*v2 + hY[5]*v2*v2 + hY[6]*v2*v2*v2 + hY[7]*v2*v2*v2*v2
                 + hY[8]*ax*v2 + hY[9]*ax*ax*v2 + hY[10]*ax*ax*ax*v2 );

     return new double[]{xp, yp};
 }

 private static double dot(double[] a, double[] b) {
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
