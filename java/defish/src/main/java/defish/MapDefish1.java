package defish;

//HemiPairsDefish.java
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;

/**
* Applies a field fitted from your point-pair mapping (334 pairs).
* Normalization is by half-height (preserves aspect ratio).
*
* Usage:
*   javac HemiPairsDefish.java
*   java HemiPairsDefish input.jpg output.png [outW outH]
*
* Fit quality (against your pairs): ~4.76 px X, 1.57 px Y RMS at H=1197.
*/
public class MapDefish1 {

 // Basis Φ(x,y): [1, r2, r2^2, r2^3, |x|, |x|^2, |x|^3,
 //                y^2, y^4, y^6,
 //                |x|y^2, |x|y^4, |x|y^6,
 //                |x|^2 y^2, |x|^2 y^4, |x|^2 y^6,
 //                |x|^3 y^2, |x|^3 y^4, |x|^3 y^6]
 // Horizontal: x' = x*(aX·Φ) + (hX0 + hX2*y^2 + hX4*y^4)
 // Vertical:   y' = y*(aY·Φ) + (hY0 + hY1|x| + hY2|x|^2 + hY3|x|^3
 //                               + hY4*y^2 + hY5*y^4 + hY6|x|y^2 + hY7|x|^2 y^2)

 private static final double[] aX = {
     0.8229629634436182, 0.15771353270406835, 0.0619970741012078, -0.004332170903247524,
     0.040252249966374706, -0.09842569080641993, -0.0828400220152391, 0.2561392235052598,
     -1.0613629360949446, 1.021226761182936, -1.0406296360026555, 4.827909326644998,
     -4.600359456662002, 1.0894651034831662, -6.048402614189964, 5.997219673871758,
     -0.3841903584142022, 2.310227248335149, -2.3798606645441343
 };
 private static final double[] hX = {
     -0.0005996671893741431, 0.001654356146693126, -0.00019956548081666038
 };
 private static final double[] aY = {
     0.8338486422427539, 0.10142636911812475, 0.1693436161613262, -0.02102642390807416,
     -0.02086495575119012, 0.03748628412354785, -0.2684159153759164, 0.06394008499285586,
     -0.26029241298686917, 0.12101141419138596, -0.13009855787645092, 0.5592677537158074,
     -0.46709577543810227, -0.15788555025990197, -0.9617914126870865, 0.8554809183738011,
     0.014855355817333774, 0.493727498146218, -0.42218383234222184
 };
 private static final double[] hY = {
     0.018181644204259523, 0.0014545236031722945, -0.004883212588492344,
     0.002302570993848372, -0.014873007919680364, -0.0012917863289990876,
     -0.003237958140917172, 0.004185829168845514
 };

 public static void main(String[] args) throws Exception {
     if (args.length < 2) {
         System.out.println("Usage: java HemiPairsDefish <input.jpg> <output.png> [outW outH]");
         System.out.println("Tip: for the cloister pair, use 1769 1197 (the Hemi canvas).");
     }
     String inPath  = args.length >= 1 ? args[0] : null;
     String outPath = args.length >= 2 ? args[1] : "hemi_like.png";
     int outW = (args.length >= 4) ? Integer.parseInt(args[2]) : -1;
     int outH = (args.length >= 4) ? Integer.parseInt(args[3]) : -1;

     BufferedImage src = ImageIO.read(new File(inPath));
     if (src == null) throw new RuntimeException("Cannot read " + inPath);
     if (outW <= 0) outW = src.getWidth();
     if (outH <= 0) outH = src.getHeight();

     BufferedImage dst = remap(src, outW, outH);
     ImageIO.write(dst, outPath.toLowerCase().endsWith(".jpg") ? "jpg" : "png", new File(outPath));
     System.out.println("Wrote " + outPath + " (" + outW + "x" + outH + ")");
 }

 public static BufferedImage remap(BufferedImage src, int outW, int outH) {
     final int W = src.getWidth(), H = src.getHeight();
     final double cx = W * 0.5, cy = H * 0.5;
     final double R  = H * 0.5;     // normalize with half-height

     final double ocx = outW * 0.5, ocy = outH * 0.5;
     final double Rout = outH * 0.5;

     BufferedImage dst = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     int[] line = new int[outW];

     for (int y = 0; y < outH; y++) {
         double yn = (y - ocy) / Rout;
         for (int x = 0; x < outW; x++) {
             double xn = (x - ocx) / Rout;

             // Invert the forward field with one Newton step
             double u = xn, v = yn;                 // initial guess
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

             // Map to source pixels
             double sx = cx + u * R;
             double sy = cy + v * R;
             line[x] = sampleBilinearClamp(src, sx, sy);
         }
         dst.setRGB(0, y, outW, 1, line, 0, outW);
     }
     return dst;
 }

 // forward: source-like normalized (u,v) -> normalized output (x',y')
 private static double[] forward(double u, double v) {
     double ax = Math.abs(u), yy = v*v, r2 = u*u + v*v;

     double[] phi = new double[19];
     int k = 0;
     phi[k++] = 1.0;
     phi[k++] = r2;             phi[k++] = r2*r2;            phi[k++] = r2*r2*r2;
     phi[k++] = ax;             phi[k++] = ax*ax;            phi[k++] = ax*ax*ax;
     phi[k++] = yy;             phi[k++] = yy*yy;            phi[k++] = yy*yy*yy;
     phi[k++] = ax*yy;          phi[k++] = ax*yy*yy;         phi[k++] = ax*yy*yy*yy;
     phi[k++] = ax*ax*yy;       phi[k++] = ax*ax*yy*yy;      phi[k++] = ax*ax*yy*yy*yy;
     phi[k++] = ax*ax*ax*yy;    phi[k++] = ax*ax*ax*yy*yy;   phi[k++] = ax*ax*ax*yy*yy*yy;

     double scaleX = dot(aX, phi);
     double xp = u * scaleX + (hX[0] + hX[1]*yy + hX[2]*yy*yy);

     double scaleY = dot(aY, phi);
     double yp = v * scaleY + (hY[0] + hY[1]*ax + hY[2]*ax*ax + hY[3]*ax*ax*ax
                               + hY[4]*yy + hY[5]*yy*yy + hY[6]*ax*yy + hY[7]*ax*ax*yy);
     return new double[]{xp, yp};
 }

 private static double dot(double[] a, double[] b) {
     double s = 0; for (int i=0;i<a.length;i++) s += a[i]*b[i]; return s;
 }

 // Bilinear sampling with clamp
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
