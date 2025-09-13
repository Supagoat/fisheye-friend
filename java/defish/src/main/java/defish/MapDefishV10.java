package defish;

//MapDefishV10.java
//Highest-quality mapping-points fisheye to hemi warp using:
//- Inverse Thin-Plate Splines (target->source) for geometry
//- EWA (Elliptical Weighted Average) Gaussian resampling (default) or Catmull–Rom bicubic
//
//JSON format: matches HemiPointMapper (src/dst with x,y and optional xn,yn).
//
//Build:
//javac MapDefishV10.java
//Run:
//java MapDefishV10 <mapping.json> <sourceImage> <output.png>
//[--lambda=1e-5] [--sigma=1.0] [--resample=ewa|bicubic] [--outW=W --outH=H] [--clampEdge]

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.*;
 

// This one is broken because it uses pixel and not relative coordinates
public class MapDefishV10 {

 public static void main(String[] args) throws Exception {
     if (args.length < 3) {
         System.out.println("Usage: java MapDefishV10 <mapping.json> <sourceImage> <output.png> [--lambda=1e-5] [--sigma=1.0] [--resample=ewa|bicubic] [--outW=W --outH=H] [--clampEdge]");
         return;
     }
     String mapPath = args[0];
     String srcPath = args[1];
     String outPath = args[2];

     double lambda = 1e-5;
     double sigma = 1.0;
     boolean clampEdge = false;
     String resample = "ewa"; // or "bicubic"
     Integer outWArg = null, outHArg = null;

     for (int i = 3; i < args.length; i++) {
         String s = args[i];
         if (s.startsWith("--lambda=")) lambda = Double.parseDouble(s.substring(9));
         else if (s.startsWith("--sigma=")) sigma = Double.parseDouble(s.substring(8));
         else if (s.equals("--clampEdge")) clampEdge = true;
         else if (s.startsWith("--resample=")) resample = s.substring(11).toLowerCase(Locale.US);
         else if (s.startsWith("--outW=")) outWArg = Integer.parseInt(s.substring(7));
         else if (s.startsWith("--outH=")) outHArg = Integer.parseInt(s.substring(7));
     }

     Mapping map = Mapping.read(new File(mapPath));
     BufferedImage srcIn = ImageIO.read(new File(srcPath));
     if (srcIn == null) throw new IOException("Cannot read source image: " + srcPath);
     BufferedImage src = toIntARGB(srcIn);

     // Choose output size: mapping target unless overridden
     int outW = (outWArg != null) ? outWArg : map.tW;
     int outH = (outHArg != null) ? outHArg : map.tH;

     // Fit inverse warp (target->source) in pixel space
     System.out.println("Fitting TPS (target->source), lambda=" + lambda + " ...");
     TPS2D tps = TPS2D.fitInverse(map, lambda);

     // Render
     System.out.println("Rendering " + outW + "x" + outH + " with " + (resample.equals("ewa") ? "EWA Gaussian (sigma="+sigma+")" : "Catmull–Rom bicubic") + (clampEdge ? ", clampEdge" : ", alpha outside"));
     BufferedImage out = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     if ("ewa".equals(resample)) {
         Renderer.renderEWA(src, out, tps, sigma, clampEdge);
     } else if ("bicubic".equals(resample)) {
         Renderer.renderBicubic(src, out, tps, clampEdge);
     } else {
         throw new IllegalArgumentException("Unknown --resample= " + resample);
     }

     ImageIO.write(out, outPath.toLowerCase(Locale.US).endsWith(".jpg") ? "jpg" : "png", new File(outPath));
     System.out.println("Wrote " + outPath);
 }

 // ---------------- Mapping / JSON ----------------
 static class Mapping {
     int sW, sH, tW, tH;
     ArrayList<Pair> pairs = new ArrayList<>();
     static class Pair { double sx, sy, dx, dy; }

     static Mapping read(File f) throws IOException {
         String text = readAll(f);
         Mapping m = new Mapping();
         m.sW = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.sH = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");
         m.tW = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.tH = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");

         // Pairs: prefer normalized if present (xn,yn); fallback to absolute x,y
         Pattern p = Pattern.compile(
           "\\{\\s*\"id\"\\s*:\\s*\\d+[^}]*?\"src\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([\\-0-9.]+))?[^}]*?\\}\\s*,\\s*\"dst\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([\\-0-9.]+))?",
           Pattern.DOTALL);
         Matcher mtr = p.matcher(text);
         int count = 0;
         while (mtr.find()) {
             Pair pr = new Pair();
             double sx = Double.parseDouble(mtr.group(1));
             double sy = Double.parseDouble(mtr.group(2));
             String sxn = mtr.group(3), syn = mtr.group(4);
             double dx = Double.parseDouble(mtr.group(5));
             double dy = Double.parseDouble(mtr.group(6));
             String dxn = mtr.group(7), dyn = mtr.group(8);

             // normalize-aware absolute positions
             if (sxn != null && syn != null) {
                 pr.sx = Double.parseDouble(sxn) * m.sW;
                 pr.sy = Double.parseDouble(syn) * m.sH;
             } else {
                 pr.sx = sx; pr.sy = sy;
             }
             if (dxn != null && dyn != null) {
                 pr.dx = Double.parseDouble(dxn) * m.tW;
                 pr.dy = Double.parseDouble(dyn) * m.tH;
             } else {
                 pr.dx = dx; pr.dy = dy;
             }
             m.pairs.add(pr); count++;
         }
         if (count == 0) throw new IOException("No pairs found in mapping JSON.");
         System.out.println("Mapping loaded: " + count + " pairs. Source " + m.sW + "x" + m.sH + "  Target " + m.tW + "x" + m.tH);
         return m;
     }
     static String readAll(File f) throws IOException {
         try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(f), StandardCharsets.UTF_8))) {
             StringBuilder sb = new StringBuilder(1<<20);
             String line; while ((line = br.readLine()) != null) sb.append(line).append('\n');
             return sb.toString();
         }
     }
     static int readInt(String text, String regex) throws IOException {
         Matcher m = Pattern.compile(regex).matcher(text);
         if (!m.find()) throw new IOException("Missing field: " + regex);
         return Integer.parseInt(m.group(1));
     }
 }

 // ---------------- TPS (2D: X and Y channels) ----------------
 static class TPS2D {
     final ThinPlateSpline tx; // X(dst) -> srcX
     final ThinPlateSpline ty; // Y(dst) -> srcY

     TPS2D(ThinPlateSpline tx, ThinPlateSpline ty) { this.tx=tx; this.ty=ty; }

     static TPS2D fitInverse(Mapping m, double lambda) {
         int n = m.pairs.size();
         double[] xd = new double[n], yd = new double[n], xs = new double[n], ys = new double[n];
         for (int i=0;i<n;i++){
             Mapping.Pair p = m.pairs.get(i);
             xd[i]=p.dx; yd[i]=p.dy; xs[i]=p.sx; ys[i]=p.sy;
         }
         ThinPlateSpline tpsX = ThinPlateSpline.fit(xd, yd, xs, lambda);
         ThinPlateSpline tpsY = ThinPlateSpline.fit(xd, yd, ys, lambda);
         return new TPS2D(tpsX, tpsY);
     }

     // Map a dst point to src point
     void map(double x, double y, double[] out) {
         out[0] = tx.eval(x,y);
         out[1] = ty.eval(x,y);
     }

     // Jacobian d(src)/d(dst) : 2x2 matrix
     void jacobian(double x, double y, double[][] J) {
         tx.grad(x,y, J[0]); // (du/dx, du/dy)
         ty.grad(x,y, J[1]); // (dv/dx, dv/dy)
     }
 }

 static class ThinPlateSpline {
     final double[] x, y;   // control points (domain)
     final double[] w;      // radial weights
     final double a0, ax, ay; // affine

     private ThinPlateSpline(double[] x, double[] y, double[] w, double a0, double ax, double ay) {
         this.x=x; this.y=y; this.w=w; this.a0=a0; this.ax=ax; this.ay=ay;
     }

     static ThinPlateSpline fit(double[] xi, double[] yi, double[] zi, double lambda) {
         int n = xi.length;
         int dim = n + 3;
         double[][] A = new double[dim][dim];
         double[]   b = new double[dim];

         // K + lambda*I
         for (int i=0;i<n;i++){
             for (int j=i;j<n;j++){
                 double r = Math.hypot(xi[i]-xi[j], yi[i]-yi[j]);
                 double v = U(r);
                 if (i==j) v += lambda;
                 A[i][j]=v; A[j][i]=v;
             }
         }
         // P blocks
         for (int i=0;i<n;i++){
             A[i][n+0]=1;       A[n+0][i]=1;
             A[i][n+1]=xi[i];   A[n+1][i]=xi[i];
             A[i][n+2]=yi[i];   A[n+2][i]=yi[i];
         }
         // RHS
         for (int i=0;i<n;i++) b[i]=zi[i];

         // Solve
         double[] coeff = solve(A,b);
         double[] w = new double[n];
         System.arraycopy(coeff, 0, w, 0, n);
         double a0 = coeff[n+0], ax = coeff[n+1], ay = coeff[n+2];
         return new ThinPlateSpline(xi.clone(), yi.clone(), w, a0, ax, ay);
     }

     double eval(double px, double py) {
         double s = a0 + ax*px + ay*py;
         for (int i=0;i<x.length;i++){
             double dx = px - x[i], dy = py - y[i];
             double r = Math.hypot(dx, dy);
             s += w[i] * U(r);
         }
         return s;
     }

     // Gradient: (diff/diffx, diff/diffy)
     void grad(double px, double py, double[] out2) {
         double gx = ax, gy = ay;
         for (int i=0;i<x.length;i++){
             double dx = px - x[i], dy = py - y[i];
             double r2 = dx*dx + dy*dy;
             if (r2 > 1e-24) {
                 double r = Math.sqrt(r2);
                 double f = (2.0*Math.log(r) + 1.0);
                 gx += w[i] * f * dx;
                 gy += w[i] * f * dy;
             }
             // for r approx 0, contribution tends to 0; skip
         }
         out2[0]=gx; out2[1]=gy;
     }

     static double U(double r) {
         if (r <= 1e-12) return 0.0;
         double r2 = r*r;
         return r2 * Math.log(r);
     }

     static double[] solve(double[][] A, double[] b) {
         int n = b.length;
         double[][] M = new double[n][n];
         for (int i=0;i<n;i++) System.arraycopy(A[i],0,M[i],0,n);
         double[] x = new double[n];
         double[] rhs = b.clone();

         for (int k=0;k<n;k++){
             int piv = k; double max = Math.abs(M[k][k]);
             for (int i=k+1;i<n;i++){
                 double v = Math.abs(M[i][k]);
                 if (v>max){ max=v; piv=i; }
             }
             if (max < 1e-14) throw new RuntimeException("Singular TPS system at k="+k);
             if (piv!=k){
                 double[] tmp = M[k]; M[k]=M[piv]; M[piv]=tmp;
                 double t = rhs[k]; rhs[k]=rhs[piv]; rhs[piv]=t;
             }
             double akk = M[k][k];
             for (int i=k+1;i<n;i++){
                 double f = M[i][k]/akk;
                 rhs[i] -= f*rhs[k];
                 for (int j=k;j<n;j++) M[i][j] -= f*M[k][j];
             }
         }
         for (int i=n-1;i>=0;i--){
             double s = rhs[i];
             for (int j=i+1;j<n;j++) s -= M[i][j]*x[j];
             x[i] = s / M[i][i];
         }
         return x;
     }
 }

 // ---------------- Rendering (EWA Gaussian or Bicubic) ----------------
 static class Renderer {

     static void renderEWA(BufferedImage src, BufferedImage dst, TPS2D tps, double sigma, boolean clampEdge) {
         final int sw = src.getWidth(), sh = src.getHeight();
         final int dw = dst.getWidth(), dh = dst.getHeight();
         final int[] sdata = ((DataBufferInt)src.getRaster().getDataBuffer()).getData();
         final int[] ddata = ((DataBufferInt)dst.getRaster().getDataBuffer()).getData();

         double[] map = new double[2];
         double[][] J = new double[2][2];
         final double kRadius = 3.0; // ~3sigma support

         for (int y=0; y<dh; y++) {
             int row = y*dw;
             for (int x=0; x<dw; x++) {
                 // 1) Inverse map: dst (x,y) -> src (sx,sy)
                 tps.map(x, y, map);
                 double sx = map[0], sy = map[1];

                 // 2) Local Jacobian J = d(src)/d(dst)
                 tps.jacobian(x, y, J);
                 double a = J[0][0], b = J[0][1]; // du/dx, du/dy
                 double c = J[1][0], d = J[1][1]; // dv/dx, dv/dy

                 // sigma = sigma^2 * J * J^T
                 double sxx = sigma*sigma*(a*a + b*b);
                 double sxy = sigma*sigma*(a*c + b*d);
                 double syy = sigma*sigma*(c*c + d*d);

                 // Invert sigma (Q = sigma^{-1}) for ellipse weight; fallback to isotropic if near-singular
                 double det = sxx*syy - sxy*sxy;
                 double A, B, C; // weight exponent = 0.5 * [A dx^2 + 2B dx dy + C dy^2]
                 if (det <= 1e-12 || !Double.isFinite(det)) {
                     // tiny footprint -> use small isotropic Gaussian
                     double inv = 1.0 / (sigma*sigma + 1e-9);
                     A = inv; B = 0; C = inv;
                 } else {
                     double invDet = 1.0 / det;
                     // Q = sigma^{-1}
                     double qxx =  syy * invDet;
                     double qxy = -sxy * invDet;
                     double qyy =  sxx * invDet;
                     A = qxx; B = qxy; C = qyy;
                 }

                 // 3) Conservative bounding box (axis-aligned)
                 int rx = (int)Math.ceil(kRadius * Math.sqrt(Math.max(sxx, 1e-12)));
                 int ry = (int)Math.ceil(kRadius * Math.sqrt(Math.max(syy, 1e-12)));
                 int xmin = (int)Math.floor(sx) - rx, xmax = (int)Math.ceil(sx) + rx;
                 int ymin = (int)Math.floor(sy) - ry, ymax = (int)Math.ceil(sy) + ry;

                 double ws = 0, wr = 0, wg = 0, wb = 0, wa = 0; // premultiplied accumulators

                 for (int yy=ymin; yy<=ymax; yy++) {
                     double dy0 = (yy + 0.5) - sy;
                     for (int xx=xmin; xx<=xmax; xx++) {
                         double dx0 = (xx + 0.5) - sx;
                         // Gaussian weight: exp(-0.5 * [A dx^2 + 2B dx dy + C dy^2])
                         double quad = A*dx0*dx0 + 2*B*dx0*dy0 + C*dy0*dy0;
                         if (quad > kRadius*kRadius*2.0) continue; // tiny weight beyond 3sigma (threshold)
                         double w = Math.exp(-0.5 * quad);

                         int px, py;
                         if (clampEdge) {
                             px = clamp(xx, 0, sw-1);
                             py = clamp(yy, 0, sh-1);
                         } else {
                             if (xx < 0 || yy < 0 || xx >= sw || yy >= sh) continue;
                             px = xx; py = yy;
                         }
                         int c0 = sdata[py*sw + px];
                         double a0 = ((c0 >>> 24) & 255) / 255.0;
                         double r0 = ((c0 >>> 16) & 255) / 255.0;
                         double g0 = ((c0 >>> 8) & 255) / 255.0;
                         double b0 = ( c0         & 255) / 255.0;

                         ws += w;
                         wa += w * a0;
                         wr += w * r0 * a0;
                         wg += w * g0 * a0;
                         wb += w * b0 * a0;
                     }
                 }

                 int outARGB;
                 if (wa <= 1e-12 || ws <= 1e-12) {
                     outARGB = 0x00000000;
                 } else {
                     double a1 = wa / ws;
                     double r1 = (wr / ws) / (a1 + 1e-12);
                     double g1 = (wg / ws) / (a1 + 1e-12);
                     double b1 = (wb / ws) / (a1 + 1e-12);
                     int ia = clamp((int)Math.round(a1 * 255.0), 0, 255);
                     int ir = clamp((int)Math.round(r1 * 255.0), 0, 255);
                     int ig = clamp((int)Math.round(g1 * 255.0), 0, 255);
                     int ib = clamp((int)Math.round(b1 * 255.0), 0, 255);
                     outARGB = (ia<<24)|(ir<<16)|(ig<<8)|ib;
                 }
                 ddata[row + x] = outARGB;
             }
         }
     }

     static void renderBicubic(BufferedImage src, BufferedImage dst, TPS2D tps, boolean clampEdge) {
         final int sw = src.getWidth(), sh = src.getHeight();
         final int dw = dst.getWidth(), dh = dst.getHeight();
         final int[] sdata = ((DataBufferInt)src.getRaster().getDataBuffer()).getData();
         final int[] ddata = ((DataBufferInt)dst.getRaster().getDataBuffer()).getData();

         double[] map = new double[2];
         for (int y=0; y<dh; y++) {
             int row = y*dw;
             for (int x=0; x<dw; x++) {
                 tps.map(x, y, map);
                 int argb = sampleBicubicARGB(sdata, sw, sh, map[0], map[1], clampEdge);
                 ddata[row + x] = argb;
             }
         }
     }

     // Catmull–Rom bicubic (no bilinear anywhere)
     static int sampleBicubicARGB(int[] data, int w, int h, double x, double y, boolean clamp) {
         if (!clamp && (x < -0.5 || y < -0.5 || x > w-0.5 || y > h-0.5)) return 0x00000000;

         int x1 = (int)Math.floor(x);
         int y1 = (int)Math.floor(y);
         double tx = x - x1;
         double ty = y - y1;

         int[] xs = new int[4], ys = new int[4];
         for (int i=-1;i<=2;i++){ xs[i+1] = clamp? clamp(x1+i,0,w-1) : (x1+i); }
         for (int j=-1;j<=2;j++){ ys[j+1] = clamp? clamp(y1+j,0,h-1) : (y1+j); }

         double[] wx = catmullRomWeights(tx);
         double[] wy = catmullRomWeights(ty);

         double a=0,r=0,g=0,b=0,W=0;
         for (int j=0;j<4;j++){
             int yy = ys[j]; if (!clamp && (yy<0 || yy>=h)) continue;
             for (int i=0;i<4;i++){
                 int xx = xs[i]; if (!clamp && (xx<0 || xx>=w)) continue;
                 double wgt = wx[i]*wy[j];
                 int c = data[yy*w + xx];
                 double ca = ((c>>>24)&255)/255.0;
                 double cr = ((c>>>16)&255)/255.0;
                 double cg = ((c>>>8 )&255)/255.0;
                 double cb = ( c      &255)/255.0;
                 a += wgt*ca; r += wgt*cr*ca; g += wgt*cg*ca; b += wgt*cb*ca; W += wgt;
             }
         }
         if (W <= 1e-12 || a <= 1e-12) return 0x00000000;
         double A = a / W;
         double R = (r / W) / (A+1e-12);
         double G = (g / W) / (A+1e-12);
         double B = (b / W) / (A+1e-12);
         int ia = clamp((int)Math.round(A*255),0,255);
         int ir = clamp((int)Math.round(R*255),0,255);
         int ig = clamp((int)Math.round(G*255),0,255);
         int ib = clamp((int)Math.round(B*255),0,255);
         return (ia<<24)|(ir<<16)|(ig<<8)|ib;
     }

     static double[] catmullRomWeights(double t) {
         double a=-0.5, t2=t*t, t3=t2*t;
         double w0 = a*(-t3 + 2*t2 - t);
         double w1 = (a+2)*t3 + (-a-3)*t2 + 1;
         double w2 = (a+2)*(-t3) + (2*a+3)*t2 + (-a)*t;
         double w3 = a*(t3 - t2);
         return new double[]{w0,w1,w2,w3};
     }
 }

 // ---------------- Utils ----------------
 static BufferedImage toIntARGB(BufferedImage in) {
     if (in.getType() == BufferedImage.TYPE_INT_ARGB) return in;
     BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
     Graphics2D g = out.createGraphics();
     g.setComposite(AlphaComposite.Src);
     g.drawImage(in, 0, 0, null);
     g.dispose();
     return out;
 }
 static int clamp(int v,int lo,int hi){ return v<lo?lo:(v>hi?hi:v); }
}
