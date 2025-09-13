package defish;
//MapDefishV13.java
//Fisheye to Hemi using mapping points, TPS in normalized coords, and
//non-separable 2D radial Lanczos-3 resampling (no bilinear).
//
//Build:  javac MapDefishV13.java
//Run:    java MapDefishV13 <mapping.json> <sourceImage> <output.png>
//      [--lambda=1e-5] [--outW=W --outH=H] [--clampEdge]

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.*;

public class MapDefishV13 {

 public static void main(String[] args) throws Exception {
     if (args.length < 3) {
         System.out.println("Usage: java MapDefishV13 <mapping.json> <sourceImage> <output.png> [--lambda=1e-5] [--outW=W --outH=H] [--clampEdge]");
         return;
     }
     String mapPath = args[0];
     String srcPath = args[1];
     String outPath = args[2];

     double lambda = 1e-5;
     Integer outWArg = null, outHArg = null;
     boolean clampEdge = false;

     for (int i=3;i<args.length;i++){
         String s = args[i];
         if (s.startsWith("--lambda=")) lambda = Double.parseDouble(s.substring(9));
         else if (s.startsWith("--outW=")) outWArg = Integer.parseInt(s.substring(7));
         else if (s.startsWith("--outH=")) outHArg = Integer.parseInt(s.substring(7));
         else if (s.equals("--clampEdge")) clampEdge = true;
     }

     Mapping map = Mapping.read(new File(mapPath));

     BufferedImage srcIn = ImageIO.read(new File(srcPath));
     if (srcIn == null) throw new IOException("Cannot read source image: " + srcPath);
     BufferedImage src = toIntARGB(srcIn);
     int sw = src.getWidth(), sh = src.getHeight();

     int outW = (outWArg != null) ? outWArg : sw;
     int outH = (outHArg != null) ? outHArg : sh;

     System.out.println("Fitting TPS (normalized) lambda=" + lambda);
     TPS2D tps = TPS2D.fitInverseNormalized(map, lambda);

     System.out.println("Rendering " + outW + "x" + outH + " with radial Lanczos-3" + (clampEdge?", clampEdge":""));
     BufferedImage out = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
     Renderer.renderLanczos3Radial(src, out, tps, clampEdge);

     ImageIO.write(out, outPath.toLowerCase(Locale.US).endsWith(".jpg")?"jpg":"png", new File(outPath));
     System.out.println("Wrote " + outPath);
 }

 // ---------- Mapping (read JSON; normalized pairs) ----------
 static class Mapping {
     int sW, sH, tW, tH;
     ArrayList<Pair> pairs = new ArrayList<>();
     static class Pair { double sxN, syN, dxN, dyN; }

     static Mapping read(File f) throws IOException {
         String text = readAll(f);
         Mapping m = new Mapping();
         m.sW = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.sH = readInt(text, "\"source\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");
         m.tW = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
         m.tH = readInt(text, "\"target\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");

         Pattern p = Pattern.compile(
             "\\{\\s*\"id\"\\s*:\\s*\\d+[^}]*?\"src\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([\\-0-9.]+))?[^}]*?\\}\\s*,\\s*\"dst\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"y\"\\s*:\\s*([\\-0-9.]+)(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)[^}]*?\"yn\"\\s*:\\s*([\\-0-9.]+))?",
             Pattern.DOTALL);
         Matcher mtr = p.matcher(text);
         int n=0;
         while (mtr.find()) {
             double sx = Double.parseDouble(mtr.group(1));
             double sy = Double.parseDouble(mtr.group(2));
             String sxn = mtr.group(3), syn = mtr.group(4);
             double dx = Double.parseDouble(mtr.group(5));
             double dy = Double.parseDouble(mtr.group(6));
             String dxn = mtr.group(7), dyn = mtr.group(8);

             Pair pr = new Pair();
             pr.sxN = (sxn!=null? Double.parseDouble(sxn) : sx / m.sW);
             pr.syN = (syn!=null? Double.parseDouble(syn) : sy / m.sH);
             pr.dxN = (dxn!=null? Double.parseDouble(dxn) : dx / m.tW);
             pr.dyN = (dyn!=null? Double.parseDouble(dyn) : dy / m.tH);
             m.pairs.add(pr); n++;
         }
         if (n==0) throw new IOException("No pairs found in "+f);
         System.out.println("Loaded "+n+" pairs (normalized).");
         return m;
     }

     static String readAll(File f) throws IOException {
         try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(f), StandardCharsets.UTF_8))){
             StringBuilder sb = new StringBuilder(1<<20);
             String line; while((line=br.readLine())!=null) sb.append(line).append('\n');
             return sb.toString();
         }
     }
     static int readInt(String text, String regex) throws IOException {
         Matcher m = Pattern.compile(regex).matcher(text);
         if (!m.find()) throw new IOException("Missing field: "+regex);
         return Integer.parseInt(m.group(1));
     }
 }

 // ---------- TPS (fit/eval in normalized coords) ----------
 static class TPS2D {
     final ThinPlateSpline tx; // (x_dn,y_dn) -> x_sn
     final ThinPlateSpline ty; // (x_dn,y_dn) -> y_sn
     TPS2D(ThinPlateSpline tx, ThinPlateSpline ty){ this.tx=tx; this.ty=ty; }

     static TPS2D fitInverseNormalized(Mapping m, double lambda){
         int n = m.pairs.size();
         double[] xd = new double[n], yd = new double[n], xs = new double[n], ys = new double[n];
         for (int i=0;i<n;i++){
             Mapping.Pair p = m.pairs.get(i);
             xd[i]=p.dxN; yd[i]=p.dyN; xs[i]=p.sxN; ys[i]=p.syN;
         }
         return new TPS2D(ThinPlateSpline.fit(xd, yd, xs, lambda),
                          ThinPlateSpline.fit(xd, yd, ys, lambda));
     }
     void mapNorm(double xn, double yn, double[] out2){
         out2[0] = tx.eval(xn, yn);
         out2[1] = ty.eval(xn, yn);
     }
 }

 static class ThinPlateSpline {
     final double[] x, y, w; final double a0, ax, ay;
     private ThinPlateSpline(double[] x,double[] y,double[] w,double a0,double ax,double ay){this.x=x;this.y=y;this.w=w;this.a0=a0;this.ax=ax;this.ay=ay;}

     static ThinPlateSpline fit(double[] xi, double[] yi, double[] zi, double lambda){
         int n=xi.length, dim=n+3;
         double[][] A = new double[dim][dim];
         double[] b = new double[dim];

         for (int i=0;i<n;i++){
             for (int j=i;j<n;j++){
                 double r = Math.hypot(xi[i]-xi[j], yi[i]-yi[j]);
                 double v = U(r) + (i==j? lambda:0);
                 A[i][j]=v; A[j][i]=v;
             }
         }
         for (int i=0;i<n;i++){
             A[i][n]=1;   A[n][i]=1;
             A[i][n+1]=xi[i]; A[n+1][i]=xi[i];
             A[i][n+2]=yi[i]; A[n+2][i]=yi[i];
         }
         for (int i=0;i<n;i++) b[i]=zi[i];

         double[] sol = solve(A,b);
         double[] w = Arrays.copyOf(sol, n);
         return new ThinPlateSpline(xi.clone(), yi.clone(), w, sol[n], sol[n+1], sol[n+2]);
     }

     double eval(double px,double py){
         double s = a0 + ax*px + ay*py;
         for (int i=0;i<x.length;i++){
             double dx=px-x[i], dy=py-y[i];
             double r = Math.hypot(dx,dy);
             s += w[i]*U(r);
         }
         return s;
     }
     static double U(double r){ if (r<=1e-12) return 0; double r2=r*r; return r2*Math.log(r); }

     static double[] solve(double[][] A,double[] b){
         int n=b.length;
         double[][] M=new double[n][n]; for(int i=0;i<n;i++) System.arraycopy(A[i],0,M[i],0,n);
         double[] x=new double[n]; double[] rhs=b.clone();
         for(int k=0;k<n;k++){
             int piv=k; double mx=Math.abs(M[k][k]);
             for(int i=k+1;i<n;i++){ double v=Math.abs(M[i][k]); if(v>mx){mx=v;piv=i;} }
             if (mx<1e-14) throw new RuntimeException("Singular TPS system");
             if (piv!=k){ double[] tmp=M[k]; M[k]=M[piv]; M[piv]=tmp; double t=rhs[k]; rhs[k]=rhs[piv]; rhs[piv]=t; }
             double akk=M[k][k];
             for(int i=k+1;i<n;i++){
                 double f=M[i][k]/akk; rhs[i]-=f*rhs[k];
                 for(int j=k;j<n;j++) M[i][j]-=f*M[k][j];
             }
         }
         for(int i=n-1;i>=0;i--){
             double s=rhs[i]; for(int j=i+1;j<n;j++) s-=M[i][j]*x[j];
             x[i]=s/M[i][i];
         }
         return x;
     }
 }

 // ---------- Renderer: non-separable radial Lanczos-3 ----------
 static class Renderer {
     static void renderLanczos3Radial(BufferedImage src, BufferedImage dst, TPS2D tps, boolean clampEdge){
         final int sw=src.getWidth(), sh=src.getHeight();
         final int dw=dst.getWidth(), dh=dst.getHeight();
         final int[] sdata=((DataBufferInt)src.getRaster().getDataBuffer()).getData();
         final int[] ddata=((DataBufferInt)dst.getRaster().getDataBuffer()).getData();

         final double a = 3.0; // Lanczos radius
         double[] map = new double[2];

         for (int y=0; y<dh; y++){
             int row = y*dw;
             double yn = (y + 0.5) / dh;
             for (int x=0; x<dw; x++){
                 double xn = (x + 0.5) / dw;

                 // inverse map in normalized coords
                 tps.mapNorm(xn, yn, map);
                 double sx = map[0]*sw - 0.5;
                 double sy = map[1]*sh - 0.5;

                 // support window
                 int xmin = (int)Math.floor(sx - a + 1);
                 int xmax = (int)Math.ceil (sx + a);
                 int ymin = (int)Math.floor(sy - a + 1);
                 int ymax = (int)Math.ceil (sy + a);

                 double W=0, Aa=0, Ar=0, Ag=0, Ab=0; // premultiplied accumulation
                 for (int yy=ymin; yy<=ymax; yy++){
                     for (int xx=xmin; xx<=xmax; xx++){
                         double dx = (xx + 0.5) - sx;
                         double dy = (yy + 0.5) - sy;
                         double r = Math.hypot(dx, dy);
                         if (r >= a) continue;
                         double w = lanczos(r, a);
                         int px=xx, py=yy;
                         if (clampEdge) {
                             if (px<0) px=0; else if (px>=sw) px=sw-1;
                             if (py<0) py=0; else if (py>=sh) py=sh-1;
                         } else {
                             if (px<0||py<0||px>=sw||py>=sh) continue;
                         }
                         int c = sdata[py*sw + px];
                         double a0=((c>>>24)&255)/255.0;
                         double r0=((c>>>16)&255)/255.0;
                         double g0=((c>>>8 )&255)/255.0;
                         double b0=( c      &255)/255.0;

                         W  += w;
                         Aa += w*a0;
                         Ar += w*r0*a0;
                         Ag += w*g0*a0;
                         Ab += w*b0*a0;
                     }
                 }
                 int outARGB;
                 if (W<=1e-12 || Aa<=1e-12) {
                     outARGB = 0x00000000;
                 } else {
                     double A1 = Aa/W;
                     double R1 = (Ar/W)/(A1+1e-12);
                     double G1 = (Ag/W)/(A1+1e-12);
                     double B1 = (Ab/W)/(A1+1e-12);
                     int ia=clamp((int)Math.round(A1*255),0,255);
                     int ir=clamp((int)Math.round(R1*255),0,255);
                     int ig=clamp((int)Math.round(G1*255),0,255);
                     int ib=clamp((int)Math.round(B1*255),0,255);
                     outARGB = (ia<<24)|(ir<<16)|(ig<<8)|ib;
                 }
                 ddata[row + x] = outARGB;
             }
         }
     }

     static double sinc(double t){
         if (Math.abs(t) < 1e-8) return 1.0;
         double x = Math.PI * t;
         return Math.sin(x) / x;
     }
     static double lanczos(double r, double a){
         // non-separable 2D radial Lanczos: sinc(r) * sinc(r/a), |r|<a
         return sinc(r) * sinc(r / a);
     }
 }

 // ---------- Utils ----------
 static BufferedImage toIntARGB(BufferedImage in){
     if (in.getType()==BufferedImage.TYPE_INT_ARGB) return in;
     BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
     Graphics2D g = out.createGraphics();
     g.setComposite(AlphaComposite.Src);
     g.drawImage(in,0,0,null);
     g.dispose();
     return out;
 }
 static int clamp(int v,int lo,int hi){ return v<lo?lo:(v>hi?hi:v); }
}
