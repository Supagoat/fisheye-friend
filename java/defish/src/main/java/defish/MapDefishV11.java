package defish;

// MapDefishV11.java
// High-quality fisheye to hemi via mapping points using TPS in normalized space.
// Resampling: Catmull–Rom bicubic (default, high contrast) or EWA Gaussian (antialiasing).
//
// Usage:
//   javac MapDefishV11.java
//   java MapDefishV11 <mapping.json> <sourceImage> <output.png>
//     [--lambda=1e-5] [--resample=bicubic|ewa] [--sigma=0.9] [--outW=W --outH=H] [--clampEdge]

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.*;

// This one is good but somewhat low contrast

public class MapDefishV11 {

    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.out.println("Usage: java MapDefishV11 <mapping.json> <sourceImage> <output.png> [--lambda=1e-5] [--resample=bicubic|ewa] [--sigma=0.9] [--outW=W --outH=H] [--clampEdge]");
            return;
        }
        String mapPath = args[0];
        String srcPath = args[1];
        String outPath = args[2];

        double lambda = 1e-5;
        String resample = "bicubic";  // crisper than bilinear
        double sigma = 0.9;
        boolean clampEdge = false;
        Integer outWArg = null, outHArg = null;

        for (int i=3;i<args.length;i++){
            String s = args[i];
            if (s.startsWith("--lambda=")) lambda = Double.parseDouble(s.substring(9));
            else if (s.startsWith("--resample=")) resample = s.substring(11).toLowerCase(Locale.US);
            else if (s.startsWith("--sigma=")) sigma = Double.parseDouble(s.substring(8));
            else if (s.startsWith("--outW=")) outWArg = Integer.parseInt(s.substring(7));
            else if (s.startsWith("--outH=")) outHArg = Integer.parseInt(s.substring(7));
            else if (s.equals("--clampEdge")) clampEdge = true;
        }

        Mapping map = Mapping.read(new File(mapPath));

        BufferedImage srcIn = ImageIO.read(new File(srcPath));
        if (srcIn == null) throw new IOException("Cannot read source image: " + srcPath);
        BufferedImage src = toIntARGB(srcIn);
        int sw = src.getWidth(), sh = src.getHeight();

        // Output size defaults to the source size (so you can drop it into your flow)
        int outW = (outWArg != null) ? outWArg : sw;
        int outH = (outHArg != null) ? outHArg : sh;

        // Fit inverse TPS in normalized domain: (dst_n) -> (src_n)
        System.out.println("Fitting TPS in normalized space (lambda="+lambda+") …");
        TPS2D tps = TPS2D.fitInverseNormalized(map, lambda);

        // Render
        System.out.println("Rendering " + outW + "x" + outH + "  resampler=" + resample + (resample.equals("ewa")?(" sigma="+sigma):"") + (clampEdge?", clampEdge":""));
        BufferedImage out = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);
        if ("ewa".equals(resample)) {
            Renderer.renderEWA(src, out, tps, sigma, clampEdge);
        } else if ("bicubic".equals(resample)) {
            Renderer.renderBicubic(src, out, tps, clampEdge);
        } else {
            throw new IllegalArgumentException("Unknown --resample=" + resample);
        }

        ImageIO.write(out, outPath.toLowerCase(Locale.US).endsWith(".jpg")?"jpg":"png", new File(outPath));
        System.out.println("Wrote " + outPath);
    }

    // ---------- Mapping (reads JSON; builds normalized pairs) ----------
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

    // ---------- TPS (fit and evaluate in normalized coords) ----------
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

        // Evaluate normalized map
        void mapNorm(double xn, double yn, double[] out2){
            out2[0] = tx.eval(xn, yn);
            out2[1] = ty.eval(xn, yn);
        }
        // Jacobian in normalized domain
        void jacobianNorm(double xn, double yn, double[][] J2){
            tx.grad(xn, yn, J2[0]); // d(x_sn)/d(x_dn), d(x_sn)/d(y_dn)
            ty.grad(xn, yn, J2[1]); // d(y_sn)/d(x_dn), d(y_sn)/d(y_dn)
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
        // Gradient (diff/diffx, diff/diffy)
        void grad(double px,double py,double[] out2){
            double gx=ax, gy=ay;
            for (int i=0;i<x.length;i++){
                double dx=px-x[i], dy=py-y[i];
                double r2 = dx*dx+dy*dy;
                if (r2>1e-24){
                    double r = Math.sqrt(r2);
                    double f = (2.0*Math.log(r) + 1.0); // d/dr of r^2 log r projected on x,y
                    gx += w[i]*f*dx;
                    gy += w[i]*f*dy;
                }
            }
            out2[0]=gx; out2[1]=gy;
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

    // ---------- Renderers ----------
    static class Renderer {

        static void renderBicubic(BufferedImage src, BufferedImage dst, TPS2D tps, boolean clampEdge){
            int sw=src.getWidth(), sh=src.getHeight();
            int dw=dst.getWidth(), dh=dst.getHeight();
            int[] sdata=((DataBufferInt)src.getRaster().getDataBuffer()).getData();
            int[] ddata=((DataBufferInt)dst.getRaster().getDataBuffer()).getData();
            double[] map = new double[2];

            for(int y=0;y<dh;y++){
                int row=y*dw;
                double yn = ((y+0.5)/dh);
                for(int x=0;x<dw;x++){
                    double xn = ((x+0.5)/dw);
                    tps.mapNorm(xn, yn, map);
                    double sx = map[0] * sw - 0.5;
                    double sy = map[1] * sh - 0.5;
                    ddata[row+x] = sampleBicubicARGB(sdata, sw, sh, sx, sy, clampEdge);
                }
            }
        }

        static void renderEWA(BufferedImage src, BufferedImage dst, TPS2D tps, double sigma, boolean clampEdge){
            int sw=src.getWidth(), sh=src.getHeight();
            int dw=dst.getWidth(), dh=dst.getHeight();
            int[] sdata=((DataBufferInt)src.getRaster().getDataBuffer()).getData();
            int[] ddata=((DataBufferInt)dst.getRaster().getDataBuffer()).getData();

            double[] map = new double[2];
            double[][] Jn = new double[2][2];
            final double kRad = 3.0;

            for(int y=0;y<dh;y++){
                int row=y*dw;
                double ydn = ((y+0.5)/dh);
                for(int x=0;x<dw;x++){
                    double xdn = ((x+0.5)/dw);

                    // map & Jacobian in normalized domain
                    tps.mapNorm(xdn, ydn, map);
                    tps.jacobianNorm(xdn, ydn, Jn);

                    // convert to source pixel units:  J_pix = S_src * J_norm * S_dst^{-1}
                    double a = Jn[0][0] * sw / dw;  // du/dx
                    double b = Jn[0][1] * sw / dh;  // du/dy
                    double c = Jn[1][0] * sh / dw;  // dv/dx
                    double d = Jn[1][1] * sh / dh;  // dv/dy

                    double sx = map[0]*sw - 0.5;
                    double sy = map[1]*sh - 0.5;

                    // sigma = sigma^2 * J * J^T
                    double sxx = sigma*sigma*(a*a + b*b);
                    double sxy = sigma*sigma*(a*c + b*d);
                    double syy = sigma*sigma*(c*c + d*d);

                    double det = sxx*syy - sxy*sxy;
                    double A,B,C;
                    if (det <= 1e-14 || !Double.isFinite(det)) {
                        double inv = 1.0 / (sigma*sigma + 1e-9);
                        A=inv; B=0; C=inv;
                    } else {
                        double invDet = 1.0/det;
                        double qxx =  syy*invDet;
                        double qxy = -sxy*invDet;
                        double qyy =  sxx*invDet;
                        A=qxx; B=qxy; C=qyy;
                    }

                    // conservative bbox around ellipse
                    int rx = (int)Math.ceil(kRad * Math.sqrt(Math.max(sxx, 1e-12)));
                    int ry = (int)Math.ceil(kRad * Math.sqrt(Math.max(syy, 1e-12)));
                    int xmin = (int)Math.floor(sx) - rx, xmax = (int)Math.ceil(sx) + rx;
                    int ymin = (int)Math.floor(sy) - ry, ymax = (int)Math.ceil(sy) + ry;

                    double W=0, Ra=0, Rr=0, Rg=0, Rb=0;
                    for (int yy=ymin; yy<=ymax; yy++){
                        double dy0 = (yy+0.5) - sy;
                        for (int xx=xmin; xx<=xmax; xx++){
                            double dx0 = (xx+0.5) - sx;
                            double quad = A*dx0*dx0 + 2*B*dx0*dy0 + C*dy0*dy0;
                            if (quad > kRad*kRad*2.0) continue;
                            double w = Math.exp(-0.5*quad);

                            int px=xx, py=yy;
                            if (clampEdge) {
                                px = clamp(px,0,sw-1);
                                py = clamp(py,0,sh-1);
                            } else {
                                if (px<0||py<0||px>=sw||py>=sh) continue;
                            }
                            int c0 = sdata[py*sw + px];
                            double a0 = ((c0>>>24)&255)/255.0;
                            double r0 = ((c0>>>16)&255)/255.0;
                            double g0 = ((c0>>>8 )&255)/255.0;
                            double b0 = ( c0      &255)/255.0;

                            W  += w;
                            Ra += w*a0;
                            Rr += w*r0*a0;
                            Rg += w*g0*a0;
                            Rb += w*b0*a0;
                        }
                    }
                    int outARGB;
                    if (W<=1e-12 || Ra<=1e-12) {
                        outARGB = 0x00000000;
                    } else {
                        double A1 = Ra/W;
                        double R1 = (Rr/W)/(A1+1e-12);
                        double G1 = (Rg/W)/(A1+1e-12);
                        double B1 = (Rb/W)/(A1+1e-12);
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

        // Catmull–Rom bicubic sampler (premultiplied alpha)
        static int sampleBicubicARGB(int[] data, int w, int h, double x, double y, boolean clamp){
            if (!clamp && (x < -0.5 || y < -0.5 || x > w-0.5 || y > h-0.5)) return 0x00000000;
            int x1 = (int)Math.floor(x), y1 = (int)Math.floor(y);
            double tx = x - x1, ty = y - y1;

            int[] xs = new int[4], ys = new int[4];
            for (int i=-1;i<=2;i++) xs[i+1] = clamp? clamp(x1+i,0,w-1) : (x1+i);
            for (int j=-1;j<=2;j++) ys[j+1] = clamp? clamp(y1+j,0,h-1) : (y1+j);

            double[] wx = catmullRomWeights(tx);
            double[] wy = catmullRomWeights(ty);

            double W=0, Aa=0, Ar=0, Ag=0, Ab=0;
            for (int j=0;j<4;j++){
                int yy = ys[j]; if (!clamp && (yy<0||yy>=h)) continue;
                for (int i=0;i<4;i++){
                    int xx = xs[i]; if (!clamp && (xx<0||xx>=w)) continue;
                    double wgt = wx[i]*wy[j];
                    int c = data[yy*w + xx];
                    double a = ((c>>>24)&255)/255.0;
                    double r = ((c>>>16)&255)/255.0;
                    double g = ((c>>>8 )&255)/255.0;
                    double b = ( c      &255)/255.0;
                    W += wgt; Aa += wgt*a; Ar += wgt*r*a; Ag += wgt*g*a; Ab += wgt*b*a;
                }
            }
            if (W<=1e-12 || Aa<=1e-12) return 0x00000000;
            double A1 = Aa/W;
            double R1 = (Ar/W)/(A1+1e-12);
            double G1 = (Ag/W)/(A1+1e-12);
            double B1 = (Ab/W)/(A1+1e-12);
            int ia=clamp((int)Math.round(A1*255),0,255);
            int ir=clamp((int)Math.round(R1*255),0,255);
            int ig=clamp((int)Math.round(G1*255),0,255);
            int ib=clamp((int)Math.round(B1*255),0,255);
            return (ia<<24)|(ir<<16)|(ig<<8)|ib;
        }

        static double[] catmullRomWeights(double t){
            double a=-0.5, t2=t*t, t3=t2*t;
            double w0 = a*(-t3 + 2*t2 - t);
            double w1 = (a+2)*t3 + (-a-3)*t2 + 1;
            double w2 = (a+2)*(-t3) + (2*a+3)*t2 + (-a)*t;
            double w3 = a*(t3 - t2);
            return new double[]{w0,w1,w2,w3};
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
