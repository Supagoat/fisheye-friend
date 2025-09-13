package defish;


/* defaults assume your two sample files and the given seed
java AutoHemiPointMapper gridColor.png gridColor_Hemi_Fullframe.png hemi_pairs_auto.json

# useful options:
#   --seedSrc=399,100     (override source seed)
#   --seedDst=591,140     (override target seed)
#   --win=21              (patch size is (2*win+1))
#   --search=32           (search radius in target, pixels)
#   --minDx=30 --maxDx=90 --minDy=30 --maxDy=90  (period search ranges, source & target)
#   --minDxT=30 --maxDxT=90 --minDyT=30 --maxDyT=90  (period search ranges, target)
#   --stride=1            (grid subsampling; use 1 for every grid node, 2 to skip every other)
#   --maxNodes=4000       (cap)
*/

//AutoHemiPointMapper.java
//Build point mappings between two grid images by robustly matching grid intersections.
//Emphasizes horizontal/vertical edges and down-weights diagonals (watermark).
//
//Usage:
//javac AutoHemiPointMapper.java
//java AutoHemiPointMapper <source.png> <target.png> <out.json> [flags...]
//
//Flags (all optional):
//--seedSrc=xs,ys       default 399,100
//--seedDst=xd,yd       default 591,140
//--win=N               half window size (default 15 -> 31x31)
//--search=R            search radius in target (default 32)
//--minDx=30 --maxDx=90 --minDy=30 --maxDy=90     (source grid period search)
//--minDxT=30 --maxDxT=90 --minDyT=30 --maxDyT=90 (target grid period search)
//--stride=S            sample every S-th grid intersection (default 1)
//--maxNodes=M          safety cap (default 4000)

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class AutoHemiPointMapper {

 public static void main(String[] args) throws Exception {
     if (args.length < 3) {
         System.out.println("Usage: java AutoHemiPointMapper <source.png> <target.png> <out.json> [flags...]");
         return;
     }
     String srcPath = args[0];
     String dstPath = args[1];
     String outPath = args[2];

     Params P = new Params();
     for (int i = 3; i < args.length; i++) P.parseFlag(args[i]);

     // Read images and convert to grayscale + ARGB buffer
     BufferedImage srcImgIn = ImageIO.read(new File(srcPath));
     BufferedImage dstImgIn = ImageIO.read(new File(dstPath));
     if (srcImgIn == null || dstImgIn == null) throw new IOException("Cannot read input images.");

     BufferedImage srcImg = toIntARGB(srcImgIn);
     BufferedImage dstImg = toIntARGB(dstImgIn);

     Gray gSrc = Gray.fromARGB(srcImg);
     Gray gDst = Gray.fromARGB(dstImg);

     // Grid-emphasis maps (down-weight diagonals to ignore watermark)
     Gray wSrc = gridEmphasis(gSrc);
     Gray wDst = gridEmphasis(gDst);

     // Estimate grid periods in source and target, near seed
     double dxS = estimatePeriodX(wSrc, P.seedSrcX, P.seedSrcY, P.minDx, P.maxDx, P.win);
     double dyS = estimatePeriodY(wSrc, P.seedSrcX, P.seedSrcY, P.minDy, P.maxDy, P.win);
     double dxT = estimatePeriodX(wDst, P.seedDstX, P.seedDstY, P.minDxT, P.maxDxT, P.win);
     double dyT = estimatePeriodY(wDst, P.seedDstX, P.seedDstY, P.minDyT, P.maxDyT, P.win);

     System.out.printf(Locale.US,
             "Estimated periods  Source: dx=%.2f dy=%.2f   Target: dx=%.2f dy=%.2f%n",
             dxS, dyS, dxT, dyT);

     // Generate source grid (lattice of intersections)
     ArrayList<PointI> srcNodes = buildSourceGrid(srcImg.getWidth(), srcImg.getHeight(),
                                                  P.seedSrcX, P.seedSrcY, dxS, dyS, P.stride, P.maxNodes);

     // For each source node, predict target location and refine with weighted NCC
     ArrayList<Pair> pairs = new ArrayList<>();
     pairs.add(new Pair(0, P.seedSrcX, P.seedSrcY, P.seedDstX, P.seedDstY)); // include seed

     // Build a template extractor around each source node (from wSrc map)
     for (int idx = 0; idx < srcNodes.size(); idx++) {
         PointI s = srcNodes.get(idx);
         // skip the seed itself (already added)
         if (s.x == P.seedSrcX && s.y == P.seedSrcY) continue;

         // Predict in target using period offsets from seed
         int di = (int)Math.round((s.x - P.seedSrcX) / dxS);
         int dj = (int)Math.round((s.y - P.seedSrcY) / dyS);
         double xPred = P.seedDstX + di * dxT;
         double yPred = P.seedDstY + dj * dyT;

         // Refine with weighted NCC around predicted location
         MatchResult mr = matchWNCC(wSrc, wDst, s.x, s.y, xPred, yPred, P.win, P.search);

         if (mr != null && mr.score > 0.35) { // conservative threshold
             pairs.add(new Pair(pairs.size(), s.x, s.y, mr.x, mr.y));
         }
     }

     // Write JSON
     writeJson(outPath, srcImg.getWidth(), srcImg.getHeight(), dstImg.getWidth(), dstImg.getHeight(), pairs);
     System.out.println("Wrote " + outPath + " with " + pairs.size() + " pairs.");
 }

 // ---------- Parameters ----------
 static class Params {
     int seedSrcX = 399, seedSrcY = 100;
     int seedDstX = 591, seedDstY = 140;
     int win = 15;      // half window -> 31x31
     int search = 32;   // search radius in target
     int minDx = 30, maxDx = 90, minDy = 30, maxDy = 90;
     int minDxT = 30, maxDxT = 90, minDyT = 30, maxDyT = 90;
     int stride = 1;
     int maxNodes = 4000;

     void parseFlag(String f) {
         if (f.startsWith("--seedSrc=")) {
             String[] p = f.substring(10).split(",");
             seedSrcX = Integer.parseInt(p[0]); seedSrcY = Integer.parseInt(p[1]);
         } else if (f.startsWith("--seedDst=")) {
             String[] p = f.substring(10).split(",");
             seedDstX = Integer.parseInt(p[0]); seedDstY = Integer.parseInt(p[1]);
         } else if (f.startsWith("--win=")) win = Integer.parseInt(f.substring(6));
         else if (f.startsWith("--search=")) search = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--minDx=")) minDx = Integer.parseInt(f.substring(8));
         else if (f.startsWith("--maxDx=")) maxDx = Integer.parseInt(f.substring(8));
         else if (f.startsWith("--minDy=")) minDy = Integer.parseInt(f.substring(8));
         else if (f.startsWith("--maxDy=")) maxDy = Integer.parseInt(f.substring(8));
         else if (f.startsWith("--minDxT=")) minDxT = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--maxDxT=")) maxDxT = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--minDyT=")) minDyT = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--maxDyT=")) maxDyT = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--stride=")) stride = Integer.parseInt(f.substring(9));
         else if (f.startsWith("--maxNodes=")) maxNodes = Integer.parseInt(f.substring(11));
     }
 }

 // ---------- Core data ----------
 static class Pair {
     int id; int sx, sy; double dx, dy;
     Pair(int id, int sx, int sy, double dx, double dy){ this.id=id; this.sx=sx; this.sy=sy; this.dx=dx; this.dy=dy; }
 }
 static class PointI { int x,y; PointI(int x,int y){this.x=x;this.y=y;} }
 static class MatchResult { double x,y,score; MatchResult(double x,double y,double s){this.x=x;this.y=y;this.score=s;} }

 // ---------- Gray image & gradients ----------
 static class Gray {
     final int w,h;
     final double[] v;  // grayscale
     Gray(int w,int h,double[] v){this.w=w;this.h=h;this.v=v;}
     static Gray fromARGB(BufferedImage img){
         int w = img.getWidth(), h = img.getHeight();
         int[] data = ((DataBufferInt)img.getRaster().getDataBuffer()).getData();
         double[] g = new double[w*h];
         for (int i=0;i<w*h;i++){
             int c = data[i];
             int r=(c>>>16)&255, gg=(c>>>8)&255, b=c&255;
             g[i] = 0.2989*r + 0.5870*gg + 0.1140*b;
         }
         return new Gray(w,h,g);
     }
     double get(int x,int y){
         if (x<0||y<0||x>=w||y>=h) return 0;
         return v[y*w+x];
     }
 }

 /** Build “grid emphasis” map that favors vertical/horizontal edges and penalizes diagonals (watermark). */
 static Gray gridEmphasis(Gray g){
     int w=g.w,h=g.h;
     double[] out = new double[w*h];
     // Sobel
     int[][] kx = {{-1,0,1},{-2,0,2},{-1,0,1}};
     int[][] ky = {{-1,-2,-1},{0,0,0},{1,2,1}};
     for (int y=1;y<h-1;y++){
         for (int x=1;x<w-1;x++){
             double gx=0, gy=0;
             for (int j=-1;j<=1;j++) for(int i=-1;i<=1;i++){
                 double p=g.get(x+i,y+j);
                 gx += kx[j+1][i+1]*p;
                 gy += ky[j+1][i+1]*p;
             }
             // axis preference (|gx| and |gy|) vs. diagonal penalty (at 45°)
             double ax = Math.abs(gx), ay = Math.abs(gy);
             double diag1 = Math.abs((gx+gy)*0.70710678);
             double diag2 = Math.abs((gx-gy)*0.70710678);
             double axisScore = Math.pow(ax, 2.0) + Math.pow(ay, 2.0);
             double diagPenalty = 0.6*(Math.pow(diag1, 2.0) + Math.pow(diag2, 2.0));
             double s = axisScore - diagPenalty;
             if (s < 0) s = 0;
             out[y*w+x] = s;
         }
     }
     // Normalize
     double max=0; for(double v:out) if (v>max) max=v;
     if (max>0) for(int i=0;i<out.length;i++) out[i]/=max;
     return new Gray(w,h,out);
 }

 // ---------- Period estimation near a point via 1D matching ----------
 static double estimatePeriodX(Gray w, int x0, int y0, int min, int max, int win){
     double bestD=min, bestS=-1;
     double[] tpl = extractRow(w, x0, y0, win);
     for (int d=min; d<=max; d++){
         double s = rowCorr(tpl, extractRow(w, x0+d, y0, win));
         if (s > bestS) { bestS=s; bestD=d; }
     }
     return bestD;
 }
 static double estimatePeriodY(Gray w, int x0, int y0, int min, int max, int win){
     double bestD=min, bestS=-1;
     double[] tpl = extractCol(w, x0, y0, win);
     for (int d=min; d<=max; d++){
         double s = colCorr(tpl, extractCol(w, x0, y0+d, win));
         if (s > bestS) { bestS=s; bestD=d; }
     }
     return bestD;
 }
 static double[] extractRow(Gray w,int cx,int cy,int win){
     int L=2*win+1;
     double[] a=new double[L];
     for(int i=-win;i<=win;i++) a[i+win]=w.get(cx+i,cy);
     return a;
 }
 static double[] extractCol(Gray w,int cx,int cy,int win){
     int L=2*win+1;
     double[] a=new double[L];
     for(int j=-win;j<=win;j++) a[j+win]=w.get(cx,cy+j);
     return a;
 }
 static double rowCorr(double[] a,double[] b){
     double ma=mean(a), mb=mean(b);
     double sa=0,sb=0,ab=0;
     for(int i=0;i<a.length;i++){
         double da=a[i]-ma, db=b[i]-mb;
         sa+=da*da; sb+=db*db; ab+=da*db;
     }
     double denom = Math.sqrt(sa*sb)+1e-9;
     return ab/denom;
 }
 static double colCorr(double[] a,double[] b){ return rowCorr(a,b); }
 static double mean(double[] a){ double s=0; for(double v:a)s+=v; return s/a.length; }

 // ---------- Build the source lattice ----------
 static ArrayList<PointI> buildSourceGrid(int w,int h,int x0,int y0,double dx,double dy,int stride,int cap){
     ArrayList<PointI> pts = new ArrayList<>();
     pts.add(new PointI(x0,y0));
     // horizontal to the right
     for (int i=stride; ; i+=stride){
         int x = (int)Math.round(x0 + i*dx);
         if (x >= w-1) break;
         pts.add(new PointI(x,y0));
         if (pts.size() >= cap) break;
     }
     // horizontal to the left
     for (int i=-stride; ; i-=stride){
         int x = (int)Math.round(x0 + i*dx);
         if (x <= 0) break;
         pts.add(new PointI(x,y0));
         if (pts.size() >= cap) break;
     }
     // create rows by stepping vertically
     ArrayList<PointI> all = new ArrayList<>();
     for (PointI base : pts) {
         // up
         for (int j=stride;; j+=stride){
             int y = (int)Math.round(y0 + j*dy);
             if (y >= h-1) break;
             all.add(new PointI(base.x, y));
             if (all.size() >= cap) break;
         }
         if (all.size() >= cap) break;
         // down
         for (int j=-stride;; j-=stride){
             int y = (int)Math.round(y0 + j*dy);
             if (y <= 0) break;
             all.add(new PointI(base.x, y));
             if (all.size() >= cap) break;
         }
         all.add(base); // include the base row points themselves
     }
     // deduplicate
     HashSet<Long> seen = new HashSet<>();
     ArrayList<PointI> uniq = new ArrayList<>();
     for (PointI p : all){
         long key = (((long)p.x)<<32) ^ (p.y & 0xffffffffL);
         if (seen.add(key)) uniq.add(p);
         if (uniq.size() >= cap) break;
     }
     // sort for consistent IDs
     uniq.sort(Comparator.comparingInt((PointI p)->p.y).thenComparingInt(p->p.x));
     return uniq;
 }

 // ---------- Weighted NCC matcher (watermark-robust) ----------
 static MatchResult matchWNCC(Gray wSrc, Gray wDst, int sx, int sy, double xPred, double yPred, int win, int search){
     // extract source patch & weights
     Patch src = Patch.extract(wSrc, sx, sy, win);
     if (src == null) return null;

     double best = -1; double bx=xPred, by=yPred;
     int xmin=(int)Math.round(xPred)-search, xmax=(int)Math.round(xPred)+search;
     int ymin=(int)Math.round(yPred)-search, ymax=(int)Math.round(yPred)+search;

     for (int y=ymin; y<=ymax; y++){
         for (int x=xmin; x<=xmax; x++){
             Patch tgt = Patch.extract(wDst, x, y, win);
             if (tgt == null) continue;
             double s = weightedNCC(src, tgt);
             if (s > best){ best = s; bx = x; by = y; }
         }
     }
     return new MatchResult(bx, by, best);
 }

 static class Patch {
     final int w,h, win;   // size = 2*win+1
     final double[] a;     // intensities (from emphasis map)
     final double[] wt;    // weights
     Patch(int w,int h,int win,double[] a,double[] wt){this.w=w;this.h=h;this.win=win;this.a=a;this.wt=wt;}
     static Patch extract(Gray g,int cx,int cy,int win){
         int W=2*win+1;
         if (cx-win<0||cy-win<0||cx+win>=g.w||cy+win>=g.h) return null;
         double[] a=new double[W*W];
         double[] w=new double[W*W];
         int k=0;
         for (int j=-win;j<=win;j++){
             for (int i=-win;i<=win;i++){
                 double v=g.get(cx+i,cy+j);
                 a[k]=v;
                 // weight uses emphasis value itself; amplify center a bit
                 double rad2 = (i*i + j*j) / (double)(win*win+1);
                 double centerBoost = 1.0 - 0.3*rad2;
                 w[k]=Math.max(0.0, v)*centerBoost;
                 k++;
             }
         }
         return new Patch(W,W,win,a,w);
     }
 }

 static double weightedNCC(Patch A, Patch B){
     // combine weights to suppress features present only in one patch (watermark)
     double num=0, denA=0, denB=0;
     double sumW=0, mA=0, mB=0;
     for (int i=0;i<A.a.length;i++){
         double w = Math.sqrt(A.wt[i]*B.wt[i]) + 1e-9;
         sumW += w;
         mA += w*A.a[i];
         mB += w*B.a[i];
     }
     mA/=sumW; mB/=sumW;
     for (int i=0;i<A.a.length;i++){
         double w = Math.sqrt(A.wt[i]*B.wt[i]) + 1e-9;
         double da = A.a[i]-mA;
         double db = B.a[i]-mB;
         num += w*da*db;
         denA += w*da*da;
         denB += w*db*db;
     }
     double denom = Math.sqrt(denA*denB) + 1e-9;
     return num/denom;
 }

 // ---------- JSON writer ----------
 static void writeJson(String outPath, int sW, int sH, int tW, int tH, ArrayList<Pair> pairs) throws IOException {
     try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(outPath), StandardCharsets.UTF_8))) {
         pw.println("{");
         pw.printf("  \"source\": {\"width\": %d, \"height\": %d},%n", sW, sH);
         pw.printf("  \"target\": {\"width\": %d, \"height\": %d},%n", tW, tH);
         pw.println("  \"pairs\": [");
         for (int i=0;i<pairs.size();i++){
             Pair p = pairs.get(i);
             pw.printf(Locale.US,
                 "    {\"id\": %d, \"src\": {\"x\": %.3f, \"y\": %.3f}, \"dst\": {\"x\": %.3f, \"y\": %.3f}}%s%n",
                 p.id, (double)p.sx, (double)p.sy, p.dx, p.dy, (i==pairs.size()-1?"":","));
         }
         pw.println("  ]");
         pw.println("}");
     }
 }

 // ---------- Utilities ----------
 static BufferedImage toIntARGB(BufferedImage in){
     if (in.getType() == BufferedImage.TYPE_INT_ARGB) return in;
     BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
     Graphics2D g = out.createGraphics();
     g.setComposite(AlphaComposite.Src);
     g.drawImage(in, 0, 0, null);
     g.dispose();
     return out;
 }
}
