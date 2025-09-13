package defish;

import com.google.gson.*;
import com.google.gson.annotations.SerializedName;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Apply a Fisheye-Hemi style warp using a Thin-Plate Spline fit
 * from destination (hemi output) -> source (fisheye input), then
 * render by inverse mapping (for each output pixel, sample input).
 *
 * Usage:
 *   java -cp .;gson-2.11.0.jar HemiWarp in.jpg out.png hemi_pairs_321.json [--lambda 1e-6] [--grid 8] [--zoom 1.00]
 *
 * Notes
 *  - Control points are consumed in normalized coords (0..1), so this generalizes to other sizes.
 *  - A coarse grid (default 8 px) is used to evaluate the TPS once per grid node, then bilinearly
 *    upsampled to every pixel for speed (hundreds of times faster than per-pixel TPS evaluation).
 *  - Regularization (--lambda) stabilizes the TPS solution.
 */
public class MapDefishV9 {

    // -------------- Data models for the uploaded JSON -----------------
    static class ImageMeta {
        String path;
        int width;
        int height;
    }
    static class Vec {
        double x, y;
        Double xn, yn; // normalized (0..1) if present
    }
    static class Pair {
        int id;
        Vec src; // fisheye input (original)
        Vec dst; // hemi output
    }
    static class PairsFile {
        ImageMeta source;
        ImageMeta target;
        List<Pair> pairs;
    }

    // -------------- Thin-Plate Spline (2D -> 1D) ---------------------
    static class TPS {
        // TPS maps (x,y) -> value, trained on control points (xi, yi) -> vi
        private final double[] xi, yi, vi;
        private final int n;
        private final double lambda;
        // Solution coefficients: w_i (size n) and affine a0 + a1*x + a2*y
        private double[] w;
        private double a0, a1, a2;

        TPS(double[] xi, double[] yi, double[] vi, double lambda) {
            this.xi = xi; this.yi = yi; this.vi = vi;
            this.n = xi.length;
            this.lambda = lambda;
            solve();
        }

        private static double U(double r2) {
            if (r2 <= 1e-300) return 0.0;
            return r2 * Math.log(r2);
        }

        private void solve() {
            // Build L matrix of size (n+3) x (n+3)
            final int m = n + 3;
            double[][] L = new double[m][m];
            // Top-left K + thetaI
            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    double dx = xi[i] - xi[j];
                    double dy = yi[i] - yi[j];
                    double r2 = dx*dx + dy*dy;
                    double val = U(r2);
                    if (i == j) val += lambda;
                    L[i][j] = val;
                    L[j][i] = val;
                }
            }
            // Top-right P
            for (int i = 0; i < n; i++) {
                L[i][n + 0] = 1.0;
                L[i][n + 1] = xi[i];
                L[i][n + 2] = yi[i];
            }
            // Bottom-left P^T
            for (int j = 0; j < n; j++) {
                L[n + 0][j] = 1.0;
                L[n + 1][j] = xi[j];
                L[n + 2][j] = yi[j];
            }
            // Bottom-right zeros already

            // Right-hand side: [v; 0; 0; 0]
            double[] b = new double[m];
            System.arraycopy(vi, 0, b, 0, n);

            // Solve L * coeff = b
            double[] coeff = solveSymmetricSystem(L, b); // partial pivoting

            // Extract
            w = new double[n];
            System.arraycopy(coeff, 0, w, 0, n);
            a0 = coeff[n + 0];
            a1 = coeff[n + 1];
            a2 = coeff[n + 2];
        }

        public double evaluate(double x, double y) {
            double sum = a0 + a1 * x + a2 * y;
            for (int i = 0; i < n; i++) {
                double dx = x - xi[i];
                double dy = y - yi[i];
                sum += w[i] * U(dx*dx + dy*dy);
            }
            return sum;
        }

        // Simple Gaussian elimination with partial pivoting (works for our sizes)
        private static double[] solveSymmetricSystem(double[][] A, double[] b) {
            int n = b.length;
            for (int p = 0; p < n; p++) {
                // Pivot
                int max = p;
                double best = Math.abs(A[p][p]);
                for (int i = p + 1; i < n; i++) {
                    double v = Math.abs(A[i][p]);
                    if (v > best) { best = v; max = i; }
                }
                if (best < 1e-15) {
                    throw new RuntimeException("Ill-conditioned TPS system; try larger --lambda");
                }
                // Swap rows
                if (max != p) {
                    double[] tmp = A[p]; A[p] = A[max]; A[max] = tmp;
                    double tb = b[p]; b[p] = b[max]; b[max] = tb;
                }
                // Eliminate
                double diag = A[p][p];
                for (int j = p; j < n; j++) A[p][j] /= diag;
                b[p] /= diag;

                for (int i = 0; i < n; i++) {
                    if (i == p) continue;
                    double factor = A[i][p];
                    if (factor == 0) continue;
                    for (int j = p; j < n; j++) A[i][j] -= factor * A[p][j];
                    b[i] -= factor * b[p];
                }
            }
            return b;
        }
    }

    // -------------- Utility: bilinear sampling/clamp -----------------
    static int sampleBilinear(BufferedImage img, double x, double y) {
        int w = img.getWidth(), h = img.getHeight();
        // clamp to [0, w-1], [0, h-1]
        if (x < 0) x = 0; if (x > w - 1) x = w - 1;
        if (y < 0) y = 0; if (y > h - 1) y = h - 1;
        int x0 = (int) Math.floor(x), y0 = (int) Math.floor(y);
        int x1 = Math.min(x0 + 1, w - 1), y1 = Math.min(y0 + 1, h - 1);
        double dx = x - x0, dy = y - y0;

        int c00 = img.getRGB(x0, y0);
        int c10 = img.getRGB(x1, y0);
        int c01 = img.getRGB(x0, y1);
        int c11 = img.getRGB(x1, y1);

        int r = (int) lerp(lerp((c00 >> 16) & 255, (c10 >> 16) & 255, dx),
                           lerp((c01 >> 16) & 255, (c11 >> 16) & 255, dx), dy);
        int g = (int) lerp(lerp((c00 >> 8) & 255, (c10 >> 8) & 255, dx),
                           lerp((c01 >> 8) & 255, (c11 >> 8) & 255, dx), dy);
        int b = (int) lerp(lerp(c00 & 255, c10 & 255, dx),
                           lerp(c01 & 255, c11 & 255, dx), dy);
        return (0xFF << 24) | (r << 16) | (g << 8) | b;
    }
    static double lerp(double a, double b, double t) { return a + (b - a) * t; }

    // -------------- Main ---------------------------------------------
    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.out.println("Usage: java -cp .;gson.jar HemiWarp <inputImage> <outputImage> <pairs.json> [--lambda 1e-6] [--grid 8] [--zoom 1.00]");
            return;
        }
        String inPath = args[0];
        String outPath = args[1];
        String jsonPath = args[2];

        double lambda = 1e-6;   // TPS regularization
        int grid = 8;           // grid step in pixels (compute TPS at grid nodes)
        double zoom = 1.0;      // simple pre-zoom around center

        for (int i = 3; i < args.length; i++) {
            if (args[i].equals("--lambda") && i + 1 < args.length) lambda = Double.parseDouble(args[++i]);
            else if (args[i].equals("--grid") && i + 1 < args.length) grid = Integer.parseInt(args[++i]);
            else if (args[i].equals("--zoom") && i + 1 < args.length) zoom = Double.parseDouble(args[++i]);
        }

        // Load image
        BufferedImage src = ImageIO.read(new File(inPath));
        if (src == null) throw new IOException("Could not read input image: " + inPath);

        int W = src.getWidth();
        int H = src.getHeight();

        // Load pairs JSON (dst -> src), consume normalized coords
        PairsFile pf = loadPairs(jsonPath);
        List<Pair> list = pf.pairs;

        // Build training arrays in normalized space
        int n = list.size();
        double[] dx = new double[n], dy = new double[n]; // dst normalized
        double[] sx = new double[n], sy = new double[n]; // src normalized
        for (int i = 0; i < n; i++) {
            Pair p = list.get(i);
            double tx = (p.dst.xn != null) ? p.dst.xn : (p.dst.x / (double) pf.target.width);
            double ty = (p.dst.yn != null) ? p.dst.yn : (p.dst.y / (double) pf.target.height);
            double qx = (p.src.xn != null) ? p.src.xn : (p.src.x / (double) pf.source.width);
            double qy = (p.src.yn != null) ? p.src.yn : (p.src.y / (double) pf.source.height);
            dx[i] = clamp01(tx);
            dy[i] = clamp01(ty);
            sx[i] = clamp01(qx);
            sy[i] = clamp01(qy);
        }

        // Fit two TPS maps: (dst.x, dst.y) -> src.x  and  -> src.y
        TPS tpsX = new TPS(dx, dy, sx, lambda);
        TPS tpsY = new TPS(dx, dy, sy, lambda);

        // Precompute forward map on a coarse grid in *output* space (normalized)
        int gw = Math.max(2, (W + grid - 1) / grid + 1);
        int gh = Math.max(2, (H + grid - 1) / grid + 1);
        double[][] mapX = new double[gh][gw];
        double[][] mapY = new double[gh][gw];

        for (int gy = 0; gy < gh; gy++) {
            for (int gx = 0; gx < gw; gx++) {
                double u = ((gx * grid) / (double) (W - 1));   // 0..1
                double v = ((gy * grid) / (double) (H - 1));   // 0..1

                // optional center zoom before TPS evaluation
                if (zoom != 1.0) {
                    u = (u - 0.5) / zoom + 0.5;
                    v = (v - 0.5) / zoom + 0.5;
                }

                double sxn = tpsX.evaluate(u, v);
                double syn = tpsY.evaluate(u, v);
                mapX[gy][gx] = sxn;
                mapY[gy][gx] = syn;
            }
        }

        // Render
        BufferedImage out = new BufferedImage(W, H, BufferedImage.TYPE_INT_ARGB);
        for (int y = 0; y < H; y++) {
            double v = y / (double) (H - 1);
            // grid cell vertical indices
            int gy0 = Math.min((y / grid), gh - 2);
            int gy1 = gy0 + 1;
            double tv = ((y - gy0 * grid) / (double) Math.max(1, grid));

            for (int x = 0; x < W; x++) {
                double u = x / (double) (W - 1);
                if (zoom != 1.0) {
                    u = (u - 0.5) / zoom + 0.5;
                    v = (v - 0.5) / zoom + 0.5;
                }
                // grid cell horizontal indices
                int gx0 = Math.min((x / grid), gw - 2);
                int gx1 = gx0 + 1;
                double tu = ((x - gx0 * grid) / (double) Math.max(1, grid));

                // Bilinear in the precomputed map
                double s00x = mapX[gy0][gx0], s10x = mapX[gy0][gx1];
                double s01x = mapX[gy1][gx0], s11x = mapX[gy1][gx1];
                double s00y = mapY[gy0][gx0], s10y = mapY[gy0][gx1];
                double s01y = mapY[gy1][gx0], s11y = mapY[gy1][gx1];

                double sxn = lerp(lerp(s00x, s10x, tu), lerp(s01x, s11x, tu), tv);
                double syn = lerp(lerp(s00y, s10y, tu), lerp(s01y, s11y, tu), tv);

                // to source pixel coords
                double sxPix = clamp((W - 1) * sxn, 0, W - 1);
                double syPix = clamp((H - 1) * syn, 0, H - 1);

                int rgba = sampleBilinear(src, sxPix, syPix);
                out.setRGB(x, y, rgba);
            }
        }

        // Write
        String fmt = guessFormat(outPath);
        ImageIO.write(out, fmt, new File(outPath));
        System.out.println("Wrote: " + outPath);
    }

    private static String guessFormat(String p) {
        int i = p.lastIndexOf('.');
        String ext = (i >= 0) ? p.substring(i + 1).toLowerCase(Locale.ROOT) : "png";
        if (ext.equals("jpg")) return "jpg";
        if (ext.equals("jpeg")) return "jpg";
        if (ext.equals("bmp")) return "bmp";
        if (ext.equals("gif")) return "gif";
        return "png";
    }

    private static double clamp01(double v) { return v < 0 ? 0 : (v > 1 ? 1 : v); }
    private static double clamp(double v, double a, double b) { return Math.max(a, Math.min(b, v)); }

    // ---------------- JSON loader ----------------
    private static PairsFile loadPairs(String path) throws IOException {
        String json = Files.readString(Paths.get(path));
        Gson gson = new GsonBuilder().create();
        return gson.fromJson(json, PairsFile.class);
    }
}
