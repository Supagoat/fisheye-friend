package defish;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.Color;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * HemiClone — reverse engineer a Fisheye Hemi style warp using landmark mappings.
 *
 * What it does
 * ============
 * 1) Reads an input image.
 * 2) Loads landmark pairs (src->dst) from a JSON mapping file in the format you supplied
 *    (this program extracts the normalized fields `src.xn, src.yn, dst.xn, dst.yn`).
 *    It also supports a simple CSV with columns: src_xn,src_yn,dst_xn,dst_yn.
 * 3) Fits two Thin Plate Splines (TPS):
 *       forward  (src -> dst)
 *       inverse  (dst -> src)   <-- used for resampling so every output pixel pulls a source sample
 * 4) Builds a coarse grid of inverse mappings for speed, then bilinearly interpolates per pixel.
 * 5) Writes the result as PNG.
 *
 * Why TPS?
 * ========
 * TPS is a smooth, well behaved warp that honors all provided correspondences while
 * regularization (lambda) damps oscillations ("wavy corners"). With your landmarks this
 * closely tracks the real Fisheye Hemi transform without hard coding a projection model.
 *
 * Usage
 * =====
 *   javac HemiClone.java
 *   java HemiClone \
 *     --input in.jpg \
 *     --output out.png \
 *     --pairs hemi_pairs_384.json \
 *     [--lambda 1e-3] [--grid 96x64] [--outW 0] [--outH 0] [--background edge|black|transparent]
 *
 * Notes
 * =====
 * - Coordinates in the mapping file are interpreted as normalized [0,1] across width/height.
 * - If outW/outH are 0, the output keeps the input size.
 * - Grid controls the acceleration quality/speed tradeoff (more cells => crisper, slower precompute).
 *
 * Public domain / MIT style. No external libraries required.
 */
public class MapDefishV7 {

    public static void main(String[] args) throws Exception {
        Map<String,String> cli = parseArgs(args);
        if (!cli.containsKey("--input") || !cli.containsKey("--output") || !cli.containsKey("--pairs")) {
            printHelp();
            return;
        }

        String inPath = cli.get("--input");
        String outPath = cli.get("--output");
        String pairsPath = cli.get("--pairs");
        double lambda = cli.containsKey("--lambda") ? Double.parseDouble(cli.get("--lambda")) : 1e-3;
        String gridStr = cli.getOrDefault("--grid", "96x64");
        int gridX = 96, gridY = 64;
        try {
            String[] g = gridStr.toLowerCase(Locale.ROOT).split("x");
            gridX = Math.max(8, Integer.parseInt(g[0]));
            gridY = Math.max(8, Integer.parseInt(g[1]));
        } catch (Exception ignored) {}
        int outW = cli.containsKey("--outW") ? Integer.parseInt(cli.get("--outW")) : 0;
        int outH = cli.containsKey("--outH") ? Integer.parseInt(cli.get("--outH")) : 0;
        String bg = cli.getOrDefault("--background", "edge").toLowerCase(Locale.ROOT);

        BufferedImage src = ImageIO.read(new File(inPath));
        if (src == null) throw new IOException("Failed to read input image: " + inPath);
        if (outW <= 0) outW = src.getWidth();
        if (outH <= 0) outH = src.getHeight();

        List<Pair> pairs = loadPairsFlexible(pairsPath);
        if (pairs.isEmpty()) throw new IllegalArgumentException("No landmark pairs parsed from: " + pairsPath);

        // Build TPS models
        ThinPlateSpline2D tpsSrcToDst = ThinPlateSpline2D.fit(
                pairs.stream().map(p -> new double[]{p.sx, p.sy}).toArray(double[][]::new),
                pairs.stream().mapToDouble(p -> p.dx).toArray(),
                pairs.stream().mapToDouble(p -> p.dy).toArray(),
                lambda
        );
        ThinPlateSpline2D tpsDstToSrc = ThinPlateSpline2D.fit(
                pairs.stream().map(p -> new double[]{p.dx, p.dy}).toArray(double[][]::new),
                pairs.stream().mapToDouble(p -> p.sx).toArray(),
                pairs.stream().mapToDouble(p -> p.sy).toArray(),
                lambda
        );

        // Build inverse grid (dst->src) for speed
        InverseGrid grid = InverseGrid.build(tpsDstToSrc, gridX, gridY);

        // Warp
        BufferedImage out = new BufferedImage(outW, outH, bg.equals("transparent") ? BufferedImage.TYPE_INT_ARGB : BufferedImage.TYPE_INT_RGB);
        Sampler sampler = new Sampler(src, bg);

        for (int y = 0; y < outH; y++) {
            double yn = (outH == 1) ? 0.0 : (y / (double)(outH - 1));
            for (int x = 0; x < outW; x++) {
                double xn = (outW == 1) ? 0.0 : (x / (double)(outW - 1));

                // approximate inverse mapping via grid bilinear interpolation
                double[] srcN = grid.map(xn, yn);

                // convert to source pixel space
                double sx = srcN[0] * (src.getWidth() - 1);
                double sy = srcN[1] * (src.getHeight() - 1);

                int rgba = sampler.sampleBilinearRGBA(sx, sy);
                out.setRGB(x, y, rgba);
            }
        }

        ImageIO.write(out, "png", new File(outPath));
        System.out.println("Wrote: " + outPath);
    }

    private static void printHelp() {
        System.out.println("\nHemiClone — TPS warp from landmark pairs (src<->dst)\n" +
                "\nRequired:\n  --input <file>     input image (jpg/png/etc)\n  --output <file>    output PNG path\n  --pairs <file>     mapping file (JSON with src.xn/src.yn/dst.xn/dst.yn, or CSV)\n" +
                "\nOptional:\n  --lambda <val>     regularization (default 1e-3). Higher => smoother, fewer ripples\n  --grid <WxH>       grid resolution for inverse map (default 96x64)\n  --outW <int>       output width (default: input width)\n  --outH <int>       output height (default: input height)\n  --background <mode> edge|black|transparent (default edge)\n");
    }

    // ------------------------------------------------------------
    // Input parsing helpers
    // ------------------------------------------------------------

    private static Map<String,String> parseArgs(String[] args) {
        Map<String,String> m = new LinkedHashMap<>();
        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if (a.startsWith("--")) {
                String v = (i+1 < args.length && !args[i+1].startsWith("--")) ? args[++i] : "true";
                m.put(a, v);
            }
        }
        return m;
    }

    static class Pair {
        final double sx, sy; // source normalized [0,1]
        final double dx, dy; // dest   normalized [0,1]
        Pair(double sx, double sy, double dx, double dy) {
            this.sx = sx; this.sy = sy; this.dx = dx; this.dy = dy;
        }
    }

    private static List<Pair> loadPairsFlexible(String path) throws IOException {
        String lower = path.toLowerCase(Locale.ROOT);
        if (lower.endsWith(".csv")) return loadPairsCSV(path);
        return loadPairsFromJsonLike(path);
    }

    private static List<Pair> loadPairsCSV(String path) throws IOException {
        List<Pair> list = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] t = line.split(",");
                if (t.length < 4) continue;
                double sx = Double.parseDouble(t[0].trim());
                double sy = Double.parseDouble(t[1].trim());
                double dx = Double.parseDouble(t[2].trim());
                double dy = Double.parseDouble(t[3].trim());
                list.add(new Pair(sx, sy, dx, dy));
            }
        }
        return list;
    }

    // Very lightweight extractor tailored to the provided JSON structure; avoids external libs.
    private static List<Pair> loadPairsFromJsonLike(String path) throws IOException {
        String s = new String(Files.readAllBytes(Paths.get(path)), StandardCharsets.UTF_8);
        // Greedy but robust enough for the supplied file: capture src.xn, src.yn, dst.xn, dst.yn in each pair object
      // DOESNT COMPILE, SO IGNORING Pattern pat = Pattern.compile("\\"src\\"\\s*:\\s*\\{[^}]*?\\"xn\\"\\s*:\\s*([\\-0-9.eE]+)\\s*,\\s*\\"yn\\"\\s*:\\s*([\\-0-9.eE]+)[^}]*?}[^}]*?\\"dst\\"\\s*:\\s*\\{[^}]*?\\"xn\\"\\s*:\\s*([\\-0-9.eE]+)\\s*,\\s*\\"yn\\"\\s*:\\s*([\\-0-9.eE]+)",Pattern.DOTALL);
        Pattern pat = null;
        Matcher m = pat.matcher(s);
        List<Pair> list = new ArrayList<>();
        while (m.find()) {
            double sx = Double.parseDouble(m.group(1));
            double sy = Double.parseDouble(m.group(2));
            double dx = Double.parseDouble(m.group(3));
            double dy = Double.parseDouble(m.group(4));
            list.add(new Pair(sx, sy, dx, dy));
        }
        return list;
    }

    // ------------------------------------------------------------
    // Thin Plate Spline (2D) with simple regularization
    // ------------------------------------------------------------

    static class ThinPlateSpline2D {
        private final double[][] pts; // N x 2 (control points in domain)
        private final double[] wX;    // weights for x mapping (N + 3)
        private final double[] wY;    // weights for y mapping (N + 3)

        private ThinPlateSpline2D(double[][] pts, double[] wX, double[] wY) {
            this.pts = pts; this.wX = wX; this.wY = wY;
        }

        public static ThinPlateSpline2D fit(double[][] srcPts, double[] dstX, double[] dstY, double lambda) {
            int N = srcPts.length;
            int M = N + 3; // affine terms
            double[][] A = new double[M][M];

            // Fill K block (N x N)
            for (int i = 0; i < N; i++) {
                A[i][i] = lambda; // regularization on diagonal
                for (int j = i + 1; j < N; j++) {
                    double r = dist(srcPts[i], srcPts[j]);
                    double u = U(r);
                    A[i][j] = u;
                    A[j][i] = u;
                }
            }

            // Fill P (N x 3) and P^T
            for (int i = 0; i < N; i++) {
                double x = srcPts[i][0];
                double y = srcPts[i][1];
                A[i][N + 0] = 1.0;
                A[i][N + 1] = x;
                A[i][N + 2] = y;

                A[N + 0][i] = 1.0;
                A[N + 1][i] = x;
                A[N + 2][i] = y;
            }
            // Lower right 3x3 is zeros

            // Build right hand sides
            double[] bx = new double[M];
            double[] by = new double[M];
            for (int i = 0; i < N; i++) {
                bx[i] = dstX[i];
                by[i] = dstY[i];
            }
            // affine RHS entries already 0

            double[] wx = solveSymmetric(A, bx);
            double[] wy = solveSymmetric(A, by);

            return new ThinPlateSpline2D(srcPts, wx, wy);
        }

        public double[] map(double x, double y) {
            int N = pts.length;
            double fx = wX[N + 0] + wX[N + 1] * x + wX[N + 2] * y;
            double fy = wY[N + 0] + wY[N + 1] * x + wY[N + 2] * y;
            for (int i = 0; i < N; i++) {
                double u = U(dist2(x, y, pts[i][0], pts[i][1]));
                fx += wX[i] * u;
                fy += wY[i] * u;
            }
            return new double[]{fx, fy};
        }

        private static double dist(double[] a, double[] b) {
            double dx = a[0] - b[0];
            double dy = a[1] - b[1];
            return Math.sqrt(dx * dx + dy * dy);
        }
        private static double dist2(double x1, double y1, double x2, double y2) {
            double dx = x1 - x2, dy = y1 - y2; return Math.sqrt(dx*dx + dy*dy);
        }
        private static double U(double r) {
            if (r <= 1e-12) return 0.0; // limit r^2 log r^2 -> 0
            double r2 = r * r;
            return r2 * Math.log(r2);
        }

        // Solve A x = b with partial pivoting Gaussian elimination
        private static double[] solveSymmetric(double[][] A, double[] b) {
            int n = b.length;
            double[][] M = new double[n][n+1];
            for (int i = 0; i < n; i++) {
                System.arraycopy(A[i], 0, M[i], 0, n);
                M[i][n] = b[i];
            }
            // Forward elim with pivoting
            for (int p = 0; p < n; p++) {
                int max = p;
                double vmax = Math.abs(M[p][p]);
                for (int r = p + 1; r < n; r++) {
                    double v = Math.abs(M[r][p]);
                    if (v > vmax) { vmax = v; max = r; }
                }
                if (Math.abs(M[max][p]) < 1e-18) throw new RuntimeException("Singular matrix at pivot " + p);
                if (max != p) {
                    double[] tmp = M[p]; M[p] = M[max]; M[max] = tmp;
                }
                double pivot = M[p][p];
                for (int c = p; c <= n; c++) M[p][c] /= pivot;
                for (int r = 0; r < n; r++) if (r != p) {
                    double f = M[r][p];
                    if (f != 0.0) {
                        for (int c = p; c <= n; c++) M[r][c] -= f * M[p][c];
                    }
                }
            }
            double[] x = new double[n];
            for (int i = 0; i < n; i++) x[i] = M[i][n];
            return x;
        }
    }

    // ------------------------------------------------------------
    // Inverse mapping grid (dst -> src) for speed
    // ------------------------------------------------------------

    static class InverseGrid {
        final int gx, gy;          // grid nodes in x and y
        final double[][] sx, sy;   // [gy][gx] source normalized coords at each dst grid node

        private InverseGrid(int gx, int gy) {
            this.gx = gx; this.gy = gy;
            this.sx = new double[gy][gx];
            this.sy = new double[gy][gx];
        }

        static InverseGrid build(ThinPlateSpline2D inv, int gx, int gy) {
            InverseGrid g = new InverseGrid(gx, gy);
            for (int j = 0; j < gy; j++) {
                double yn = (gy == 1) ? 0.0 : (j / (double)(gy - 1));
                for (int i = 0; i < gx; i++) {
                    double xn = (gx == 1) ? 0.0 : (i / (double)(gx - 1));
                    double[] s = inv.map(xn, yn);
                    g.sx[j][i] = s[0];
                    g.sy[j][i] = s[1];
                }
            }
            return g;
        }

        // Bilinear interpolate the grid
        double[] map(double xn, double yn) {
            double fx = (gx - 1) * xn;
            double fy = (gy - 1) * yn;
            int x0 = clamp((int)Math.floor(fx), 0, gx - 2);
            int y0 = clamp((int)Math.floor(fy), 0, gy - 2);
            double tx = fx - x0;
            double ty = fy - y0;

            double s00x = sx[y0][x0];
            double s10x = sx[y0][x0 + 1];
            double s01x = sx[y0 + 1][x0];
            double s11x = sx[y0 + 1][x0 + 1];

            double s00y = sy[y0][x0];
            double s10y = sy[y0][x0 + 1];
            double s01y = sy[y0 + 1][x0];
            double s11y = sy[y0 + 1][x0 + 1];

            double sxInterp = lerp(lerp(s00x, s10x, tx), lerp(s01x, s11x, tx), ty);
            double syInterp = lerp(lerp(s00y, s10y, tx), lerp(s01y, s11y, tx), ty);
            return new double[]{sxInterp, syInterp};
        }

        private static int clamp(int v, int lo, int hi) { return (v < lo) ? lo : (v > hi) ? hi : v; }
        private static double lerp(double a, double b, double t) { return a + (b - a) * t; }
    }

    // ------------------------------------------------------------
    // Image sampling
    // ------------------------------------------------------------

    static class Sampler {
        final BufferedImage img;
        final String mode; // edge | black | transparent
        final boolean hasAlpha;
        Sampler(BufferedImage img, String mode) {
            this.img = img; this.mode = mode; this.hasAlpha = mode.equals("transparent");
        }
        int sampleBilinearRGBA(double x, double y) {
            int w = img.getWidth(), h = img.getHeight();
            if (Double.isNaN(x) || Double.isNaN(y)) return fallback();
            if (x < 0 || y < 0 || x > w - 1 || y > h - 1) {
                switch (mode) {
                    case "black": return 0xFF000000;
                    case "transparent": return 0x00000000;
                    default: // edge clamp
                        x = clamp(x, 0, w - 1);
                        y = clamp(y, 0, h - 1);
                }
            }
            int x0 = (int)Math.floor(x), y0 = (int)Math.floor(y);
            int x1 = Math.min(x0 + 1, w - 1), y1 = Math.min(y0 + 1, h - 1);
            double tx = x - x0, ty = y - y0;

            int c00 = img.getRGB(x0, y0);
            int c10 = img.getRGB(x1, y0);
            int c01 = img.getRGB(x0, y1);
            int c11 = img.getRGB(x1, y1);

            double[] a = rgba(c00), b = rgba(c10), c = rgba(c01), d = rgba(c11);
            double[] top = mix(a, b, tx);
            double[] bot = mix(c, d, tx);
            double[] m = mix(top, bot, ty);
            int A = (int)clamp(Math.round(m[0]), 0, 255);
            int R = (int)clamp(Math.round(m[1]), 0, 255);
            int G = (int)clamp(Math.round(m[2]), 0, 255);
            int B = (int)clamp(Math.round(m[3]), 0, 255);
            return (A << 24) | (R << 16) | (G << 8) | B;
        }
        private static double[] rgba(int c) {
            return new double[]{(c >>> 24) & 0xFF, (c >>> 16) & 0xFF, (c >>> 8) & 0xFF, c & 0xFF};
        }
        private static double[] mix(double[] a, double[] b, double t) {
            return new double[]{
                    a[0] + (b[0] - a[0]) * t,
                    a[1] + (b[1] - a[1]) * t,
                    a[2] + (b[2] - a[2]) * t,
                    a[3] + (b[3] - a[3]) * t
            };
        }
        private static double clamp(double v, double lo, double hi) { return v < lo ? lo : Math.min(v, hi); }
        private int fallback() { return mode.equals("transparent") ? 0x00000000 : 0xFF000000; }
    }
}
