package defish;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.*;

/** Center-origin radial model:
 *  r_dst(r,phi) = f0(r) + c2(r)*|sin(phi)|^2 + c4(r)*|sin(phi)|^4
 *  All radii normalized by half-height (H/2). Angles in radians.
 */
public class HemiRadialModel {

    private final double[] rGrid, f0, c2, c4; // monotone vs r
    public final double originX = 0.5, originY = 0.5; // fraction of width/height

    public HemiRadialModel(double[] rGrid, double[] f0, double[] c2, double[] c4) {
        this.rGrid = rGrid; this.f0 = f0; this.c2 = c2; this.c4 = c4;
        if (rGrid.length != f0.length || f0.length != c2.length || c2.length != c4.length)
            throw new IllegalArgumentException("Arrays must have equal length");
    }

    /** Load from the provided JSON produced above (no external libs). */
    public static HemiRadialModel load(File jsonFile) throws IOException {
        String txt = new String(java.nio.file.Files.readAllBytes(jsonFile.toPath()), StandardCharsets.UTF_8);
        double[] r = readArray(txt, "\"r_grid\"");
        double[] f0 = readArray(txt, "\"f0\"");
        double[] c2 = readArray(txt, "\"c2\"");
        double[] c4 = readArray(txt, "\"c4\"");
        return new HemiRadialModel(r,f0,c2,c4);
    }

    private static double[] readArray(String txt, String key) {
        Matcher m = Pattern.compile(key + "\\s*:\\s*\\[(.*?)\\]", Pattern.DOTALL).matcher(txt);
        if (!m.find()) throw new IllegalArgumentException("Missing array: " + key);
        String body = m.group(1);
        String[] parts = body.split("[,\\s]+");
        ArrayList<Double> vals = new ArrayList<>();
        for (String p : parts) if (!p.isEmpty()) vals.add(Double.parseDouble(p));
        double[] out = new double[vals.size()];
        for (int i=0;i<out.length;i++) out[i]=vals.get(i);
        return out;
    }

    /** Piecewise-linear interpolation y(x) with clamp. */
    private static double interp(double x, double[] xs, double[] ys) {
        if (x <= xs[0]) return ys[0];
        int n = xs.length;
        if (x >= xs[n-1]) return ys[n-1];
        int lo = 0, hi = n-1;
        while (hi - lo > 1) {
            int mid = (lo+hi)>>>1;
            if (xs[mid] <= x) lo = mid; else hi = mid;
        }
        double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + t * (ys[hi] - ys[lo]);
    }

    /** Evaluate r_dst for given normalized radius r and angle phi (radians). */
    public double mapRadius(double r, double phi) {
        if (r < 0) r = 0;
        double s = Math.abs(Math.sin(phi));
        double base = interp(r, rGrid, f0);
        double a2 = interp(r, rGrid, c2);
        double a4 = interp(r, rGrid, c4);
        double rd = base + a2*(s*s) + a4*(s*s*s*s);
        // enforce monotone / non-negative
        return rd < 0 ? 0 : rd;
    }

    /** Generate one curve (r_src -> r_dst) with N radial samples up to rMax. */
    public double[][] generateCurve(int samples, double phi, double rMax) {
        double[][] out = new double[2][samples]; // [0]=r_src, [1]=r_dst
        double step = rMax / (samples - 1);
        double last = 0;
        for (int i=0;i<samples;i++) {
            double rs = i * step;
            double rd = mapRadius(rs, phi);
            // enforce monotone along this curve
            if (rd < last) rd = last;
            out[0][i] = rs; out[1][i] = rd; last = rd;
        }
        return out;
    }

    /** Example CLI: emit K curves (angles) with S samples to CSV. */
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.out.println("Usage: java HemiRadialModel <radial_model.json> [curves=4000] [samples=129] [rMax=1.7] [out=curves.csv]");
            return;
        }
        HemiRadialModel model = HemiRadialModel.load(new File(args[0]));
        int curves = args.length>1 ? Integer.parseInt(args[1]) : 4000;
        int samples= args.length>2 ? Integer.parseInt(args[2]) : 129;
        double rMax = args.length>3 ? Double.parseDouble(args[3]) : 1.7;
        String out = args.length>4 ? args[4] : "curves.csv";

        try (PrintWriter pw = new PrintWriter(out, StandardCharsets.UTF_8)) {
            pw.println("angle_deg,r_src,r_dst");
            for (int j=0;j<curves;j++) {
                double phi = -Math.PI + (2*Math.PI)*j/curves;
                double angDeg = Math.toDegrees(phi);
                double[][] curve = model.generateCurve(samples, phi, rMax);
                for (int i=0;i<samples;i++)
                    pw.printf(Locale.US,"%.6f,%.6f,%.6f%n", angDeg, curve[0][i], curve[1][i]);
            }
        }
        System.out.println("Wrote " + out);
    }
}
