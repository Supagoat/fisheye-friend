package defish;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import javax.imageio.ImageIO;

/**
 * HemiCurveLab — a lightweight Swing tool to help reverse engineer Fisheye Hemi mappings.
 *
 * Features
 *  - Load BEFORE and AFTER images and view them side-by-side (scrollable at 1:1 scale).
 *  - On the BEFORE image, add full-span vertical or horizontal lines with a click.
 *  - For every BEFORE line, an AFTER-side cubic Bézier curve is created (initially straight).
 *  - On the AFTER image, drag the Bézier control points to adjust the curve.
 *  - Save/Load the curve set to/from JSON (no external libs; minimal JSON writer/reader included).
 *  - Export a JSON file of BEFORE grid intersections mapped onto AFTER Bézier intersections.
 *
 * Notes
 *  - Curves for vertical lines are parameterized by t = y / afterHeight.
 *    Curves for horizontal lines are parameterized by s = x / afterWidth.
 *  - When exporting intersections, the AFTER point for (vertical i, horizontal j)
 *    is computed by evaluating the vertical curve i at the normalized y of the
 *    BEFORE horizontal j (mapped to AFTER using y_n * afterHeight). This keeps
 *    the Y in AFTER consistent with the BEFORE intersection while letting the
 *    curve supply the X (and vice versa logic can be added if desired).
 */
public class HemiCurveLab extends JFrame {
    // Images
    private BufferedImage beforeImg, afterImg;
    private File beforeFile, afterFile;

    // Data structures
    private final java.util.List<VLine> vlines = new ArrayList<>();
    private final java.util.List<HLine> hlines = new ArrayList<>();
    private int nextId = 1;

    // UI
    private final BeforePanel beforePanel = new BeforePanel();
    private final AfterPanel afterPanel = new AfterPanel();

    private Mode mode = Mode.ADD_VERTICAL; // default
    private JLabel status = new JLabel(" ");

    enum Mode { ADD_VERTICAL, ADD_HORIZONTAL, SELECT_AFTER }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            HemiCurveLab app = new HemiCurveLab();
            app.setVisible(true);
        });
    }

    public HemiCurveLab() {
        super("HemiCurveLab — BEFORE lines to AFTER Bézier curves");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Menus
        setJMenuBar(buildMenuBar());

        // Toolbar
        JToolBar tb = new JToolBar();
        JButton vb = new JButton("Add Vertical (V)");
        vb.addActionListener(e -> setMode(Mode.ADD_VERTICAL));
        JButton hb = new JButton("Add Horizontal (H)");
        hb.addActionListener(e -> setMode(Mode.ADD_HORIZONTAL));
        JButton sb = new JButton("Select/Drag After (S)");
        sb.addActionListener(e -> setMode(Mode.SELECT_AFTER));
        tb.add(vb); tb.add(hb); tb.add(sb);
        add(tb, BorderLayout.NORTH);

        // Split pane with scroll panes
        JScrollPane left = new JScrollPane(beforePanel);
        JScrollPane right = new JScrollPane(afterPanel);
        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, left, right);
        split.setResizeWeight(0.5);
        add(split, BorderLayout.CENTER);

        // Status bar
        status.setBorder(BorderFactory.createEmptyBorder(4,8,4,8));
        add(status, BorderLayout.SOUTH);

        // Shortcuts
        setupKeyBindings();

        setSize(1280, 800);
        setLocationRelativeTo(null);
        setMode(Mode.ADD_VERTICAL);
    }

    private JMenuBar buildMenuBar() {
        JMenuBar mb = new JMenuBar();

        JMenu file = new JMenu("File");
        JMenuItem loadBefore = new JMenuItem("Load BEFORE image...");
        loadBefore.addActionListener(e -> doLoadBefore());
        JMenuItem loadAfter = new JMenuItem("Load AFTER image...");
        loadAfter.addActionListener(e -> doLoadAfter());
        JMenuItem saveCurves = new JMenuItem("Save curves to JSON...");
        saveCurves.addActionListener(e -> doSaveCurves());
        JMenuItem loadCurves = new JMenuItem("Load curves from JSON...");
        loadCurves.addActionListener(e -> doLoadCurves());
        JMenuItem export = new JMenuItem("Export intersections JSON...");
        export.addActionListener(e -> doExportIntersections());
        JMenuItem quit = new JMenuItem("Quit");
        quit.addActionListener(e -> System.exit(0));
        file.add(loadBefore);
        file.add(loadAfter);
        file.addSeparator();
        file.add(saveCurves);
        file.add(loadCurves);
        file.addSeparator();
        file.add(export);
        file.addSeparator();
        file.add(quit);

        JMenu help = new JMenu("Help");
        JMenuItem about = new JMenuItem("Usage Tips");
        about.addActionListener(e -> showTips());
        help.add(about);

        mb.add(file);
        mb.add(help);
        return mb;
    }

    private void showTips() {
        String msg = String.join("\n",
            "Workflow:",
            "  1) Load BEFORE and AFTER images.",
            "  2) Press V and click on BEFORE to add vertical lines (full height).",
            "  3) Press H and click on BEFORE to add horizontal lines (full width).",
            "  4) Press S, then drag Bézier control points on the AFTER view.",
            "  5) File to Save curves as JSON / Load curves from JSON.",
            "  6) File to Export intersections JSON to generate BEFORE between AFTER mappings.",
            "Notes:",
            "  • Vertical curve param is t = y_after / afterHeight.",
            "  • Horizontal curve param is s = x_after / afterWidth.",
            "  • Intersections use BEFORE normalized positions mapped onto AFTER.");
        JOptionPane.showMessageDialog(this, msg, "HemiCurveLab", JOptionPane.INFORMATION_MESSAGE);
    }

    private void setupKeyBindings() {
        JComponent root = (JComponent)getContentPane();
        root.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke('V'), "mode_v");
        root.getActionMap().put("mode_v", new AbstractAction(){public void actionPerformed(ActionEvent e){setMode(Mode.ADD_VERTICAL);}});
        root.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke('H'), "mode_h");
        root.getActionMap().put("mode_h", new AbstractAction(){public void actionPerformed(ActionEvent e){setMode(Mode.ADD_HORIZONTAL);}});
        root.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke('S'), "mode_s");
        root.getActionMap().put("mode_s", new AbstractAction(){public void actionPerformed(ActionEvent e){setMode(Mode.SELECT_AFTER);}});
    }

    private void setMode(Mode m) {
        mode = m;
        switch (m) {
            case ADD_VERTICAL: status.setText("Mode: Add Vertical — click on BEFORE to add a vertical line"); break;
            case ADD_HORIZONTAL: status.setText("Mode: Add Horizontal — click on BEFORE to add a horizontal line"); break;
            case SELECT_AFTER: status.setText("Mode: Select/Drag on AFTER — drag curve control points"); break;
        }
        beforePanel.repaint();
        afterPanel.repaint();
    }

    // UI actions
    private void doLoadBefore() {
        File f = pickImageFile("Select BEFORE image");
        if (f == null) return;
        try {
            beforeImg = ImageIO.read(f);
            beforeFile = f;
            beforePanel.setPreferredSize(new Dimension(beforeImg.getWidth(), beforeImg.getHeight()));
            beforePanel.revalidate();
            beforePanel.repaint();
        } catch (IOException ex) { showErr(ex); }
    }

    private void doLoadAfter() {
        File f = pickImageFile("Select AFTER image");
        if (f == null) return;
        try {
            afterImg = ImageIO.read(f);
            afterFile = f;
            afterPanel.setPreferredSize(new Dimension(afterImg.getWidth(), afterImg.getHeight()));
            afterPanel.revalidate();
            afterPanel.repaint();
        } catch (IOException ex) { showErr(ex); }
    }

    private File pickImageFile(String title) {
        JFileChooser fc = new JFileChooser();
        fc.setDialogTitle(title);
        fc.setFileFilter(new FileNameExtensionFilter("Images", "png", "jpg", "jpeg", "bmp"));
        int r = fc.showOpenDialog(this);
        return (r == JFileChooser.APPROVE_OPTION) ? fc.getSelectedFile() : null;
    }

    private void doSaveCurves() {
        if (!imagesLoaded()) return;
        JFileChooser fc = new JFileChooser();
        fc.setDialogTitle("Save curves JSON");
        int r = fc.showSaveDialog(this);
        if (r != JFileChooser.APPROVE_OPTION) return;
        File f = fc.getSelectedFile();
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(f), StandardCharsets.UTF_8))) {
            pw.print(serializeCurves());
        } catch (IOException ex) { showErr(ex); }
    }

    private void doLoadCurves() {
        JFileChooser fc = new JFileChooser();
        fc.setDialogTitle("Load curves JSON");
        int r = fc.showOpenDialog(this);
        if (r != JFileChooser.APPROVE_OPTION) return;
        File f = fc.getSelectedFile();
        try {
            String json = new String(java.nio.file.Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
            deserializeCurves(json);
            beforePanel.repaint();
            afterPanel.repaint();
        } catch (IOException ex) { showErr(ex); }
    }

    private void doExportIntersections() {
        if (!imagesLoaded()) return;
        JFileChooser fc = new JFileChooser();
        fc.setDialogTitle("Export intersections JSON");
        int r = fc.showSaveDialog(this);
        if (r != JFileChooser.APPROVE_OPTION) return;
        File f = fc.getSelectedFile();
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(f), StandardCharsets.UTF_8))) {
            pw.print(serializeIntersections());
        } catch (IOException ex) { showErr(ex); }
    }

    private boolean imagesLoaded() {
        if (beforeImg == null || afterImg == null) {
            JOptionPane.showMessageDialog(this, "Please load BOTH BEFORE and AFTER images first.", "Missing image", JOptionPane.WARNING_MESSAGE);
            return false;
        }
        return true;
    }

    private void showErr(Exception ex) {
        ex.printStackTrace();
        JOptionPane.showMessageDialog(this, ex.toString(), "Error", JOptionPane.ERROR_MESSAGE);
    }

    // ---------------- Panels -----------------
    private class BeforePanel extends JPanel {
        BeforePanel() {
            setBackground(Color.DARK_GRAY);
            addMouseListener(new MouseAdapter() {
                @Override public void mousePressed(MouseEvent e) {
                    if (beforeImg == null) return;
                    if (mode == Mode.ADD_VERTICAL) {
                        int x = clamp(e.getX(), 0, beforeImg.getWidth()-1);
                        addVerticalLine(x);
                    } else if (mode == Mode.ADD_HORIZONTAL) {
                        int y = clamp(e.getY(), 0, beforeImg.getHeight()-1);
                        addHorizontalLine(y);
                    }
                }
            });
        }
        @Override public Dimension getPreferredSize() {
            return (beforeImg != null) ? new Dimension(beforeImg.getWidth(), beforeImg.getHeight()) : new Dimension(640,480);
        }
        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            if (beforeImg != null) g2.drawImage(beforeImg, 0,0,null);
            // draw lines
            g2.setStroke(new BasicStroke(1f));
            g2.setColor(new Color(255, 100, 100));
            for (VLine vl : vlines) {
                g2.drawLine(vl.x, 0, vl.x, beforeImg.getHeight());
            }
            g2.setColor(new Color(100, 180, 255));
            for (HLine hl : hlines) {
                g2.drawLine(0, hl.y, beforeImg.getWidth(), hl.y);
            }
        }
    }

    private class AfterPanel extends JPanel {
        CurveHandle active; // current drag target
        AfterPanel() {
            setBackground(Color.DARK_GRAY);
            MouseAdapter ma = new MouseAdapter() {
                @Override public void mousePressed(MouseEvent e) {
                    if (afterImg == null) return;
                    if (mode != Mode.SELECT_AFTER) return;
                    active = findHandleNear(e.getPoint());
                }
                @Override public void mouseDragged(MouseEvent e) {
                    if (afterImg == null) return;
                    if (active != null) {
                        active.dragTo(e.getPoint());
                        repaint();
                    }
                }
                @Override public void mouseReleased(MouseEvent e) { active = null; }
            };
            addMouseListener(ma);
            addMouseMotionListener(ma);
        }
        @Override public Dimension getPreferredSize() {
            return (afterImg != null) ? new Dimension(afterImg.getWidth(), afterImg.getHeight()) : new Dimension(640,480);
        }
        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            if (afterImg != null) g2.drawImage(afterImg, 0,0,null);

            // Draw curves
            if (afterImg == null) return;
            int W = afterImg.getWidth();
            int H = afterImg.getHeight();

            // Vertical curves in red
            g2.setStroke(new BasicStroke(2f));
            g2.setColor(new Color(220,60,60));
            for (VLine vl : vlines) {
                Cubic c = vl.curve;
                if (c == null) continue;
                drawCubicV(g2, c, H);
                drawHandles(g2, c, true);
            }

            // Horizontal curves in blue
            g2.setColor(new Color(80,140,240));
            for (HLine hl : hlines) {
                Cubic c = hl.curve;
                if (c == null) continue;
                drawCubicH(g2, c, W);
                drawHandles(g2, c, false);
            }
        }

        private void drawCubicV(Graphics2D g2, Cubic c, int H) {
            Path2D path = new Path2D.Double();
            path.moveTo(c.p0.x, c.p0.y);
            path.curveTo(c.p1.x, c.p1.y, c.p2.x, c.p2.y, c.p3.x, c.p3.y);
            g2.draw(path);
        }
        private void drawCubicH(Graphics2D g2, Cubic c, int W) {
            Path2D path = new Path2D.Double();
            path.moveTo(c.p0.x, c.p0.y);
            path.curveTo(c.p1.x, c.p1.y, c.p2.x, c.p2.y, c.p3.x, c.p3.y);
            g2.draw(path);
        }
        private void drawHandles(Graphics2D g2, Cubic c, boolean vertical) {
            final int r = 4;
            // guides
            g2.setStroke(new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 1f, new float[]{5f,5f}, 0f));
            g2.setColor(new Color(0,0,0,80));
            g2.draw(new Line2D.Double(c.p0, c.p1));
            g2.draw(new Line2D.Double(c.p2, c.p3));
            // handles
            g2.setStroke(new BasicStroke(1.5f));
            g2.setColor(new Color(30,30,30));
            for (int i=0;i<4;i++) {
                Point2D p = c.get(i);
                g2.fill(new Ellipse2D.Double(p.getX()-r, p.getY()-r, r*2, r*2));
            }
        }

        private CurveHandle findHandleNear(Point p) {
            final double tol = 8.0;
            if (afterImg == null) return null;
            for (VLine vl : vlines) {
                if (vl.curve == null) continue;
                for (int i=0;i<4;i++) {
                    Point2D pt = vl.curve.get(i);
                    if (pt.distance(p) <= tol) return new CurveHandle(vl.curve, i, true, afterImg.getWidth(), afterImg.getHeight());
                }
            }
            for (HLine hl : hlines) {
                if (hl.curve == null) continue;
                for (int i=0;i<4;i++) {
                    Point2D pt = hl.curve.get(i);
                    if (pt.distance(p) <= tol) return new CurveHandle(hl.curve, i, false, afterImg.getWidth(), afterImg.getHeight());
                }
            }
            return null;
        }
    }

    // ----------- Model helpers -----------
    private static class Cubic {
        Point2D.Double p0, p1, p2, p3; // order matters
        Cubic(Point2D.Double p0, Point2D.Double p1, Point2D.Double p2, Point2D.Double p3){
            this.p0=p0; this.p1=p1; this.p2=p2; this.p3=p3;
        }
        Point2D get(int idx){
            switch(idx){
                case 0: return p0; case 1: return p1; case 2: return p2; case 3: return p3; default: return p0;
            }
        }
        // Evaluate cubic at parameter t in [0,1]
        Point2D.Double eval(double t){
            double u = 1-t;
            double b0 = u*u*u;
            double b1 = 3*u*u*t;
            double b2 = 3*u*t*t;
            double b3 = t*t*t;
            double x = b0*p0.x + b1*p1.x + b2*p2.x + b3*p3.x;
            double y = b0*p0.y + b1*p1.y + b2*p2.y + b3*p3.y;
            return new Point2D.Double(x,y);
        }
    }

    private static class CurveHandle {
        final Cubic cubic; final int idx; final boolean vertical; final int W,H;
        CurveHandle(Cubic c, int idx, boolean vertical, int W, int H){ this.cubic=c; this.idx=idx; this.vertical=vertical; this.W=W; this.H=H; }
        void dragTo(Point p){
            // Constrain endpoints to edges so curves always span full dimension
            if (vertical) {
                if (idx==0) { cubic.p0.x = clamp(p.x, 0, W); cubic.p0.y = 0; }
                else if (idx==3) { cubic.p3.x = clamp(p.x, 0, W); cubic.p3.y = H; }
                else if (idx==1) { cubic.p1.x = clamp(p.x, 0, W); cubic.p1.y = clamp(p.y, 0, H); }
                else if (idx==2) { cubic.p2.x = clamp(p.x, 0, W); cubic.p2.y = clamp(p.y, 0, H); }
            } else {
                if (idx==0) { cubic.p0.x = 0; cubic.p0.y = clamp(p.y, 0, H); }
                else if (idx==3) { cubic.p3.x = W; cubic.p3.y = clamp(p.y, 0, H); }
                else if (idx==1) { cubic.p1.x = clamp(p.x, 0, W); cubic.p1.y = clamp(p.y, 0, H); }
                else if (idx==2) { cubic.p2.x = clamp(p.x, 0, W); cubic.p2.y = clamp(p.y, 0, H); }
            }
        }
    }

    private class VLine {
        int id; int x; double xn; // BEFORE
        Cubic curve;              // AFTER
    }
    private class HLine {
        int id; int y; double yn; // BEFORE
        Cubic curve;              // AFTER
    }

    private void addVerticalLine(int xBefore) {
        if (!imagesLoaded()) return;
        VLine vl = new VLine();
        vl.id = nextId++;
        vl.x = xBefore;
        vl.xn = (double)xBefore / (double)beforeImg.getWidth();
        // Create straight vertical cubic on AFTER at same normalized x
        int W = afterImg.getWidth();
        int H = afterImg.getHeight();
        double xA = vl.xn * W;
        vl.curve = new Cubic(
            new Point2D.Double(xA, 0),
            new Point2D.Double(xA, H/3.0),
            new Point2D.Double(xA, 2*H/3.0),
            new Point2D.Double(xA, H)
        );
        vlines.add(vl);
        beforePanel.repaint();
        afterPanel.repaint();
    }

    private void addHorizontalLine(int yBefore) {
        if (!imagesLoaded()) return;
        HLine hl = new HLine();
        hl.id = nextId++;
        hl.y = yBefore;
        hl.yn = (double)yBefore / (double)beforeImg.getHeight();
        int W = afterImg.getWidth();
        int H = afterImg.getHeight();
        double yA = hl.yn * H;
        hl.curve = new Cubic(
            new Point2D.Double(0, yA),
            new Point2D.Double(W/3.0, yA),
            new Point2D.Double(2*W/3.0, yA),
            new Point2D.Double(W, yA)
        );
        hlines.add(hl);
        beforePanel.repaint();
        afterPanel.repaint();
    }

    // ---------- Serialization (no external libs) ----------
    private String serializeCurves() {
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"source\": {\"path\": \"").append(escape(beforeFile!=null?beforeFile.getAbsolutePath():""))
          .append("\", \"width\": ").append(beforeImg!=null?beforeImg.getWidth():0)
          .append(", \"height\": ").append(beforeImg!=null?beforeImg.getHeight():0).append("},\n");
        sb.append("  \"target\": {\"path\": \"").append(escape(afterFile!=null?afterFile.getAbsolutePath():""))
          .append("\", \"width\": ").append(afterImg!=null?afterImg.getWidth():0)
          .append(", \"height\": ").append(afterImg!=null?afterImg.getHeight():0).append("},\n");
        sb.append("  \"verticals\": [\n");
        for (int i=0;i<vlines.size();i++) {
            VLine vl = vlines.get(i);
            sb.append("    {");
            sb.append("\"id\": ").append(vl.id).append(", ");
            sb.append("\"before\": {\"x\": ").append(vl.x).append(", \"xn\": ").append(fmt(vl.xn)).append("}, ");
            appendCurve(sb, vl.curve, afterImg.getWidth(), afterImg.getHeight());
            sb.append("}");
            if (i < vlines.size()-1) sb.append(",");
            sb.append("\n");
        }
        sb.append("  ],\n");
        sb.append("  \"horizontals\": [\n");
        for (int i=0;i<hlines.size();i++) {
            HLine hl = hlines.get(i);
            sb.append("    {");
            sb.append("\"id\": ").append(hl.id).append(", ");
            sb.append("\"before\": {\"y\": ").append(hl.y).append(", \"yn\": ").append(fmt(hl.yn)).append("}, ");
            appendCurve(sb, hl.curve, afterImg.getWidth(), afterImg.getHeight());
            sb.append("}");
            if (i < hlines.size()-1) sb.append(",");
            sb.append("\n");
        }
        sb.append("  ]\n");
        sb.append("}\n");
        return sb.toString();
    }

    private void appendCurve(StringBuilder sb, Cubic c, int W, int H) {
        sb.append("\"curve\": {");
        sb.append("\"p0\": ").append(ptJson(c.p0, W, H)).append(", ");
        sb.append("\"p1\": ").append(ptJson(c.p1, W, H)).append(", ");
        sb.append("\"p2\": ").append(ptJson(c.p2, W, H)).append(", ");
        sb.append("\"p3\": ").append(ptJson(c.p3, W, H));
        sb.append("}");
    }

    private String ptJson(Point2D.Double p, int W, int H) {
        return "{" +
            "\"x\": " + fmt(p.x) + ", \"y\": " + fmt(p.y) + ", " +
            "\"xn\": " + fmt(W>0 ? p.x/W : 0) + ", \"yn\": " + fmt(H>0 ? p.y/H : 0) +
            "}";
    }

    private void deserializeCurves(String json) {
        // Very light-weight loader for files produced by serializeCurves().
        // Clears current lines and rebuilds from JSON content.
        vlines.clear(); hlines.clear(); nextId = 1;
        int Wt = afterImg!=null?afterImg.getWidth():0;
        int Ht = afterImg!=null?afterImg.getHeight():0;

        // Parse verticals
        String verts = extractArray(json, "\"verticals\"");
        if (verts != null) {
            java.util.List<String> objs = splitObjects(verts);
            for (String o : objs) {
                int id = (int)extractNumber(o, "\"id\"");
                double x = extractNumber(o, "\"before\"\\s*:\\s*\\{[^}] *\\\"x\\\"\s*:\s*", true);
                double xn = extractNumber(o, "\"before\"\\s*:\\s*\\{[^}] *\\\"xn\\\"\s*:\s*", true);
                VLine vl = new VLine();
                vl.id = id; nextId = Math.max(nextId, id+1);
                vl.x = (int)Math.round(x);
                vl.xn = xn;
                vl.curve = parseCurve(o, Wt, Ht);
                vlines.add(vl);
            }
        }
        // Parse horizontals
        String hors = extractArray(json, "\"horizontals\"");
        if (hors != null) {
            java.util.List<String> objs = splitObjects(hors);
            for (String o : objs) {
                int id = (int)extractNumber(o, "\"id\"");
                double y = extractNumber(o, "\"before\"\\s*:\\s*\\{[^}] *\\\"y\\\"\s*:\s*", true);
                double yn = extractNumber(o, "\"before\"\\s*:\\s*\\{[^}] *\\\"yn\\\"\s*:\s*", true);
                HLine hl = new HLine();
                hl.id = id; nextId = Math.max(nextId, id+1);
                hl.y = (int)Math.round(y);
                hl.yn = yn;
                hl.curve = parseCurve(o, Wt, Ht);
                hlines.add(hl);
            }
        }
    }

    private Cubic parseCurve(String obj, int W, int H) {
        Point2D.Double p0 = parsePoint(obj, "p0", W, H);
        Point2D.Double p1 = parsePoint(obj, "p1", W, H);
        Point2D.Double p2 = parsePoint(obj, "p2", W, H);
        Point2D.Double p3 = parsePoint(obj, "p3", W, H);
        return new Cubic(p0,p1,p2,p3);
    }

    private Point2D.Double parsePoint(String obj, String key, int W, int H) {
        String blk = extractObject(obj, "\""+key+"\"\s*:\s*\\{");
        double x = extractNumber(blk, "\"x\"");
        double y = extractNumber(blk, "\"y\"");
        return new Point2D.Double(x,y);
    }

    // ---------- Intersections export ----------
    private String serializeIntersections() {
        int Ws = beforeImg.getWidth();
        int Hs = beforeImg.getHeight();
        int Wt = afterImg.getWidth();
        int Ht = afterImg.getHeight();

        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"created_at\": \"").append(new java.util.Date()).append("\",\n");
        sb.append("  \"source\": {\"path\": \"").append(escape(beforeFile!=null?beforeFile.getAbsolutePath():""))
          .append("\", \"width\": ").append(Ws).append(", \"height\": ").append(Hs).append("},\n");
        sb.append("  \"target\": {\"path\": \"").append(escape(afterFile!=null?afterFile.getAbsolutePath():""))
          .append("\", \"width\": ").append(Wt).append(", \"height\": ").append(Ht).append("},\n");
        sb.append("  \"pairs\": [\n");
        int id = 1; // fresh ids for pairs
        for (VLine vl : vlines) {
            double xn = (double)vl.x / Ws;
            for (HLine hl : hlines) {
                double yn = (double)hl.y / Hs;
                // BEFORE intersection
                double sx = vl.x;
                double sy = hl.y;
                // AFTER intersection from vertical curve at same normalized y
                double yAfter = yn * Ht;
                double t = yAfter / Ht; // equals yn
                Point2D.Double p = vl.curve.eval(t);
                double dx = p.x;
                double dy = yAfter; // keep exact y from normalized value

                sb.append("    {\"id\": ").append(id++).append(", ")
                  .append("\"src\": {\"x\": ").append((int)Math.round(sx)).append(", \"y\": ").append((int)Math.round(sy))
                  .append(", \"xn\": ").append(fmt(xn)).append(", \"yn\": ").append(fmt(yn)).append("}, ")
                  .append("\"dst\": {\"x\": ").append((int)Math.round(dx)).append(", \"y\": ").append((int)Math.round(dy))
                  .append(", \"xn\": ").append(fmt(dx/Wt)).append(", \"yn\": ").append(fmt(dy/Ht)).append("}}");
                if (!(vl == vlines.get(vlines.size()-1) && hl == hlines.get(hlines.size()-1))) sb.append(",");
                sb.append("\n");
            }
        }
        sb.append("  ]\n}");
        return sb.toString();
    }

    // ----------- Tiny JSON helpers (very limited) -----------
    private static String escape(String s) { return s.replace("\\", "\\\\").replace("\"", "\\\""); }
    private static String fmt(double v) { return String.format(java.util.Locale.US, "%.6f", v); }
    private static int clamp(int v, int lo, int hi) { return Math.max(lo, Math.min(hi, v)); }

    private static String extractArray(String json, String key) {
        int k = json.indexOf(key);
        if (k < 0) return null;
        int lb = json.indexOf('[', k);
        if (lb < 0) return null;
        int depth = 0;
        for (int i = lb; i < json.length(); i++) {
            char c = json.charAt(i);
            if (c=='[') depth++;
            else if (c==']') { depth--; if (depth==0) return json.substring(lb+1, i); }
        }
        return null;
    }
    private static String extractObject(String json, String startKeyRegex) {
        // find start using regex, then return {...}
        java.util.regex.Pattern p = java.util.regex.Pattern.compile(startKeyRegex);
        java.util.regex.Matcher m = p.matcher(json);
        if (!m.find()) return "{}";
        int lb = json.indexOf('{', m.end()-1);
        if (lb < 0) return "{}";
        int depth = 0;
        for (int i = lb; i < json.length(); i++) {
            char c = json.charAt(i);
            if (c=='{') depth++;
            else if (c=='}') { depth--; if (depth==0) return json.substring(lb, i+1); }
        }
        return "{}";
    }
    private static java.util.List<String> splitObjects(String inner) {
        java.util.List<String> out = new ArrayList<>();
        int i = 0; int n = inner.length();
        while (i < n) {
            int lb = inner.indexOf('{', i);
            if (lb < 0) break;
            int depth = 0; int j = lb;
            for (; j < n; j++) {
                char c = inner.charAt(j);
                if (c=='{') depth++;
                else if (c=='}') { depth--; if (depth==0) { j++; break; } }
            }
            if (j > lb) out.add(inner.substring(lb, j));
            i = j + 1;
        }
        return out;
    }
    private static double extractNumber(String json, String key) { return extractNumber(json, key, false); }
    private static double extractNumber(String json, String key, boolean isRegexPrefix) {
        String regex = isRegexPrefix ? key + "([+-]?(?:\\d+\\.\\d*|\\d*\\.\\d+|\\d+))" : "\\\""+key+"\\\"\\s*:\\s*([+-]?(?:\\d+\\.\\d*|\\d*\\.\\d+|\\d+))";
        java.util.regex.Pattern p = java.util.regex.Pattern.compile(regex);
        java.util.regex.Matcher m = p.matcher(json);
        if (m.find()) {
            try { return Double.parseDouble(m.group(1)); } catch (Exception ignored) {}
        }
        return 0.0;
    }
}
