
package defish;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.function.IntConsumer;

/**
 * HemiPointMapper
 * ----------------
 * A two-pane point-pair annotation tool for reverse-engineering Fisheye-Hemi.
 * - Open Source (fisheye) and Target (Hemi output) images
 * - Place paired nodes: click Source, then Target -> creates pair #n
 * - Edit mode: drag points to adjust
 * - Hold LEFT SHIFT for 4x magnifier under cursor (per-pane)
 * - Right-click near a point -> delete that pair
 * - Save mapping to JSON with absolute and normalized coords (easy to parse later)
 */
public class HemiPointMapper extends JFrame {

    // --- App state ---
    private final PairManager pairs = new PairManager();
    private final ImagePane srcPane = new ImagePane("Source", Side.SOURCE, pairs);
    private final ImagePane dstPane = new ImagePane("Target", Side.TARGET, pairs);

    private final JLabel status = new JLabel("Open both images to begin.");
    private final JToggleButton placeModeBtn = new JToggleButton("Placement", true);
    private final JToggleButton editModeBtn  = new JToggleButton("Edit");

    private File srcPath = null, dstPath = null;

    public static void main(String[] args) {
        SwingUtilities.invokeLater(HemiPointMapper::new);
    }

    public HemiPointMapper() {
        super("Hemi Point Mapper — pair control points (Shift = 4× loupe)");
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setLayout(new BorderLayout(6,6));

        // Menu
        JMenuBar mb = new JMenuBar();
        JMenu file = new JMenu("File");
        JMenuItem openSrc = new JMenuItem("Open Source Image…");
        JMenuItem openDst = new JMenuItem("Open Target Image…");
        JMenuItem saveMap = new JMenuItem("Save Mapping…");
        JMenuItem exit    = new JMenuItem("Exit");
        openSrc.addActionListener(e -> openImage(Side.SOURCE));
        openDst.addActionListener(e -> openImage(Side.TARGET));
        saveMap.addActionListener(e -> saveMapping());
        exit.addActionListener(e -> dispose());
        file.add(openSrc); file.add(openDst); file.addSeparator(); file.add(saveMap); file.addSeparator(); file.add(exit);
        mb.add(file);

        JMenu view = new JMenu("View");
        JCheckBoxMenuItem showGrid = new JCheckBoxMenuItem("Show crosshair", true);
        showGrid.addActionListener(e -> { srcPane.setShowCrosshair(showGrid.isSelected()); dstPane.setShowCrosshair(showGrid.isSelected()); });
        view.add(showGrid);
        mb.add(view);

        setJMenuBar(mb);

        // Top controls
        JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 6));
        ButtonGroup modes = new ButtonGroup();
        modes.add(placeModeBtn); modes.add(editModeBtn);
        placeModeBtn.addActionListener(e -> setMode(Mode.PLACE));
        editModeBtn.addActionListener(e -> setMode(Mode.EDIT));
        top.add(placeModeBtn); top.add(editModeBtn);

        JButton clearBtn = new JButton("Clear All");
        clearBtn.addActionListener(e -> {
            if (JOptionPane.showConfirmDialog(this, "Clear all pairs?", "Confirm", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                pairs.clear();
                srcPane.repaint(); dstPane.repaint(); updateStatus();
            }
        });
        top.add(clearBtn);

        add(top, BorderLayout.NORTH);

        // Center panes side-by-side with scroll
        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                new JScrollPane(srcPane), new JScrollPane(dstPane));
        split.setResizeWeight(0.5);
        add(split, BorderLayout.CENTER);

        // Status
        status.setBorder(new EmptyBorder(4,8,4,8));
        add(status, BorderLayout.SOUTH);

        // Global key handling for SHIFT-loupe
        KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new KeyEventDispatcher() {
            @Override public boolean dispatchKeyEvent(KeyEvent e) {
                boolean down = (e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0;
                srcPane.setShiftDown(down);
                dstPane.setShiftDown(down);
                return false;
            }
        });
        JMenuItem loadItem = new JMenuItem("Load Map JSON…");
        loadItem.addActionListener(ev -> {
            javax.swing.JFileChooser fc = new javax.swing.JFileChooser();
            fc.setFileFilter(new javax.swing.filechooser.FileNameExtensionFilter("JSON", "json"));
            if (fc.showOpenDialog(this) == javax.swing.JFileChooser.APPROVE_OPTION) {
                try { loadMap(fc.getSelectedFile()); }
                catch (Exception ex) { ex.printStackTrace(); javax.swing.JOptionPane.showMessageDialog(this, ex.toString(), "Load failed", javax.swing.JOptionPane.ERROR_MESSAGE); }
            }
        });
        file.add(loadItem);     // add near your existing “Save Map…”

        
        setSize(1280, 800);
        setLocationRelativeTo(null);
        setVisible(true);
        setMode(Mode.PLACE);
    }

 // === JSON load helpers ===
    private static String slurp(File f) throws IOException {
        try (java.io.BufferedReader br = new java.io.BufferedReader(
                new java.io.InputStreamReader(new java.io.FileInputStream(f), java.nio.charset.StandardCharsets.UTF_8))) {
            StringBuilder sb = new StringBuilder(1<<20);
            String line; while ((line = br.readLine()) != null) sb.append(line).append('\n');
            return sb.toString();
        }
    }
    private static Double parseNullable(String s) { return (s == null ? null : Double.valueOf(s)); }
    private static int extractInt(String text, String regex) throws IOException {
        java.util.regex.Matcher m = java.util.regex.Pattern.compile(regex).matcher(text);
        if (!m.find()) throw new IOException("Missing field for regex: " + regex);
        return Integer.parseInt(m.group(1));
    }

 // Call AFTER both images are opened (so we know current sizes).
    private void loadMap(File jsonFile) throws Exception {

        String txt = slurp(jsonFile);

        // (we read sizes for logging only)
        int sw0 = extractInt(txt, "\"source\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
        int sh0 = extractInt(txt, "\"source\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");
        int dw0 = extractInt(txt, "\"target\"\\s*:\\s*\\{[^}]*?\"width\"\\s*:\\s*(\\d+)");
        int dh0 = extractInt(txt, "\"target\"\\s*:\\s*\\{[^}]*?\"height\"\\s*:\\s*(\\d+)");

        // Each pair; supports optional xn/yn fields.
        java.util.regex.Pattern P = java.util.regex.Pattern.compile(
            "\\{\\s*\"id\"\\s*:\\s*(\\d+)\\s*,\\s*\"src\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)\\s*,\\s*\"y\"\\s*:\\s*([\\-0-9.]+)"
          + "(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)\\s*,\\s*\"yn\"\\s*:\\s*([\\-0-9.]+))?[^}]*?\\}\\s*,\\s*\"dst\"\\s*:\\s*\\{[^}]*?\"x\"\\s*:\\s*([\\-0-9.]+)\\s*,\\s*\"y\"\\s*:\\s*([\\-0-9.]+)"
          + "(?:[^}]*?\"xn\"\\s*:\\s*([\\-0-9.]+)\\s*,\\s*\"yn\"\\s*:\\s*([\\-0-9.]+))?",
            java.util.regex.Pattern.DOTALL);

        java.util.regex.Matcher m = P.matcher(txt);

        java.util.List<NodePair> loaded = new java.util.ArrayList<>();
        //int cw = srcImg.getWidth(), ch = srcImg.getHeight();
      //  int dw = dstImg.getWidth(),  dh = dstImg.getHeight();

        while (m.find()) {
            int id = Integer.parseInt(m.group(1));

            double sxAbs = Double.parseDouble(m.group(2));
            double syAbs = Double.parseDouble(m.group(3));
            //double sxn = Double.parseDouble(m.group(4));
            //double syn = Double.parseDouble(m.group(5));

            double dxAbs = Double.parseDouble(m.group(6));
            double dyAbs = Double.parseDouble(m.group(7));
           // Double dxn = parseNullable(m.group(8));
            //Double dyn = parseNullable(m.group(9));

            // Prefer normalized if present; adapts to current image sizes.
           // double sxAbs = (sxn != null ? sxn * cw : sx);
            //double syAbs = (syn != null ? syn * ch : sy);
           // double dxAbs = (dxn != null ? dxn * dw : dx);
           // double dyAbs = (dyn != null ? dyn * dh : dy);

            loaded.add(new NodePair(id, new Point2D.Double(sxAbs, syAbs), new Point2D.Double(dxAbs, dyAbs)));
        }

        // Replace current pairs, re-id sequentially
        pairs.clear();
        loaded.sort(java.util.Comparator.comparingInt(p -> p.id));
        for (int i = 0; i < loaded.size(); i++) {
            NodePair p = loaded.get(i);
           //p.id = i;
            pairs.add(p);
        }

        status.setText("Loaded " + pairs.size() + " pairs (map src " + sw0 + "x" + sh0 + ", dst " + dw0 + "x" + dh0 + ")");
        repaint();
    }

    private void setMode(Mode m) {
        srcPane.setMode(m); dstPane.setMode(m);
        placeModeBtn.setSelected(m == Mode.PLACE);
        editModeBtn.setSelected(m == Mode.EDIT);
        updateStatus();
    }

    private void updateStatus() {
        String srcInfo = srcPane.hasImage() ? srcPane.getImageInfo() : "no source";
        String dstInfo = dstPane.hasImage() ? dstPane.getImageInfo() : "no target";
        String stage = pairs.getPlacementStage().desc();
        status.setText("Pairs: " + pairs.size() + " | Placement: " + stage + " | " + srcInfo + " | " + dstInfo);
    }

    private void openImage(Side side) {
        JFileChooser fc = new JFileChooser();
        if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
            File f = fc.getSelectedFile();
            try {
                BufferedImage bi = ImageIO.read(f);
                if (bi == null) throw new IOException("Unsupported/corrupt image.");
                if (side == Side.SOURCE) { srcPane.setImage(bi); srcPath = f; }
                else { dstPane.setImage(bi); dstPath = f; }
                updateStatus();
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, "Open failed: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    private void saveMapping() {
        if (!srcPane.hasImage() || !dstPane.hasImage()) {
            JOptionPane.showMessageDialog(this, "Open both images before saving mapping.", "Info", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        if (!pairs.isAllComplete()) {
            int res = JOptionPane.showConfirmDialog(this,
                    "There are incomplete pairs. Save only complete pairs?", "Confirm",
                    JOptionPane.OK_CANCEL_OPTION);
            if (res != JOptionPane.OK_OPTION) return;
        }

        JFileChooser fc = new JFileChooser();
        fc.setSelectedFile(new File("hemi_pairs.json"));
        if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            try (OutputStream os = new FileOutputStream(fc.getSelectedFile())) {
                String json = MappingIO.toJSON(
                        srcPath == null ? "" : srcPath.getAbsolutePath(),
                        dstPath == null ? "" : dstPath.getAbsolutePath(),
                        srcPane.getImageWidth(), srcPane.getImageHeight(),
                        dstPane.getImageWidth(), dstPane.getImageHeight(),
                        pairs.getCompletePairs());
                os.write(json.getBytes(StandardCharsets.UTF_8));
                JOptionPane.showMessageDialog(this, "Saved " + pairs.getCompletePairs().size() + " pairs.", "Saved", JOptionPane.INFORMATION_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, "Save failed: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    // --- Enums & data ---
    enum Side { SOURCE, TARGET }
    enum Mode { PLACE, EDIT }
    enum Stage { EXPECT_SRC, EXPECT_DST, READY;
        String desc(){ switch(this){
            case EXPECT_SRC: return "click SOURCE for next point";
            case EXPECT_DST: return "click TARGET to complete pair";
            default: return "place next pair";
        }}}

    static class NodePair {
        final int id;
        Point2D.Double src; // in source image pixels
        Point2D.Double dst; // in target image pixels
        NodePair(int id){ this.id=id; }
        boolean complete(){ return src!=null && dst!=null; }
        
        public NodePair(int id, Point2D.Double src, Point2D.Double dst) {
        	this.id = id;
        	this.src = src;
        	this.dst = dst;
        }
    }
    


    static class PairManager {
        private final List<NodePair> list = new ArrayList<>();
        private NodePair current = null; // building a new pair

        public void add(NodePair p) {
        	list.add(p);
        }
        
        int size(){ return list.size(); }
        List<NodePair> getAll(){ return list; }
        List<NodePair> getCompletePairs(){
            List<NodePair> out = new ArrayList<>();
            for (NodePair p: list) if (p.complete()) out.add(p);
            return out;
        }
        boolean isAllComplete(){
            for (NodePair p: list) if (!p.complete()) return false;
            return true;
        }
        Stage getPlacementStage(){
            if (current != null) return current.src==null ? Stage.EXPECT_SRC : Stage.EXPECT_DST;
            return Stage.EXPECT_SRC;
        }
        NodePair getOrCreate(){
            if (current == null || current.complete()){
                current = new NodePair(list.size()+1);
                list.add(current);
            }
            return current;
        }
        void setPoint(Side side, Point2D.Double pt){
            NodePair p = getOrCreate();
            if (side == Side.SOURCE) p.src = pt; else p.dst = pt;
            if (p.complete()) current = null;
        }
        NodePair findNearest(Side side, double xImg, double yImg, double maxDistPx){
            NodePair best=null; double bestd=maxDistPx*maxDistPx;
            for (NodePair p: list){
                Point2D.Double q = (side==Side.SOURCE)? p.src : p.dst;
                if (q==null) continue;
                double dx=q.x-xImg, dy=q.y-yImg, d=dx*dx+dy*dy;
                if (d<bestd){ best=p; bestd=d; }
            }
            return best;
        }
        void delete(NodePair p){ list.remove(p); if (current==p) current=null; }
        void clear(){ list.clear(); current=null; }
    }

    // --- Image pane with interactive nodes + magnifier ---
    static class ImagePane extends JComponent {
        private final String title;
        private final Side side;
        private final PairManager pairs;
        private BufferedImage img;
        private Mode mode = Mode.PLACE;
        private boolean showCrosshair = true;
        private boolean shiftDown = false;

        private NodePair dragPair = null;
        private boolean dragging = false;
        private Point lastMouse = null;

        private static final int NODE_RADIUS = 5;
        private static final double PICK_RADIUS_IMG = 12.0; // px in image space
        private static final int MAG_SIZE = 160; // loupe size (screen px)
        private static final int MAG_ZOOM = 4;   // 4× magnification

        ImagePane(String title, Side side, PairManager pairs){
            this.title=title; this.side=side; this.pairs=pairs;
            setOpaque(true); setBackground(new Color(245,245,248));
            setToolTipText(title + " — " + side.name());

            // Mouse
            addMouseListener(new MouseAdapter() {
                @Override public void mousePressed(MouseEvent e) {
                    if (img==null) return;
                    requestFocusInWindow();
                    Point2D.Double ip = toImageCoords(e.getPoint());
                    if (mode==Mode.PLACE && SwingUtilities.isLeftMouseButton(e)) {
                        pairs.setPoint(side, clampToImage(ip));
                        repaintSiblings();
                    } else if (mode==Mode.EDIT && SwingUtilities.isLeftMouseButton(e)) {
                        NodePair p = pairs.findNearest(side, ip.x, ip.y, PICK_RADIUS_IMG);
                        if (p != null){
                            dragPair = p; dragging=true;
                        }
                    } else if (SwingUtilities.isRightMouseButton(e)) {
                        NodePair p = pairs.findNearest(side, ip.x, ip.y, PICK_RADIUS_IMG);
                        if (p != null) {
                            int res = JOptionPane.showConfirmDialog(ImagePane.this,
                                    "Delete pair #" + p.id + "?", "Confirm", JOptionPane.OK_CANCEL_OPTION);
                            if (res==JOptionPane.OK_OPTION){ pairs.delete(p); repaintSiblings(); }
                        }
                    }
                }
                @Override public void mouseReleased(MouseEvent e) {
                    dragging=false; dragPair=null;
                }
            });
            addMouseMotionListener(new MouseMotionAdapter() {
                @Override public void mouseDragged(MouseEvent e) {
                    if (img==null || !dragging || dragPair==null) return;
                    Point2D.Double ip = clampToImage(toImageCoords(e.getPoint()));
                    if (side==Side.SOURCE) dragPair.src=ip; else dragPair.dst=ip;
                    repaintSiblings();
                }
                @Override public void mouseMoved(MouseEvent e) { lastMouse = e.getPoint(); repaint(); }
            });

            // Resize to image size
            setPreferredSize(new Dimension(900, 600));
        }

        void setMode(Mode m){ this.mode=m; repaint(); }
        void setShowCrosshair(boolean b){ this.showCrosshair=b; repaint(); }
        void setShiftDown(boolean b){ this.shiftDown=b; if (lastMouse!=null) repaint(); }
        boolean hasImage(){ return img!=null; }
        String getImageInfo(){ return title + ": " + (img==null?"-":(img.getWidth()+"×"+img.getHeight())); }
        int getImageWidth(){ return img==null?0:img.getWidth(); }
        int getImageHeight(){ return img==null?0:img.getHeight(); }

        void setImage(BufferedImage bi){
            this.img = bi;
            setPreferredSize(new Dimension(bi.getWidth(), bi.getHeight()));
            revalidate(); repaint();
        }

        private void repaintSiblings(){
            repaint();
            Container p = getParent();
            if (p!=null && p.getParent() instanceof JViewport){
                Component other = findSibling();
                if (other!=null) other.repaint();
            }
            Window w = SwingUtilities.getWindowAncestor(this);
            if (w instanceof HemiPointMapper) ((HemiPointMapper)w).updateStatus();
        }
        private Component findSibling(){
            // naive: parent of JViewport -> JSplitPane -> other component
            JSplitPane split = (JSplitPane) SwingUtilities.getAncestorOfClass(JSplitPane.class, this);
            if (split==null) return null;
            Component left = ((JViewport)((JScrollPane)split.getLeftComponent()).getViewport()).getView();
            Component right= ((JViewport)((JScrollPane)split.getRightComponent()).getViewport()).getView();
            return (left==this) ? right : left;
        }

        private Point2D.Double clampToImage(Point2D.Double p){
            double x=Math.max(0, Math.min(img.getWidth()-1, p.x));
            double y=Math.max(0, Math.min(img.getHeight()-1, p.y));
            return new Point2D.Double(x,y);
        }

        private Point2D.Double toImageCoords(Point pt){
            // panel coords == image coords (no scaling); component may be inside a scroll pane
            return new Point2D.Double(pt.x, pt.y);
        }
        private Point toPanel(Point2D.Double ip){
            return new Point((int)Math.round(ip.x), (int)Math.round(ip.y));
        }

        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2=(Graphics2D)g.create();
            g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

            if (img!=null) g2.drawImage(img, 0, 0, null);

            // title
            g2.setColor(new Color(0,0,0,90));
            g2.fillRoundRect(8,8,120,22,8,8);
            g2.setColor(Color.WHITE);
            g2.drawString(title + " (" + (mode==Mode.PLACE?"Place":"Edit") + ")", 14, 24);

            // draw nodes
            if (img!=null) {
                // hover highlight
                Point2D.Double hover = null;
                if (lastMouse!=null) {
                    Point2D.Double ip = toImageCoords(lastMouse);
                    NodePair p = pairs.findNearest(side, ip.x, ip.y, PICK_RADIUS_IMG);
                    if (p!=null) hover = (side==Side.SOURCE)? p.src : p.dst;
                }

                for (NodePair p: pairs.getAll()){
                    Point2D.Double a = (side==Side.SOURCE)? p.src : p.dst;
                    if (a==null) continue;
                    Point pt = toPanel(a);
                    // circle
                    g2.setStroke(new BasicStroke(2f));
                    g2.setColor(new Color(30,144,255));
                    g2.drawOval(pt.x-NODE_RADIUS, pt.y-NODE_RADIUS, NODE_RADIUS*2, NODE_RADIUS*2);
                    // fill hover/current
                    if (hover!=null && hover==a) {
                        g2.setColor(new Color(30,144,255,120));
                        g2.fillOval(pt.x-NODE_RADIUS, pt.y-NODE_RADIUS, NODE_RADIUS*2, NODE_RADIUS*2);
                    }
                    // id label
                    g2.setColor(new Color(0,0,0,150));
                    String s = "#"+p.id;
                    g2.drawString(s, pt.x+NODE_RADIUS+2, pt.y- NODE_RADIUS - 2);
                }
            }

            // crosshair + coords
            if (img!=null && showCrosshair && lastMouse!=null) {
                g2.setColor(new Color(0,0,0,90));
                g2.drawLine(lastMouse.x, 0, lastMouse.x, getHeight());
                g2.drawLine(0, lastMouse.y, getWidth(), lastMouse.y);
                DecimalFormat df=new DecimalFormat("0");
                g2.setColor(new Color(0,0,0,150));
                g2.fillRoundRect(lastMouse.x+10, lastMouse.y+10, 110, 20, 8,8);
                g2.setColor(Color.WHITE);
                g2.drawString(df.format(lastMouse.x)+", "+df.format(lastMouse.y), lastMouse.x+16, lastMouse.y+24);
            }

            // magnifier (loupe) when SHIFT held
            if (shiftDown && img!=null && lastMouse!=null) {
                int box = MAG_SIZE;
                int half = box/2;
                int srcBox = box / MAG_ZOOM;
                int sx = Math.max(0, Math.min(img.getWidth() - srcBox, lastMouse.x - srcBox/2));
                int sy = Math.max(0, Math.min(img.getHeight()- srcBox, lastMouse.y - srcBox/2));
                int dx = Math.min(getWidth()-box-10, lastMouse.x + 20);
                int dy = Math.min(getHeight()-box-10, lastMouse.y + 20);

                // background
                g2.setColor(new Color(0,0,0,120));
                g2.fillRoundRect(dx-4, dy-4, box+8, box+8, 10,10);
                // zoomed image
                g2.drawImage(img,
                        dx, dy, dx+box, dy+box,
                        sx, sy, sx+srcBox, sy+srcBox, null);
                // crosshair in magnifier
                g2.setColor(new Color(255,255,255,160));
                g2.drawLine(dx, dy+half, dx+box, dy+half);
                g2.drawLine(dx+half, dy, dx+half, dy+box);
                g2.setColor(Color.WHITE);
                g2.drawRoundRect(dx-4, dy-4, box+8, box+8, 10,10);
            }

            g2.dispose();
        }
    }

    // --- Mapping I/O (JSON) ---
    static class MappingIO {
        static String toJSON(String srcPath, String dstPath,
                             int sw, int sh, int dw, int dh,
                             List<NodePair> pairs) {
            StringBuilder sb = new StringBuilder(1<<16);
            sb.append("{\n");
            sb.append("  \"created_at\":\"").append(escape(Instant.now().toString())).append("\",\n");
            sb.append("  \"source\":{");
            sb.append("\"path\":\"").append(escape(srcPath)).append("\",\"width\":").append(sw).append(",\"height\":").append(sh).append("},\n");
            sb.append("  \"target\":{");
            sb.append("\"path\":\"").append(escape(dstPath)).append("\",\"width\":").append(dw).append(",\"height\":").append(dh).append("},\n");
            sb.append("  \"pairs\":[\n");
            for (int i=0;i<pairs.size();i++){
                NodePair p = pairs.get(i);
                double sx = p.src.x, sy = p.src.y, dx = p.dst.x, dy = p.dst.y;
                sb.append("    {\"id\":").append(p.id).append(", ");
                sb.append("\"src\":{\"x\":").append(fmt(sx)).append(",\"y\":").append(fmt(sy))
                  .append(",\"xn\":").append(fmt(sx/sw)).append(",\"yn\":").append(fmt(sy/sh)).append("}, ");
                sb.append("\"dst\":{\"x\":").append(fmt(dx)).append(",\"y\":").append(fmt(dy))
                  .append(",\"xn\":").append(fmt(dx/dw)).append(",\"yn\":").append(fmt(dy/dh)).append("}}");
                if (i < pairs.size()-1) sb.append(",");
                sb.append("\n");
            }
            sb.append("  ]\n");
            sb.append("}\n");
            return sb.toString();
        }
        private static String fmt(double v){ return new DecimalFormat("0.######").format(v); }
        private static String escape(String s){ return s.replace("\\","\\\\").replace("\"","\\\""); }
    }
}
