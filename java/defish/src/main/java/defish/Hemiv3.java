package defish;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.DecimalFormat;
import java.util.function.IntConsumer;

public class Hemiv3 extends JFrame {
    // ===== Image + UI =====
    private final ImgPanel imgPanel = new ImgPanel();
    private final JPanel controls = new JPanel(new GridBagLayout());
    private final JButton openBtn = new JButton("Open");
    private final JButton saveBtn = new JButton("Save");
    private final JButton presetCloister = new JButton("Preset: Cloister");

    // ===== Core view / lens parameters =====
    private volatile double horizFOVDeg = 148;   // horizontal FOV of the view (deg)
    private volatile double thetaMaxDeg  = 90;   // fisheye theta at image-circle edge (~90 for 180° FF fisheye)
    private volatile double rotDeg       = 0;    // rotate the view (deg)
    private volatile double ky           = 1.10; // vertical scaling inside the cylinder (1.05–1.15 good)

    // ===== vertical squeeze =====
    private volatile double paniniD         = 1.10; // 0 = rectilinear, 1 ≈ cylindrical stereographic
    private volatile double squeezeStrength = 0.62; // 0..1, vertical squeeze vs |x| (preserves corners)
    private volatile double squeezePower    = 2.1;  // 1.5–3.0, how quickly squeeze grows toward sides

    // ===== Framing / sampling =====
    private volatile boolean autoFit = true; // binary-search zoom so nothing samples outside fisheye circle
    private volatile double  zoom    = 1.00; // used only if autoFit=false
    private volatile int     smoothCols = 16; // small horizontal blur to reduce resampling chatter

    // fisheye image center offsets (percent of radius = half height)
    private volatile double cxPct = 0.0, cyPct = 0.0;

    // Images
    private volatile BufferedImage src, out;
    private volatile boolean rendering = false;

    public static void main(String[] args) {
        SwingUtilities.invokeLater(Hemiv3::new);
    }

    public Hemiv3() {
        super("HemiPanini — vertical-straight / curved-horiz defish (Panini-Strict)");
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Menu
        JMenuBar mb = new JMenuBar();
        JMenu file = new JMenu("File");
        JMenuItem mOpen = new JMenuItem("Open…");
        JMenuItem mSave = new JMenuItem("Save…");
        mOpen.addActionListener(e -> openImage());
        mSave.addActionListener(e -> saveImage());
        file.add(mOpen); file.add(mSave);
        mb.add(file); setJMenuBar(mb);

        // Top bar
        JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
        top.add(openBtn); openBtn.addActionListener(e -> openImage());
        top.add(saveBtn); saveBtn.addActionListener(e -> saveImage());
        top.add(presetCloister); presetCloister.addActionListener(e -> applyCloisterPreset());
        add(top, BorderLayout.NORTH);

        // Controls
        controls.setBorder(new EmptyBorder(8,8,8,8));
        GridBagConstraints gc = new GridBagConstraints();
        gc.gridx=0; gc.gridy=0; gc.weightx=1; gc.fill=GridBagConstraints.HORIZONTAL; gc.insets=new Insets(4,4,4,4);

        addSlider(controls, gc, "horizFOV°", 90, 170, (int)Math.round(horizFOVDeg),
                v -> { horizFOVDeg = v; trigger(); }, false);
        addSlider(controls, gc, "thetaMax°", 70, 110, (int)Math.round(thetaMaxDeg),
                v -> { thetaMaxDeg = v; trigger(); }, false);
        addSlider(controls, gc, "rotate°", -45, 45, (int)Math.round(rotDeg),
                v -> { rotDeg = v; trigger(); }, false);
        addSlider(controls, gc, "ky", 80, 150, (int)Math.round(ky*100),
                v -> { ky = v/100.0; trigger(); }, true);

        JLabel sep1 = new JLabel("Panini + vertical squeeze"); sep1.setForeground(new Color(0x224488));
        controls.add(sep1, gc); gc.gridy++;

        addSlider(controls, gc, "paniniD", 0, 200, (int)Math.round(paniniD*100),
                v -> { paniniD = v/100.0; trigger(); }, true);
        addSlider(controls, gc, "squeezeStrength", 0, 100, (int)Math.round(squeezeStrength*100),
                v -> { squeezeStrength = v/100.0; trigger(); }, true);
        addSlider(controls, gc, "squeezePower", 100, 350, (int)Math.round(squeezePower*100),
                v -> { squeezePower = v/100.0; trigger(); }, true);

        JLabel sep2 = new JLabel("Framing / sampling"); sep2.setForeground(new Color(0x224488));
        controls.add(sep2, gc); gc.gridy++;

        JCheckBox cbAuto = new JCheckBox("Auto-fit zoom", autoFit);
        cbAuto.addActionListener(e -> { autoFit = cbAuto.isSelected(); trigger(); });
        controls.add(cbAuto, gc); gc.gridy++;

        addSlider(controls, gc, "zoom (if manual)", 50, 250, (int)Math.round(zoom*100),
                v -> { zoom = v/100.0; trigger(); }, true);
        addSlider(controls, gc, "smoothCols (px)", 0, 64, smoothCols,
                v -> { smoothCols = v; trigger(); }, false);
        addSlider(controls, gc, "centerX (%R)", -50, 50, (int)Math.round(cxPct),
                v -> { cxPct = v; trigger(); }, false);
        addSlider(controls, gc, "centerY (%R)", -50, 50, (int)Math.round(cyPct),
                v -> { cyPct = v; trigger(); }, false);

        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(imgPanel), new JScrollPane(controls));
        split.setResizeWeight(1.0);
        add(split, BorderLayout.CENTER);

        setSize(1200, 800);
        setLocationRelativeTo(null);
        setVisible(true);
    }

    private void applyCloisterPreset() {
        horizFOVDeg = 148;
        thetaMaxDeg  = 90;
        rotDeg       = 0;
        ky           = 1.10;

        paniniD         = 1.10;
        squeezeStrength = 0.62;
        squeezePower    = 2.10;

        autoFit = true; zoom = 1.00; smoothCols = 16;
        cxPct = 0; cyPct = 0;
        trigger();
    }

    // ===== File I/O =====
    private void openImage() {
        JFileChooser fc = new JFileChooser();
        if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
            try {
                BufferedImage bi = ImageIO.read(fc.getSelectedFile());
                if (bi == null) throw new Exception("Unsupported/corrupt file.");
                src = ensureARGB(bi);
                out = null; imgPanel.setImage(src);
                trigger();
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, ex.getMessage(), "Open error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }
    private void saveImage() {
        if (out == null) { JOptionPane.showMessageDialog(this, "No processed image yet."); return; }
        JFileChooser fc = new JFileChooser();
        fc.setSelectedFile(new File("hemipanini.png"));
        if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            try {
                String name = fc.getSelectedFile().getName().toLowerCase();
                String fmt = (name.endsWith(".jpg") || name.endsWith(".jpeg")) ? "jpg" : "png";
                ImageIO.write(out, fmt, fc.getSelectedFile());
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, ex.getMessage(), "Save error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    // ===== Render orchestration =====
    private void trigger() {
        if (src == null || rendering) return;
        rendering = true;
        SwingWorker<BufferedImage,Void> w = new SwingWorker<>() {
            @Override protected BufferedImage doInBackground() { return render(src); }
            @Override protected void done() {
                try { out = get(); imgPanel.setImage(out); } catch (Exception ignored) {}
                rendering = false;
            }
        };
        w.execute();
    }

    // ===== Core render: Panini-Strict + vertical-only squeeze, then inverse fisheye =====
    private BufferedImage render(BufferedImage srcImg) {
        final int W = srcImg.getWidth(), H = srcImg.getHeight();
        final double halfH = H * 0.5;                 // use half height as fisheye radius
        final double R = halfH;

        final double cx = W * 0.5 + (cxPct/100.0)*R;  // fisheye center offsets
        final double cy = H * 0.5 + (cyPct/100.0)*R;

        final double fovH = Math.toRadians(horizFOVDeg);
        final double tanHF = Math.tan(0.5 * fovH);
        final double d = paniniD;
        final double thetaMax = Math.toRadians(thetaMaxDeg);

        final double rot = Math.toRadians(rotDeg);
        final double cR = Math.cos(rot), sR = Math.sin(rot);

        // Auto-fit zoom solves for largest zoom that still stays inside circle
        final double zoomUse = autoFit
                ? solveZoomToFit_Panini(W, H, R, cx, cy, rotDeg, horizFOVDeg, thetaMaxDeg,
                ky, d, squeezeStrength, squeezePower, 0.997)
                : zoom;

        BufferedImage dst = new BufferedImage(W, H, BufferedImage.TYPE_INT_ARGB);

        for (int y = 0; y < H; y++) {
            double ny = (y - H*0.5) / (R * zoomUse);    // normalized output y

            for (int x = 0; x < W; x++) {
                double nx = (x - W*0.5) / (R * zoomUse);

                // rotate output plane
                double rx = nx * cR - ny * sR;
                double ry = nx * sR + ny * cR;

                // --- Panini forward coords (h,v) ---
                // horizontal coordinate in Panini space
                double h = rx * tanHF;

                // vertical squeeze depends on |x| only -> verticals remain straight
                double squeeze = 1.0 / (1.0 + squeezeStrength * Math.pow(Math.abs(rx), squeezePower));
                double v = (ry * tanHF) * squeeze / ky;

                // --- Invert Panini to (phi, lat) ---
                // h = S sinφ, v = S tanλ,  S = (d+1)/(d + cosφ)
                final double A = d + 1.0;
                double k = (h / A) * (h / A);
                // Quadratic in cosφ (stable form)
                double root = k*k*d*d - (k + 1.0) * (k*d*d - 1.0);
                if (root < 0) root = 0;
                double cosPhi = (-k*d + Math.sqrt(root)) / (k + 1.0);
                cosPhi = clamp(cosPhi, -1.0, 1.0);
                double S = A / (d + cosPhi);
                double sinPhi = clamp(h / S, -1.0, 1.0);
                double phi = Math.atan2(sinPhi, cosPhi);

                double lat = Math.atan( v / S );

                // 3D ray
                double cosLat = Math.cos(lat), sinLat = Math.sin(lat);
                double Xv = cosLat * Math.sin(phi);
                double Yv = sinLat;
                double Zv = cosLat * Math.cos(phi);

                // angle from forward axis
                double theta = Math.acos(clamp(Zv, -1.0, 1.0));
                double rNorm = theta / thetaMax; // equidistant fisheye

                // azimuth in image plane
                double imgAng = Math.atan2(Yv, Xv);

                // back to fisheye pixel
                double sX = cx + (rNorm * R) * Math.cos(imgAng);
                double sY = cy + (rNorm * R) * Math.sin(imgAng);

                int argb = sampleBilinearClamp(srcImg, sX, sY);
                dst.setRGB(x, y, argb);
            }
        }

        if (smoothCols > 0) dst = blurH(dst, Math.min(smoothCols, Math.max(1, W/2)));
        return dst;
    }

    // ===== Auto-fit zoom for Panini mapping =====
    private double solveZoomToFit_Panini(int W,int H,double R,double cx,double cy,double rotDeg,
                                         double fovHdeg,double thetaMaxDeg,double ky,double d,
                                         double squeezeB,double squeezePow,double safety){
        double lo=0.60, hi=3.00;
        for(int i=0;i<22;i++){
            double mid=(lo+hi)/2.0;
            if (fitsPanini(W,H,R,cx,cy,mid,rotDeg,fovHdeg,thetaMaxDeg,ky,d,squeezeB,squeezePow,safety)) lo=mid; else hi=mid;
        }
        return lo;
    }

    private boolean fitsPanini(int W,int H,double R,double cx,double cy,double zoom,double rotDeg,
                               double fovHdeg,double thetaMaxDeg,double ky,double d,
                               double squeezeB,double squeezePow,double safety){
        final double rot=Math.toRadians(rotDeg), cR=Math.cos(rot), sR=Math.sin(rot);
        final double fovH=Math.toRadians(fovHdeg), tanHF=Math.tan(0.5*fovH);
        final double A=d+1.0;
        final double thetaMax=Math.toRadians(thetaMaxDeg);

        int N=720;
        for(int i=0;i<N;i++){
            double t = (i/(double)N)*Math.PI*2.0;
            // rectangle border sample in normalized output coords
            double u = Math.cos(t)*0.999, v = Math.sin(t)*0.999;
            double nx = u * (W*0.5) / (R*zoom);
            double ny = v * (H*0.5) / (R*zoom);

            double rx = nx * cR - ny * sR;
            double ry = nx * sR + ny * cR;

            double h = rx * tanHF;
            double squeeze = 1.0 / (1.0 + squeezeB * Math.pow(Math.abs(rx), squeezePow));
            double vpan = (ry * tanHF) * squeeze / ky;

            double k = (h / A) * (h / A);
            double root = k*k*d*d - (k+1.0)*(k*d*d - 1.0);
            if (root < 0) root = 0;
            double cosPhi = (-k*d + Math.sqrt(root)) / (k + 1.0);
            cosPhi = clamp(cosPhi, -1.0, 1.0);
            double S = A / (d + cosPhi);
            double sinPhi = clamp(h / S, -1.0, 1.0);
            double phi = Math.atan2(sinPhi, cosPhi);
            double lat = Math.atan(vpan / S);

            double cosLat = Math.cos(lat), sinLat = Math.sin(lat);
            double Xv = cosLat * Math.sin(phi);
            double Yv = sinLat;
            double Zv = cosLat * Math.cos(phi);
            double theta = Math.acos(clamp(Zv, -1.0, 1.0));
            double rNorm = theta / thetaMax;
            double imgAng = Math.atan2(Yv, Xv);
            double sX = cx + (rNorm * R) * Math.cos(imgAng);
            double sY = cy + (rNorm * R) * Math.sin(imgAng);

            double rr = Math.hypot(sX - cx, sY - cy);
            if (rr > R * safety) return false;
        }
        return true;
    }

    // ===== Utilities =====
    private static double clamp(double v,double lo,double hi){ return v<lo?lo:(v>hi?hi:v); }

    private static BufferedImage ensureARGB(BufferedImage in){
        if (in.getType()==BufferedImage.TYPE_INT_ARGB) return in;
        BufferedImage out=new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
        Graphics2D g=out.createGraphics(); g.drawImage(in,0,0,null); g.dispose(); return out;
    }

    private static int sampleBilinearClamp(BufferedImage img,double x,double y){
        int w=img.getWidth(), h=img.getHeight();
        if (x<0) x=0; if (y<0) y=0; if (x>w-1) x=w-1; if (y>h-1) y=h-1;
        int x0=(int)Math.floor(x), y0=(int)Math.floor(y);
        int x1=Math.min(x0+1,w-1), y1=Math.min(y0+1,h-1);
        double fx=x-x0, fy=y-y0;
        int c00=img.getRGB(x0,y0), c10=img.getRGB(x1,y0), c01=img.getRGB(x0,y1), c11=img.getRGB(x1,y1);
        int a=lerp2((c00>>>24)&255,(c10>>>24)&255,(c01>>>24)&255,(c11>>>24)&255,fx,fy);
        int r=lerp2((c00>>>16)&255,(c10>>>16)&255,(c01>>>16)&255,(c11>>>16)&255,fx,fy);
        int g=lerp2((c00>>>8)&255,(c10>>>8)&255,(c01>>>8)&255,(c11>>>8)&255,fx,fy);
        int b=lerp2(c00&255,c10&255,c01&255,c11&255,fx,fy);
        return (a<<24)|(r<<16)|(g<<8)|b;
    }
    private static int lerp2(int c00,int c10,int c01,int c11,double fx,double fy){
        double i0=c00+(c10-c00)*fx, i1=c01+(c11-c01)*fx;
        return (int)Math.round(i0+(i1-i0)*fy);
    }

    private static BufferedImage blurH(BufferedImage src,int radius){
        if (radius<=0) return src;
        int w=src.getWidth(), h=src.getHeight();
        BufferedImage out=new BufferedImage(w,h,BufferedImage.TYPE_INT_ARGB);
        int win=radius*2+1;
        for(int y=0;y<h;y++){
            long sa=0,sr=0,sg=0,sb=0;
            for(int x=-radius;x<=radius;x++){
                int xx=Math.max(0,Math.min(w-1,x));
                int c=src.getRGB(xx,y);
                sa+=(c>>>24)&255; sr+=(c>>>16)&255; sg+=(c>>>8)&255; sb+=c&255;
            }
            for(int x=0;x<w;x++){
                int a=(int)(sa/win), r=(int)(sr/win), g=(int)(sg/win), b=(int)(sb/win);
                out.setRGB(x,y,(a<<24)|(r<<16)|(g<<8)|b);
                int xOut=x-radius, xIn=x+radius+1;
                int cOut=src.getRGB(Math.max(0,Math.min(w-1,xOut)),y);
                int cIn =src.getRGB(Math.max(0,Math.min(w-1,xIn )),y);
                sa+=((cIn>>>24)&255)-((cOut>>>24)&255);
                sr+=((cIn>>>16)&255)-((cOut>>>16)&255);
                sg+=((cIn>>>8)&255)-((cOut>>>8)&255);
                sb+=( cIn&255)-( cOut&255);
            }
        }
        return out;
    }

    // ===== Slider helper =====
    private void addSlider(JPanel p, GridBagConstraints gc, String label, int min, int max, int val,
                           IntConsumer onChange, boolean showFloat){
        gc.gridy++;
        JPanel row=new JPanel(new BorderLayout(6,0));
        JLabel name=new JLabel(label);
        JLabel read=new JLabel();
        JSlider s=new JSlider(min,max,val);
        s.addChangeListener((ChangeEvent e)->{
            int v=s.getValue();
            read.setText(showFloat ? new DecimalFormat("0.00").format(v/100.0) : String.valueOf(v));
            if (!s.getValueIsAdjusting()) onChange.accept(v);
        });
        read.setPreferredSize(new Dimension(68,24));
        row.add(name,BorderLayout.WEST); row.add(s,BorderLayout.CENTER); row.add(read,BorderLayout.EAST);
        read.setText(showFloat ? new DecimalFormat("0.00").format(val/100.0) : String.valueOf(val));
        p.add(row, gc);
    }

    // ===== Simple image panel =====
    private static class ImgPanel extends JPanel {
        private BufferedImage img;
        void setImage(BufferedImage bi){ img=bi; revalidate(); repaint(); }
        @Override public Dimension getPreferredSize(){
            return img==null ? new Dimension(900,600)
                             : new Dimension(Math.max(900,img.getWidth()), Math.max(600,img.getHeight()));
        }
        @Override protected void paintComponent(Graphics g){
            super.paintComponent(g);
            if (img==null){ g.setColor(Color.DARK_GRAY); g.drawString("Open an image…", 20,20); return; }
            Graphics2D g2=(Graphics2D)g.create();
            g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            g2.drawImage(img,0,0,null);
            g2.dispose();
        }
    }
}
