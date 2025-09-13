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

public class Hemiv2 extends JFrame {
    // --- UI ---
    private final ImgPanel imgPanel = new ImgPanel();
    private final JPanel controls = new JPanel(new GridBagLayout());
    private final JButton openBtn = new JButton("Open");
    private final JButton saveBtn = new JButton("Save");
    private final JButton presetCloister = new JButton("Preset: Cloister");

    private volatile double ky = 1.10;        // vertical strength in cylindrical space
    private volatile double horizFovDeg = 145;// cylinder half-width*2 (how much horizontal arc we see)
    private volatile double thetaMaxDeg = 90; // fisheye θ at circle edge (≈full-frame 180° fisheye)
    private volatile double rotDeg = 0;       // output rotation (deg)
    private volatile double cxPct = 0;        // fisheye center X offset (% of radius)
    private volatile double cyPct = 0;        // fisheye center Y offset (% of radius)
    private volatile double zoom = 1.00;      // manual zoom (ignored if autoFit true)
    private volatile boolean autoFit = true;  // solve largest zoom that fits circle
    private volatile int smoothCols = 16;     // small horizontal blur to tame resampling

    // 1) longitude compression grows with |y|
    private volatile double edgeCompress = 0.72; // α  (0..1)
    private volatile double edgePower = 1.20;    // p  (>=1)
    // 2) radial exponent grows with |y|
    private volatile double radialEdge = 0.36;   // γ  (0..1)
    private volatile double radialPower = 1.35;  // q  (>=1)
    // 3) gentle Y curve to avoid corner pinch
    private volatile double yCurvePow = 0.85;    // <1 expands mid/high |y| slightly
    // 4) horizon shift in cylindrical space
    private volatile double horizonShift = 0.00; // [-1..1]

    private volatile BufferedImage src, out;
    private volatile boolean rendering = false;
    
    private volatile double paniniD = 1.0;        // 0 = rectilinear, 1 ≈ cylindrical stereographic
    private volatile double squeezeStrength = 0.60;  // 0..1 (vertical squeeze vs |x|)
    private volatile double squeezePower = 2.0;      // 1.5–3.0


    public static void main(String[] args) {
        SwingUtilities.invokeLater(Hemiv2::new);
    }

    public Hemiv2() {
        super("HemiMatch — Fisheye-Hemi style remap");
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
        mb.add(file);
        setJMenuBar(mb);

        // Top bar
        JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
        top.add(openBtn); openBtn.addActionListener(e->openImage());
        top.add(saveBtn); saveBtn.addActionListener(e->saveImage());
        top.add(presetCloister); presetCloister.addActionListener(e->applyCloisterPreset());
        add(top, BorderLayout.NORTH);

        // Controls
        controls.setBorder(new EmptyBorder(8,8,8,8));
        GridBagConstraints gc = new GridBagConstraints();
        gc.gridx=0; gc.gridy=0; gc.weightx=1; gc.fill=GridBagConstraints.HORIZONTAL; gc.insets=new Insets(4,4,4,4);

        addSlider(controls, gc, "ky", 60, 180, (int)Math.round(ky*100), v->{ ky=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "horizFOV", 90, 170, (int)Math.round(horizFovDeg), v->{ horizFovDeg=v; trigger(); }, false);
        addSlider(controls, gc, "thetaMax", 70, 110, (int)Math.round(thetaMaxDeg), v->{ thetaMaxDeg=v; trigger(); }, false);
        addSlider(controls, gc, "rotate", -45, 45, (int)Math.round(rotDeg), v->{ rotDeg=v; trigger(); }, false);
        addSlider(controls, gc, "centerX (%R)", -50, 50, (int)Math.round(cxPct), v->{ cxPct=v; trigger(); }, false);
        addSlider(controls, gc, "centerY (%R)", -50, 50, (int)Math.round(cyPct), v->{ cyPct=v; trigger(); }, false);

        JCheckBox cbAuto = new JCheckBox("Auto-fit zoom", autoFit);
        cbAuto.addActionListener(e->{ autoFit=cbAuto.isSelected(); trigger(); });
        controls.add(cbAuto, gc); gc.gridy++;

        addSlider(controls, gc, "zoom (if manual)", 50, 250, (int)Math.round(zoom*100), v->{ zoom=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "smoothCols (px)", 0, 64, smoothCols, v->{ smoothCols=v; trigger(); }, false);

        // Edge/corner behavior
        JLabel sep = new JLabel("Edge & Corner Behavior"); sep.setForeground(new Color(0x224488));
        controls.add(sep, gc); gc.gridy++;
        addSlider(controls, gc, "edgeCompress α", 0, 100, (int)Math.round(edgeCompress*100), v->{ edgeCompress=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "edgePower p", 80, 200, (int)Math.round(edgePower*100), v->{ edgePower=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "radialEdge γ", 0, 100, (int)Math.round(radialEdge*100), v->{ radialEdge=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "radialPower q", 80, 200, (int)Math.round(radialPower*100), v->{ radialPower=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "yCurvePow", 70, 110, (int)Math.round(yCurvePow*100), v->{ yCurvePow=v/100.0; trigger(); }, true);
        addSlider(controls, gc, "horizonShift", -100, 100, (int)Math.round(horizonShift*100), v->{ horizonShift=v/100.0; trigger(); }, true);

        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(imgPanel), new JScrollPane(controls));
        split.setResizeWeight(1.0);
        add(split, BorderLayout.CENTER);

        setSize(1200, 800);
        setLocationRelativeTo(null);
        setVisible(true);
    }

    private void applyCloisterPreset() {
        ky = 1.10;
        horizFovDeg = 145;
        thetaMaxDeg = 90;
        rotDeg = 0;
        cxPct = 0; cyPct = 0;
        autoFit = true; zoom = 1.00; smoothCols = 16;

        edgeCompress = 0.72; edgePower = 1.20;
        radialEdge  = 0.36; radialPower = 1.35;
        yCurvePow   = 0.85; horizonShift = 0.00;
        trigger();
    }

    private void openImage() {
        JFileChooser fc = new JFileChooser();
        if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
            try {
                BufferedImage bi = ImageIO.read(fc.getSelectedFile());
                if (bi == null) throw new Exception("Unsupported/corrupt file.");
                src = ensureARGB(bi);
                out = null;
                imgPanel.setImage(src);
                trigger();
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, ex.getMessage(), "Open error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    private void saveImage() {
        if (out == null) { JOptionPane.showMessageDialog(this, "No processed image yet."); return; }
        JFileChooser fc = new JFileChooser();
        fc.setSelectedFile(new File("hemimatch.png"));
        if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            try {
                String name = fc.getSelectedFile().getName().toLowerCase();
                String fmt = (name.endsWith(".jpg")||name.endsWith(".jpeg")) ? "jpg" : "png";
                ImageIO.write(out, fmt, fc.getSelectedFile());
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, ex.getMessage(), "Save error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

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

    private BufferedImage render(BufferedImage srcImg) {
        final int W = srcImg.getWidth(), H = srcImg.getHeight();
        final double halfH = H * 0.5;                  // normalize by image half-height (stable aspect)
        final double cx = W * 0.5 + (cxPct/100.0)*halfH;
        final double cy = H * 0.5 + (cyPct/100.0)*halfH;

        // ---- Panini + fisheye params ----
        final double fovH = Math.toRadians(horizFovDeg);      // horizontal FOV of the *view* (deg slider)
        final double tanHF = Math.tan(0.5 * fovH);            // scales output x to Panini h
        final double d = paniniD;                             // NEW: Panini compression (0=perspective, 1≈cyl.stereographic)
        final double thetaMax = Math.toRadians(thetaMaxDeg);  // fisheye θ at image-circle edge (≈90° for 180° fisheye)
        final double rot = Math.toRadians(rotDeg);
        final double cR = Math.cos(rot), sR = Math.sin(rot);

        // Vertical shaping (purely vertical => verticals remain straight)
        final double kyLocal = ky;                // vertical strength (≈1.05–1.15)
        final double squeezeB = squeezeStrength;  // NEW: 0..1, more = straighter horizontals with preserved corners
        final double squeezePow = squeezePower;   // NEW: 1.5–3; how fast squeeze grows toward sides

        // choose zoom (binary search so no sampling outside the fisheye circle)
        final double zoomUse =  0;// autoFit
          //  ? solveZoomToFit_Panini(W,H,halfH,cx,cy,rotDeg,horizFovDeg,thetaMaxDeg,kyLocal,
          //                          d,squeezeB,squeezePow,0.997)
          //  : zoom;

        BufferedImage dst = new BufferedImage(W, H, BufferedImage.TYPE_INT_ARGB);

        for (int y = 0; y < H; y++) {
            // normalize by half height so x/y have same units
            double ny = (y - H*0.5) / (halfH * zoomUse);

            for (int x = 0; x < W; x++) {
                double nx = (x - W*0.5) / (halfH * zoomUse);

                // rotate output plane (just rotates the view, keeps verticality)
                double rx = nx * cR - ny * sR;
                double ry = nx * sR + ny * cR;

                // ---- Map output pixel (rx,ry) -> Panini (h,v) coordinates ----
                double h = rx * tanHF;           // horizontal coord in Panini space

                // vertical "squeeze" depends on |h| only (so verticals remain straight)
                double squeeze = 1.0 / (1.0 + squeezeB * Math.pow(Math.abs(h / tanHF), squeezePow));
                double v = (ry * tanHF) * squeeze / kyLocal;

                // ---- Invert Panini to (longitude φ, latitude λ) ----
                // Basic Panini:  h = S sinφ,  v = S tanλ,  S = (d+1)/(d + cosφ)
                // Solve cosφ from h,d (Paul Bourke derivation).
                double A = d + 1.0;
                double k = (h / A) * (h / A);
                // guard numeric issues near extremes
                double rootTerm = k*k*d*d - (k+1.0)*(k*d*d - 1.0);
                if (rootTerm < 0) rootTerm = 0;
                double cosPhi = (-k*d + Math.sqrt(rootTerm)) / (k + 1.0);
                cosPhi = Math.max(-1.0, Math.min(1.0, cosPhi));
                double S = A / (d + cosPhi);
                double sinPhi = h / S;
                sinPhi = Math.max(-1.0, Math.min(1.0, sinPhi));
                double phi = Math.atan2(sinPhi, cosPhi);          // longitude

                double lat = Math.atan( v / S );                   // latitude (cylindrical vertical, straight columns)

                // 3D ray from (phi,lat)
                double cosLat = Math.cos(lat), sinLat = Math.sin(lat);
                double Xv = cosLat * Math.sin(phi);
                double Yv = sinLat;
                double Zv = cosLat * Math.cos(phi);

                // fisheye: angle from forward axis
                double theta = Math.acos(Math.max(-1.0, Math.min(1.0, Zv)));
                double rNorm = theta / thetaMax;                   // equidistant fisheye

                // azimuth in image plane
                double imgAng = Math.atan2(Yv, Xv);
                double sX = cx + (rNorm * halfH) * Math.cos(imgAng);
                double sY = cy + (rNorm * halfH) * Math.sin(imgAng);

                int argb = sampleBilinearClamp(srcImg, sX, sY);
                dst.setRGB(x, y, argb);
            }
        }

        if (smoothCols > 0) dst = blurH(dst, Math.min(smoothCols, Math.max(1, W/2)));
        return dst;
    }



    // --- Auto-fit: find largest zoom that stays within fisheye circle ---
    private double solveZoomToFit(int W,int H,double R,double cx,double cy,double rotDeg,double fov,double ky,double horizonShift,
                                  double yPow,double eComp,double ePow,double rEdge,double rPow,double thetaMaxDeg,double safety) {
        double lo=0.60, hi=3.00;
        for(int i=0;i<22;i++){
            double mid=(lo+hi)/2.0;
            if (fits(W,H,R,cx,cy,mid,rotDeg,fov,ky,horizonShift,yPow,eComp,ePow,rEdge,rPow,thetaMaxDeg,safety)) lo=mid;
            else hi=mid;
        }
        return lo;
    }
    private boolean fits(int W,int H,double R,double cx,double cy,double zoom,double rotDeg,double fov,double ky,double horizonShift,
                         double yPow,double eComp,double ePow,double rEdge,double rPow,double thetaMaxDeg,double safety){
        final double rot=Math.toRadians(rotDeg), cR=Math.cos(rot), sR=Math.sin(rot);
        final double lonHalf=Math.toRadians(fov)*0.5;
        final double thetaMax=Math.toRadians(thetaMaxDeg);
        for(int i=0;i<720;i++){
            double t=i*(Math.PI*2/720.0);
            double u=Math.cos(t)*0.999, v=Math.sin(t)*0.999;
            double nx=u*(W*0.5)/R, ny=v*(H*0.5)/R;
            double rx=(nx/zoom)*cR - (ny/zoom)*sR;
            double ry=(nx/zoom)*sR + (ny/zoom)*cR;
            double ryc=Math.copySign(Math.pow(Math.abs(ry), yPow), ry);
            double comp=1.0 + eComp * Math.pow(Math.abs(ryc), ePow);
            double lon=(rx/comp)*lonHalf;
            double yCyl=ky*(ryc + horizonShift);
            double denom=1.0/Math.sqrt(1.0 + yCyl*yCyl);
            double Xv=Math.sin(lon)*denom, Zv=Math.cos(lon)*denom, Yv=yCyl*denom;
            double theta=Math.acos(Zv<-1?-1:(Zv>1?1:Zv));
            double exp=1.0 + rEdge * Math.pow(Math.abs(ryc), rPow);
            double rNorm=Math.pow(theta/thetaMax, exp);
            double phi=Math.atan2(Yv,Xv);
            double sx=cx + (rNorm*R)*Math.cos(phi);
            double sy=cy + (rNorm*R)*Math.sin(phi);
            double rr=Math.hypot(sx-cx, sy-cy);
            if (rr > R*safety) return false;
        }
        return true;
    }

    // --- Utilities ---
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
                out.setRGB(x,y, ((int)(sa/win)<<24)|((int)(sr/win)<<16)|((int)(sg/win)<<8)|((int)(sb/win)));
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

    // Slider helper
    private void addSlider(JPanel p, GridBagConstraints gc, String label, int min, int max, int val,
                           IntConsumer onChange, boolean showFloat){
        gc.gridy++;
        JPanel row=new JPanel(new BorderLayout(6,0));
        JLabel name=new JLabel(label);
        JLabel read=new JLabel();
        JSlider s=new JSlider(min,max,val);
        s.addChangeListener((ChangeEvent e)->{
            int v=s.getValue();
            read.setText(showFloat? new DecimalFormat("0.00").format(v/100.0) : String.valueOf(v));
            if(!s.getValueIsAdjusting()) onChange.accept(v);
        });
        read.setPreferredSize(new Dimension(64,24));
        row.add(name,BorderLayout.WEST); row.add(s,BorderLayout.CENTER); row.add(read,BorderLayout.EAST);
        read.setText(showFloat? new DecimalFormat("0.00").format(val/100.0) : String.valueOf(val));
        p.add(row, gc);
    }

    // Simple image panel
    private static class ImgPanel extends JPanel {
        private BufferedImage img;
        void setImage(BufferedImage bi){ img=bi; revalidate(); repaint();}
        @Override public Dimension getPreferredSize(){ return img==null ? new Dimension(800,600) : new Dimension(Math.max(800,img.getWidth()), Math.max(600,img.getHeight())); }
        @Override protected void paintComponent(Graphics g){
            super.paintComponent(g);
            if (img==null){ g.drawString("Open an image…", 20,20); return; }
            ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            g.drawImage(img,0,0,null);
        }
    }
}
