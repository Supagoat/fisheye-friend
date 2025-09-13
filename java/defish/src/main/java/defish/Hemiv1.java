package defish;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.DecimalFormat;
import java.util.function.DoubleConsumer;
import java.util.function.IntConsumer;

public class Hemiv1 extends JFrame {
	// UI
	private final ImagePanel imagePanel = new ImagePanel();
	private final JPanel controls = new JPanel(new GridBagLayout());
	private final JButton openBtn = new JButton("Open Image");
	private final JButton saveBtn = new JButton("Save Processed");
	private final JCheckBox normalizeY = new JCheckBox("Normalize Y", true);

	// Parameters (defaults match values you’ve liked)
	private volatile double ky = 1.20; // vertical scale in cylindrical space
	private volatile double zoom = 1.25; // overall output zoom
	private volatile double blendEdge = 0.25; // soft feather outside fisheye circle (0..1)
	private volatile double normPower = 0.75; // exponent for Y normalization
	private volatile int smoothCols = 20; // post-pass horizontal blur radius (px)

	// Extra helpers
	private volatile double horizFovDeg = 120; // total horizontal FOV mapped into frame (deg)
	private volatile double thetaMaxDeg = 90; // fisheye θ mapped to circle edge (deg)
	private volatile double rotDeg = 0; // rotate output plane before mapping (deg)
	private volatile double cxPct = 0; // lens center X offset (% of radius)
	private volatile double cyPct = 0; // lens center Y offset (% of radius)
	private volatile double horizonShift = 0; // vertical shift in cylindrical space (-1..1)

	private volatile BufferedImage src; // loaded image
	private volatile BufferedImage out; // rendered result
	private volatile boolean rendering = false; // prevent re-entrant renders

	private final JCheckBox autoFit = new JCheckBox("Auto-fit", true);
	
	private volatile double edgeCompress = 0.60; // 0..1 (0=off, 1=strong)
	private volatile double edgePower    = 1.20; // shape of compression with |y|
	private final JCheckBox expandCanvas = new JCheckBox("Expand canvas (X)", false);


	public static void main(String[] args) {
		SwingUtilities.invokeLater(Hemiv1::new);
	}

	public Hemiv1() {
		super("Fisheye-Hemi (Java) — vertical-straight / horizontal-curved defish");
		setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		setLayout(new BorderLayout());

		// Menu
		JMenuBar mb = new JMenuBar();
		JMenu file = new JMenu("File");
		JMenuItem open = new JMenuItem("Open Image…");
		JMenuItem save = new JMenuItem("Save Processed…");
		JMenuItem exit = new JMenuItem("Exit");
		open.addActionListener(e -> openImage());
		save.addActionListener(e -> saveImage());
		exit.addActionListener(e -> dispose());
		file.add(open);
		file.add(save);
		file.addSeparator();
		file.add(exit);
		mb.add(file);
		setJMenuBar(mb);

		// Top bar
		JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
		top.add(openBtn);
		top.add(saveBtn);
		openBtn.addActionListener(e -> openImage());
		saveBtn.addActionListener(e -> saveImage());

		// Controls
		controls.setBorder(new EmptyBorder(8, 8, 8, 8));
		GridBagConstraints gc = new GridBagConstraints();
		gc.gridx = 0;
		gc.gridy = 0;
		gc.weightx = 1;
		gc.fill = GridBagConstraints.HORIZONTAL;
		gc.insets = new Insets(4, 4, 4, 4);

		addSlider(controls, gc, "ky (vertical strength)", 50, 300, (int) Math.round(ky * 100), v -> {
			ky = v / 100.0;
			triggerRender();
		}, true);
		addSlider(controls, gc, "zoom", 50, 300, (int) Math.round(zoom * 100), v -> {
			zoom = v / 100.0;
			triggerRender();
		}, true);
		addSlider(controls, gc, "blendEdge", 0, 100, (int) Math.round(blendEdge * 100), v -> {
			blendEdge = v / 100.0;
			triggerRender();
		}, true);
		addSlider(controls, gc, "normPower", 25, 200, (int) Math.round(normPower * 100), v -> {
			normPower = v / 100.0;
			triggerRender();
		}, true);
		addSliderInt(controls, gc, "smoothCols (px)", 0, 100, smoothCols, v -> {
			smoothCols = v;
			triggerRender();
		});

		gc.gridy++;
		normalizeY.addActionListener(e -> {
			triggerRender();
		});
		controls.add(normalizeY, gc);
		
		addSlider(controls, gc, "edgeCompress (0..1)", 0, 100, (int)Math.round(edgeCompress*100),
		         v -> { edgeCompress = v/100.0; triggerRender(); }, true);
		addSlider(controls, gc, "edgePower", 50, 250, (int)Math.round(edgePower*100),
		         v -> { edgePower = v/100.0; triggerRender(); }, true);


		gc.gridy++;
		controls.add(autoFit, gc);
		autoFit.addActionListener(e -> triggerRender());

		addSlider(controls, gc, "horizFOV (deg)", 60, 160, (int) Math.round(horizFovDeg), v -> {
			horizFovDeg = v;
			triggerRender();
		}, false);
		addSlider(controls, gc, "thetaMax (deg)", 60, 120, (int) Math.round(thetaMaxDeg), v -> {
			thetaMaxDeg = v;
			triggerRender();
		}, false);
		addSlider(controls, gc, "rotate (deg)", -45, 45, (int) Math.round(rotDeg), v -> {
			rotDeg = v;
			triggerRender();
		}, false);
		addSlider(controls, gc, "centerX (% radius)", -50, 50, (int) Math.round(cxPct), v -> {
			cxPct = v;
			triggerRender();
		}, false);
		addSlider(controls, gc, "centerY (% radius)", -50, 50, (int) Math.round(cyPct), v -> {
			cyPct = v;
			triggerRender();
		}, false);
		addSlider(controls, gc, "horizonShift (-1..1)", -100, 100, (int) Math.round(horizonShift * 100), v -> {
			horizonShift = v / 100.0;
			triggerRender();
		}, true);

		gc.gridy++;
		controls.add(expandCanvas, gc);
		expandCanvas.addActionListener(e -> triggerRender());

		// Layout
		JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(imagePanel),
				new JScrollPane(controls));
		split.setResizeWeight(1.0);
		add(top, BorderLayout.NORTH);
		add(split, BorderLayout.CENTER);

		setSize(1200, 800);
		setLocationRelativeTo(null);
		setVisible(true);
	}

	// Helper: create labeled slider with live value readout
	private void addSlider(JPanel panel, GridBagConstraints gc, String label, int min, int max, int value,
			IntConsumer onChange, boolean showAsFloat) {
		gc.gridy++;
		JPanel row = new JPanel(new BorderLayout(6, 0));
		JLabel name = new JLabel(label);
		JLabel val = new JLabel();
		JSlider s = new JSlider(min, max, value);
		s.addChangeListener((ChangeEvent e) -> {
			int v = s.getValue();
			val.setText(showAsFloat ? fmt(v / 100.0) : String.valueOf(v));
			if (!s.getValueIsAdjusting())
				onChange.accept(v);
		});
		val.setPreferredSize(new Dimension(64, 24));
		row.add(name, BorderLayout.WEST);
		row.add(s, BorderLayout.CENTER);
		row.add(val, BorderLayout.EAST);
		val.setText(showAsFloat ? fmt(value / 100.0) : String.valueOf(value));
		panel.add(row, gc);
	}

	private void addSliderInt(JPanel panel, GridBagConstraints gc, String label, int min, int max, int value,
			IntConsumer onChange) {
		addSlider(panel, gc, label, min, max, value, onChange, false);
	}

	private static String fmt(double d) {
		return new DecimalFormat("0.00").format(d);
	}

	private void openImage() {
		JFileChooser fc = new JFileChooser();
		int r = fc.showOpenDialog(this);
		if (r == JFileChooser.APPROVE_OPTION) {
			try {
				BufferedImage img = ImageIO.read(fc.getSelectedFile());
				if (img == null)
					throw new Exception("Unsupported or corrupt image.");
				src = ensureRGBA(img);
				out = null;
				imagePanel.setImage(src);
				triggerRender();
			} catch (Exception ex) {
				JOptionPane.showMessageDialog(this, "Open failed: " + ex.getMessage(), "Error",
						JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	private void saveImage() {
		if (out == null) {
			JOptionPane.showMessageDialog(this, "No processed image yet.", "Info", JOptionPane.INFORMATION_MESSAGE);
			return;
		}
		JFileChooser fc = new JFileChooser();
		fc.setSelectedFile(new File("hemi.png"));
		int r = fc.showSaveDialog(this);
		if (r == JFileChooser.APPROVE_OPTION) {
			try {
				String name = fc.getSelectedFile().getName().toLowerCase();
				String fmt = name.endsWith(".jpg") || name.endsWith(".jpeg") ? "png" : "png";
				ImageIO.write(out, fmt, fc.getSelectedFile());
			} catch (Exception ex) {
				JOptionPane.showMessageDialog(this, "Save failed: " + ex.getMessage(), "Error",
						JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	private void triggerRender() {
		if (src == null)
			return;
		if (rendering)
			return;
		rendering = true;

		SwingWorker<BufferedImage, Void> worker = new SwingWorker<>() {
			@Override
			protected BufferedImage doInBackground() {
				return render(src);
			}

			@Override
			protected void done() {
				try {
					out = get();
					imagePanel.setImage(out);
				} catch (Exception ignored) {
				}
				rendering = false;
			}
		};
		worker.execute();
	}

	private BufferedImage render(BufferedImage source) {
		final int W = source.getWidth();
		final int H = source.getHeight();
		final int D = Math.min(W, H);
		final double radius = D * 0.5;
		final double cx = W * 0.5 + (cxPct / 100.0) * radius;
		final double cy = H * 0.5 + (cyPct / 100.0) * radius;

		// Precompute constants
		final double rot = Math.toRadians(rotDeg);
		final double cosR = Math.cos(rot), sinR = Math.sin(rot);
		final double lonHalfRange = Math.toRadians(horizFovDeg) * 0.5; // cylinder half-width in radians
		final double thetaMax = Math.toRadians(thetaMaxDeg); // fisheye θ at circle edge
		final double invThetaMax = 1.0 / thetaMax;

		int outW = expandCanvas.isSelected() ? (int)Math.round(W * 1.25) : W; // 25% wider
		int outH = H;
		BufferedImage dst = new BufferedImage(outW, outH, BufferedImage.TYPE_INT_ARGB);

		// For each output pixel, compute inverse map to source
		for (int y = 0; y < H; y++) {
			double ny = (y - H * 0.5) / (radius); // [-~1..~1] space
			double nny = ny / zoom; // apply zoom

			for (int x = 0; x < W; x++) {
				double nx = (x - W * 0.5) / (radius);
				double nnx = nx / zoom;

				// Rotate output plane if desired
				double rx = nnx * cosR - nny * sinR;
				double ry = nnx * sinR + nny * cosR;

				// Optional Y normalization (soften vertical squeeze near top/bottom)
				if (normalizeY.isSelected()) {
					ry = Math.signum(ry) * Math.pow(Math.abs(ry), normPower);
				}

				// Apply vertical scale and optional horizon shift in cylindrical space
				double yCyl = ky * (ry + horizonShift);

				// BEFORE: double lon = rx * lonHalfRange;
				// AFTER: compress lon more as |ry| grows (Panini-like)
				double c = 1.0 + edgeCompress * Math.pow(Math.abs(ry), edgePower);
				double lon = (rx / c) * lonHalfRange;

				
				// Map to 3D direction via inverse cylindrical projection
				// lon in [-lonHalfRange, +lonHalfRange] when rx in [-1..1]
				//double lon = rx * lonHalfRange;
				double denom = 1.0 / Math.sqrt(1.0 + yCyl * yCyl);

				double Xv = Math.sin(lon) * denom;
				double Zv = Math.cos(lon) * denom;
				double Yv = yCyl * denom;

				// Convert 3D ray to fisheye image sample (equidistant-ish): r = f * theta
				// theta from forward axis (Z)
				double theta = Math.acos(clamp(Zv, -1.0, 1.0)); // [0..π]
				double rNorm = theta * invThetaMax; // 0 at center, 1 at circle edge
				// azimuth around forward axis (image plane angle)
				double phi = Math.atan2(Yv, Xv);

				// Source pixel
				double srcX = cx + (rNorm * radius) * Math.cos(phi);
				double srcY = cy + (rNorm * radius) * Math.sin(phi);

				// Feather if outside the fisheye circle (blendEdge)
				// Circle (center cx,cy, radius R)
				double dx = srcX - cx, dy = srcY - cy;
				double rr = Math.sqrt(dx * dx + dy * dy);
				double alpha = 1.0;
				if (rr > radius) {
					double feather = radius * Math.max(1e-6, blendEdge);
					if (rr > radius + feather) {
						alpha = 0.0;
					} else {
						alpha = 1.0 - (rr - radius) / feather;
					}
				}

				int argb = sampleBilinearClamp(source, srcX, srcY);
				if (alpha < 1.0) {
					int a = (argb >>> 24) & 0xFF;
					a = (int) Math.round(a * alpha);
					argb = (a << 24) | (argb & 0x00FFFFFF);
				}
				dst.setRGB(x, y, argb);
			}
		}

		if (smoothCols > 0) {
			dst = horizontalBoxBlur(dst, Math.min(smoothCols, Math.max(1, W / 2)));
		}

		double zoomLocal = zoom;
		if (autoFit.isSelected()) {
			zoomLocal = solveZoomToFit(W, H, radius, cx, cy, rotDeg, horizFovDeg, ky, horizonShift,
					normalizeY.isSelected(), normPower, thetaMaxDeg, 0.995 /* safety */);
		}

		return dst;
	}

	private double solveZoomToFit(int W, int H, double radius, double cx, double cy, double rotDeg, double horizFovDeg,
			double ky, double horizonShift, boolean normY, double normPower, double thetaMaxDeg, double safety) {
		double lo = 0.6, hi = 3.0; // practical bounds
		for (int it = 0; it < 22; it++) { // ~1e-6 precision
			double mid = 0.5 * (lo + hi);
			if (borderFits(W, H, radius, cx, cy, mid, rotDeg, horizFovDeg, ky, horizonShift, normY, normPower,
					thetaMaxDeg, safety)) {
				lo = mid;
			} else {
				hi = mid;
			}
		}
		return lo;
	}

	private boolean borderFits(int W, int H, double R, double cx, double cy, double zoom, double rotDeg,
			double horizFovDeg, double ky, double horizonShift, boolean normY, double normPower, double thetaMaxDeg,
			double safety) {
		final double rot = Math.toRadians(rotDeg);
		final double cosR = Math.cos(rot), sinR = Math.sin(rot);
		final double lonHalfRange = Math.toRadians(horizFovDeg) * 0.5;
		final double thetaMax = Math.toRadians(thetaMaxDeg);
		final double invThetaMax = 1.0 / thetaMax;

// sample a ring of output pixels near the border
		for (int i = 0; i < 720; i++) {
			double t = (i / 720.0) * Math.PI * 2.0;
// parametrize a rectangle border smoothly
			double u = (Math.cos(t) * 0.999); // [-1..1]
			double v = (Math.sin(t) * 0.999);
// map to the *outer rectangle* edges
			double nx = u * (W * 0.5) / R;
			double ny = v * (H * 0.5) / R;

// inverse of the same pipeline in render()
			double nnx = nx / zoom, nny = ny / zoom;
			double rx = nnx * cosR - nny * sinR;
			double ry = nnx * sinR + nny * cosR;
			if (normY)
				ry = Math.signum(ry) * Math.pow(Math.abs(ry), normPower);
			double yCyl = ky * (ry + horizonShift);

			double lon = rx * lonHalfRange;
			double denom = 1.0 / Math.sqrt(1.0 + yCyl * yCyl);
			double Xv = Math.sin(lon) * denom;
			double Zv = Math.cos(lon) * denom;
			double Yv = yCyl * denom;

			double theta = Math.acos(Math.max(-1, Math.min(1, Zv)));
			double rNorm = theta * invThetaMax;
			double phi = Math.atan2(Yv, Xv);

			double srcX = cx + (rNorm * R) * Math.cos(phi);
			double srcY = cy + (rNorm * R) * Math.sin(phi);
			double rr = Math.hypot(srcX - cx, srcY - cy);

			if (rr > R * safety)
				return false; // would sample outside circle
		}
		return true;
	}

	private static double clamp(double v, double lo, double hi) {
		return v < lo ? lo : (v > hi ? hi : v);
	}

	private static BufferedImage ensureRGBA(BufferedImage in) {
		if (in.getType() == BufferedImage.TYPE_INT_ARGB)
			return in;
		BufferedImage out = new BufferedImage(in.getWidth(), in.getHeight(), BufferedImage.TYPE_INT_ARGB);
		Graphics2D g = out.createGraphics();
		g.drawImage(in, 0, 0, null);
		g.dispose();
		return out;
	}

	// Bilinear sampler with edge clamp
	private static int sampleBilinearClamp(BufferedImage img, double x, double y) {
		int w = img.getWidth(), h = img.getHeight();
		if (x < 0)
			x = 0;
		if (y < 0)
			y = 0;
		if (x > w - 1)
			x = w - 1;
		if (y > h - 1)
			y = h - 1;

		int x0 = (int) Math.floor(x), y0 = (int) Math.floor(y);
		int x1 = Math.min(x0 + 1, w - 1), y1 = Math.min(y0 + 1, h - 1);
		double fx = x - x0, fy = y - y0;

		int c00 = img.getRGB(x0, y0);
		int c10 = img.getRGB(x1, y0);
		int c01 = img.getRGB(x0, y1);
		int c11 = img.getRGB(x1, y1);

		int a = bilerp((c00 >>> 24) & 0xFF, (c10 >>> 24) & 0xFF, (c01 >>> 24) & 0xFF, (c11 >>> 24) & 0xFF, fx, fy);
		int r = bilerp((c00 >>> 16) & 0xFF, (c10 >>> 16) & 0xFF, (c01 >>> 16) & 0xFF, (c11 >>> 16) & 0xFF, fx, fy);
		int g = bilerp((c00 >>> 8) & 0xFF, (c10 >>> 8) & 0xFF, (c01 >>> 8) & 0xFF, (c11 >>> 8) & 0xFF, fx, fy);
		int b = bilerp(c00 & 0xFF, c10 & 0xFF, c01 & 0xFF, c11 & 0xFF, fx, fy);

		return (a << 24) | (r << 16) | (g << 8) | b;
	}

	private static int bilerp(int c00, int c10, int c01, int c11, double fx, double fy) {
		double i0 = c00 + (c10 - c00) * fx;
		double i1 = c01 + (c11 - c01) * fx;
		return (int) Math.round(i0 + (i1 - i0) * fy);
	}

	// Fast horizontal box blur on ARGB
	private static BufferedImage horizontalBoxBlur(BufferedImage src, int radius) {
		int w = src.getWidth(), h = src.getHeight();
		BufferedImage out = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
		int win = radius * 2 + 1;

		for (int y = 0; y < h; y++) {
			long sa = 0, sr = 0, sg = 0, sb = 0;

			// Initialize window
			for (int x = -radius; x <= radius; x++) {
				int xx = clampi(x, 0, w - 1);
				int c = src.getRGB(xx, y);
				sa += (c >>> 24) & 0xFF;
				sr += (c >>> 16) & 0xFF;
				sg += (c >>> 8) & 0xFF;
				sb += c & 0xFF;
			}

			for (int x = 0; x < w; x++) {
				int a = (int) (sa / win);
				int r = (int) (sr / win);
				int g = (int) (sg / win);
				int b = (int) (sb / win);
				out.setRGB(x, y, (a << 24) | (r << 16) | (g << 8) | b);

				// Slide window
				int xOut = x - radius;
				int xIn = x + radius + 1;
				int cOut = src.getRGB(clampi(xOut, 0, w - 1), y);
				int cIn = src.getRGB(clampi(xIn, 0, w - 1), y);
				sa += ((cIn >>> 24) & 0xFF) - ((cOut >>> 24) & 0xFF);
				sr += ((cIn >>> 16) & 0xFF) - ((cOut >>> 16) & 0xFF);
				sg += ((cIn >>> 8) & 0xFF) - ((cOut >>> 8) & 0xFF);
				sb += (cIn & 0xFF) - (cOut & 0xFF);
			}
		}
		return out;
	}

	private static int clampi(int v, int lo, int hi) {
		return v < lo ? lo : (v > hi ? hi : v);
	}

	// Simple image panel that fits the image to the viewport while preserving scale
	private static class ImagePanel extends JPanel {
		private BufferedImage img;

		public void setImage(BufferedImage i) {
			this.img = i;
			revalidate();
			repaint();
		}

		@Override
		public Dimension getPreferredSize() {
			if (img == null)
				return new Dimension(800, 600);
			return new Dimension(Math.max(800, img.getWidth()), Math.max(600, img.getHeight()));
		}

		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (img == null) {
				g.setColor(Color.DARK_GRAY);
				g.drawString("Open an image to start…", 20, 20);
				return;
			}
			Graphics2D g2 = (Graphics2D) g.create();
			g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
			g2.drawImage(img, 0, 0, null);
			g2.dispose();
		}
	}
}
