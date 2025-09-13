package defish;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import defish.HemiPointMapper.MappingIO;
import defish.HemiPointMapper.NodePair;
import defish.HemiPointMapper.PairManager;

public class PointReplicator {

	public static void main(String[] args) throws Exception {
		
		PointReplicator rep = new PointReplicator();
		rep.loadMap(new File(args[0]));
		
		PairManager completed = rep.replicate();
		
		rep.saveMapping(completed, new File(args[1]));
	}
	
	private final PairManager pairs = new PairManager();
	int idCtr = -1;
	
	public PairManager replicate() {
		PairManager dest = new PairManager();
		for(NodePair ip : pairs.getCompletePairs()) {
			dest.add(ip);
			NodePair nEast = new NodePair(idCtr, new Point2D.Double(width-ip.src.x, ip.src.y), new Point2D.Double(width-ip.dst.x, ip.dst.y));
			idCtr++;
			dest.add(nEast);
			
			NodePair sWest = new NodePair(idCtr, new Point2D.Double(ip.src.x, height-ip.src.y), new Point2D.Double(ip.dst.x, height-ip.dst.y));
			idCtr++;
			dest.add(sWest);
			
			NodePair sEast = new NodePair(idCtr, new Point2D.Double(width-ip.src.x, height-ip.src.y), new Point2D.Double(width-ip.dst.x, height-ip.dst.y));
			idCtr++;
			dest.add(sEast);

		}
		return dest;
	}
	
	public static final int width = 6000;
	public static final int height = 4000;
	public static final int hWidth = width/2;
	public static final int hHeight = height/2;
	 private void loadMap(File jsonFile) throws Exception {

	        String txt = slurp(jsonFile);

	        // (we read sizes for logging only)
	        int sw0 = width;
	        int sh0 = height;
	        int dw0 = width;
	        int dh0 = height;

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
	            if(idCtr < id) {
	            	idCtr = id;
	            }
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
	        idCtr++;
	    }
	 
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
	    
	    private void saveMapping(PairManager output, File outFile) {

	            try (OutputStream os = new FileOutputStream(outFile)) {
	                String json = MappingIO.toJSON(
	                        "sourcePath",
	                        "destPath",
	                        width, height,
	                        width, height,
	                        output.getCompletePairs());
	                os.write(json.getBytes(StandardCharsets.UTF_8));
	            } catch(Exception e) {
	            	e.printStackTrace();
	            }
	    }
}
