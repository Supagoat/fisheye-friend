package defish;


import com.google.gson.*;
import com.google.gson.annotations.SerializedName;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.awt.*;

public class GridDraw {
	
	public static int w = 6000;
	public static int h = 4000;
	
	public static void main(String[] args) throws Exception{
		BufferedImage i = new BufferedImage(w,h, BufferedImage.TYPE_INT_ARGB);

		
		for(int y=0;y<h;y++) {
			for(int x=0;x<w;x++) {
				i.setRGB(x,y,Color.WHITE.getRGB());
			}
		}
		int red = 10;
		int green = 60;
		int blue = 60;
		
		
		
		

		red = 200;
		green = 200;
		blue = 160;
		
		for(int x=0;x<w;x++) {
			for(int y=0;y<50;y++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			blue+=30;
		}
		
		blue = 160;
		for(int x=0;x<w;x++) {
			for(int y=h-50;y<h;y++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			blue+=30;
		}
		blue = 160;
		for(int x=0;x<50;x++) {
			for(int y=0;y<h;y++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			blue+=30;
		}
		
		blue = 160;
		for(int x=w-50;x<w;x++) {
			for(int y=0;y<h;y++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			blue+=30;
		}
		
		
		
		red = 10;
		green = 60;
		blue = 60;
		// the reds
		for(int y=0;y<h;y+=25) {
			for(int x=0;x<w;x++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			red+=60;
			if(red > 255) {
				red = 10;
			}
		}
		
		red = 40;
		green = 10;
		blue = 60;
		
		for(int x=0;x<w;x+=25) {
			for(int y=0;y<h;y++) {
				int color = 0xff000000;
				color = color | red << 16;
				color = color | green << 8;
				color = color | blue;
				i.setRGB(x,y,color);

			}
			green+=60;
			if(green > 255) {
				green = 10;
			}
		}
		
		
		
		ImageIO.write(i, "png", new File(args[0]));
	}
}