"""region.py <methods,comma>: PSNR by distance from the frame centre (normalised radius bands), both conditions."""
import os, sys, numpy as np
from PIL import Image
E=os.environ.get("FFEXP", ".")
m=np.asarray(Image.open(f"{E}/mask_1200x800.png"))>0
H,W=m.shape; yy,xx=np.mgrid[0:H,0:W]
r=np.hypot((xx+0.5-W/2)/(W/2),(yy+0.5-H/2)/(W/2))  # normalised radius
bands=[(0,0.3),(0.3,0.6),(0.6,0.9),(0.9,2)]
def rd(p): return np.asarray(Image.open(p).convert("RGB")).astype(np.float64)
for meth in sys.argv[1].split(","):
  for c in ["clean","soft"]:
    se=[np.zeros(0) for b in bands]; out=[]
    acc=[[0,0] for b in bands]
    for sc in ["a125","a130","a185","a190","fa17","f099","f1"]:
        e=((rd(f"{E}/{sc}_{c}_{meth}.png")-rd(f"{E}/{sc}_{c}_gt.png"))**2).mean(axis=2)
        for i,(a,b) in enumerate(bands):
            k=m&(r>=a)&(r<b); acc[i][0]+=e[k].sum(); acc[i][1]+=k.sum()
    print(meth,c," ".join(f"r{a}-{b}: {10*np.log10(255**2/(s/n)):.2f}" for (a,b),(s,n) in zip(bands,acc)))
