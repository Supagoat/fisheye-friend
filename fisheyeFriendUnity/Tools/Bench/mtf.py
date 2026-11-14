"""mtf.py out1.png [out2.png ...]: detail kept at each frequency, relative to an ideal lossless resample of the original.
Region: centre of nov11Street-1 (scale ~1.18). Prints sqrt(P_method / P_ideal) per band = effective MTF of the renderer.
FFSAMPLES = folder holding nov11Street-1.jpg; FFK = output scale (1 = same as photo, 1.19014 = Full detail).
Local scale and source position come from `Bench jac <pairs> 7952 5304 3990 2981`. Needs numpy and opencv-python."""
import os, cv2, numpy as np, sys
U=os.environ.get("FFSAMPLES", ".") + "/"
cv2.setNumThreads(2)
def lum(path, box=None):
    im=cv2.imread(path, cv2.IMREAD_COLOR)
    if box: x0,y0,x1,y1=box; im=im[y0:y1,x0:x1]
    im=im.astype(np.float64)
    return 0.114*im[...,0]+0.587*im[...,1]+0.299*im[...,2]
def psd(im):
    im=im-im.mean(); h,w=im.shape
    win=np.outer(np.hanning(h),np.hanning(w))
    F=np.abs(np.fft.fft2(im*win))**2/(win**2).sum()
    fy=np.fft.fftfreq(h); fx=np.fft.fftfreq(w)
    return F, fx, fy
# output region (same-size defish) and its source region from the Jacobian at the centre
import os
k=float(os.environ.get("FFK","1"))
N=768; cx,cy=int(round(3990*k)),int(round(2960*k))
ob=(cx-N//2,cy-N//2,cx+N//2,cy+N//2)
Jx,Jy=1.1848/k,1.1756/k; sx,sy=3992.68,3042.28-(2981-2960)*1.1756
Ms=int(N*1.18/k)
sb=(int(sx-Ms/2),int(sy-Ms/2),int(sx-Ms/2)+Ms,int(sy-Ms/2)+Ms)
src=lum(U+"nov11Street-1.jpg", sb)
Fs,fxs,fys=psd(src)
edges=np.array([0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])
# ideal output density at output frequency f: P_src(f/s) / (Jx Jy)
Ros=np.hypot(fxs[None,:]*Jx, fys[:,None]*Jy)   # source bins expressed in output frequency
def band(F,R):
    return np.array([F[(R>=a)&(R<b)].mean() for a,b in zip(edges[:-1],edges[1:])])
ideal=band(Fs,Ros)/(Jx*Jy)
print("band (cyc/px)  "+" ".join("%5.2f-%4.2f"%(a,b) for a,b in zip(edges[:-1],edges[1:])))
for p in sys.argv[1:]:
    o=lum(p, ob); F,fx,fy=psd(o); R=np.hypot(fx[None,:],fy[:,None])
    m=np.sqrt(band(F,R)/ideal)
    print("%-14s "%p.split('/')[-1][:14]+" ".join("%9.3f"%v for v in m))
