"""score.py <methods,comma> [scenes] [conds]: PSNR vs ground truth over the mask.

Run in the folder holding the outputs of `Bench exp-prep` / `exp-run` (or set FFEXP to it). Needs numpy and Pillow.
"""
import os, sys, numpy as np
from PIL import Image
E = os.environ.get("FFEXP", "."); methods = sys.argv[1].split(",")
scenes = sys.argv[2].split(",") if len(sys.argv) > 2 and sys.argv[2] else ["a125", "a130", "a185", "a190", "fa17", "f099", "f1"]
conds = sys.argv[3].split(",") if len(sys.argv) > 3 else ["clean", "soft"]
def rd(p): return np.asarray(Image.open(p).convert("RGB")).astype(np.float64)
m = np.asarray(Image.open(f"{E}/mask_1200x800.png")) > 0
for c in conds:
    out = []
    for meth in methods:
        v = []
        for sc in scenes:
            try: v.append(10 * np.log10(255 ** 2 / ((rd(f"{E}/{sc}_{c}_{meth}.png") - rd(f"{E}/{sc}_{c}_gt.png")) ** 2)[m].mean()))
            except FileNotFoundError: pass
        out.append(f"{meth} {np.mean(v):6.2f}" if v else f"{meth} n/a")
    print(f"{c:5s} ({len(v)} scenes)  " + "  ".join(out))
