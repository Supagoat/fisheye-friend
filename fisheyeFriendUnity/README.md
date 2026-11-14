# Fisheye Friend

A Unity (6.3) desktop app that defishes photos taken with a fisheye lens: open a photo, inspect it
pixel-perfectly, defish it, compare it with the original, and save it as a maximum-quality JPEG or PNG. A batch
mode does the same for every photo in a folder.

The defish is a model of the Hemi projection: a lens and a projection fitted to Fisheye Friend's 596 hand-measured
point pairs (within 0.9 px rms of them). Like the pairs it straightens vertical lines and takes away the fisheye's
sideways stretch while keeping a gentle horizontal curve, so people look natural. Unlike the pairs it can be
adjusted: for a stronger or milder fisheye, for more or less correction horizontally and vertically, and for how
the correction changes from the centre to the corners (see [The defish model](#the-defish-model)).

It is built like the de-squeeze app, with the same window layout, viewer, side panel, dialogs, shortcuts, file
handling and batch mode, so the two work alike. All image processing lives in a plain C# library
(`Assets/FisheyeFriend/Core`, no Unity references), so it can be tested outside Unity.

The previous version of the project (Unity 2021, uGUI, StandaloneFileBrowser, System.Drawing) is kept
untouched in `Legacy (Unity 2021)/` (its `Assets`, `Packages` and `ProjectSettings`).

---

## Running it

1. Open this folder in Unity 6.3 (6000.3.11f1).
2. Press **Play** in any scene (`Assets/Scenes/SampleScene` is the build scene). The app creates itself at
   start-up; nothing needs adding to a scene.
3. To build: **File > Build Profiles**, pick macOS (or Windows/Linux), **Build**.

On first load an editor script sets desktop-friendly player settings (windowed, resizable, 1600 x 1000,
keeps running in the background so batches don't pause, Retina on). Re-apply any time with
**Fisheye Friend > Apply Recommended Player Settings**. The same menu opens the config file and clears the
remembered settings.

## Interactive mode

| Action | How |
|---|---|
| Open | **Open…**, or type/paste a path in the toolbar and press Enter, or Ctrl/Cmd+O |
| Previous / next photo in the folder | **<** **>** in the toolbar, `[` `]` or Page Up / Page Down |
| Defish | **Defish** (Ctrl/Cmd+Enter): full resolution, on all cores, cancellable. With **Defish when opened** on (the default) every photo is defished as it opens |
| Toggle original / defished | **Original / Defished** buttons, Space or Tab. The zoom is kept and the view follows the projection, so the same part of the scene stays under the view |
| Zoom | Mouse wheel (at the cursor), **+ / −**, **Fit**, **1:1**, keys `+` `-` `1` `F` |
| Pan | Drag, arrow keys or W A S D (Shift = faster) |
| Save | Choose **JPG** or **PNG**, then **Save As…** (Ctrl/Cmd+S). Defishes first if the settings changed |

Viewing is pixel-perfect: magnification is a whole number (point-sampled), and zooming out uses exact 1/N
reductions averaged in linear light on the CPU. The status bar shows the pixel under the cursor and its values.

## Batch mode

Pick an input and an output folder. Every `.jpg/.jpeg/.png` is defished with the current settings and written
**in its own format**: JPEG in gives maximum-quality JPEG out (Q100, 4:4:4), PNG in gives PNG out (16-bit stays
16-bit). Options: include subfolders (mirrors the structure), overwrite existing, file-name suffix (default
`_defished`). The originals can never be overwritten. Each file's row shows its result size and time.

## Settings (side panel)

* **Preview.** With a photo open, the top of the panel shows the whole defished frame, redrawn as you change the
  settings below (a reduced copy of the photo through a coarse version of the projection; tens of milliseconds).
  **Defish** renders it at full size.
* **Defish.**
  * *Lens strength* (50-115%, with the lens' diagonal field of view): how strongly the fisheye bends. 100% is the
    lens the defish was measured on (177° across the diagonal of a 3:2 frame). If vertical lines near the sides
    still bow outwards after defishing, raise it; if they bow inwards, lower it. The width of the result is kept.
  * *Horizontal* and *Vertical* (0-150%): how far the defish moves pixels sideways and up and down. 0% keeps the
    photo's positions on that axis, 100% is the measured defish, more goes further. *Link horizontal and
    vertical* (on by default) moves both together, and the two curves below.
  * *Across the frame*: a curve for each axis with four points (the centre, a third and two thirds of the way
    out, the corners), 0-200% of the amounts above. Drag a point (Alt/Option: exact values, no snapping to 100%),
    double-click it for 100%. The dashed lines mark where the frame's top and side edges are. Unlinked, the
    buttons choose which curve to edit; the other is drawn dimmed. *Flatten curves* sets every point to 100%.
  * *Reset defish* restores the defaults. Settings that would fold the picture over itself (only possible with
    steep curves at high amounts) are flagged and can't be rendered.
* **Output size.** *Full detail* (default) or *Same as photo* (what the old Fisheye Friend made). The Hemi
  projection reduces the centre of the frame 1.19x and enlarges the edges up to 1.8x. At the photo's size the
  centre therefore loses its finest detail and, seen at the same scale, looks softer than the original. Full
  detail makes the result 1.19x larger each way (1.42x the pixels), so the centre keeps the photo's sharpness
  (see [Sharpness](#sharpness)). The factor follows the adjustment (from 1.04x at 50% lens strength to 1.30x at
  115%); the panel shows the size, with ≈ while it is still estimating from the preview.
* **Edges.** *Full frame* (default, as before) keeps the whole frame; along the top and bottom the projection
  reaches up to about 5 pixels past the photo (on a 42 MP frame), and those pixels repeat its edge. *Crop to
  photo* crops 0.3% each way so every pixel comes from the photo.
* **Filter.** Lanczos 4 (default), Lanczos 3, Bicubic.
* **Advanced.** Light space (Sigmoidal default), aperture correction (on), anti-aliasing of the centre (off),
  reset to defaults.

Portrait photos are defished in the sensor's landscape frame (turned the way the camera's orientation tag
says, defished, turned back), so the projection lines up with the lens. The model works in the lens' own terms,
so photos of any shape (4:3, 16:9, crops) are defished correctly; a geometry is made for each frame shape (about
half a second) and the last few are kept. (With the old point-pair file, photos that are not 3:2 are stretched
to its frame, and the panel says so.) EXIF data and ICC profiles are carried to the
output (and a PNG's text chunks to a PNG); the orientation tag is reset as the pixels are upright, and the
EXIF pixel dimensions are updated.

`Assets/StreamingAssets/fisheye-friend-config.json` holds the defaults (including `defishStrength`,
`defishHorizontal` and `defishVertical`, in percent) and the projection file: `defish-hemi-model.json`, the model.
Setting `projectionFile` back to `defish-hemi-596.json` (Fisheye Friend's 596 point pairs, kept unchanged) uses the
measured pairs directly, as fixed splines that can't be adjusted. Without either file the built-in copy of the
model is used. The config can be edited in a built app too; restart after editing.

---

## The defish model

Fisheye Friend's defish was defined by hand: 596 pairs of points, each a spot in a fisheye photo and where it
belongs in the defished picture, joined by thin-plate splines. Studying the pairs showed that they describe a
physical projection: they are symmetric left-right and top-bottom, the defished x depends only on the scene's
longitude (vertical lines come out straight, to 0.6 px), and y follows a Panini-like cylindrical projection
whose vertical scale grows towards the sides (the horizontal curve). The model is that projection:

* **The lens** (fisheye photo): a ray at angle θ from the axis lands at distance r from the centre, with
  θ(r) = t1·r + t3·r³ + t5·r⁵ + t7·r⁷ (r in photo widths). Fitted: 177° across the diagonal of a 3:2 frame.
* **The projection** (defished picture), for a ray at longitude λ and latitude φ:
  x = k·(P+1)·sin λ / (P + cos λ) · (1 + h1·sin²λ + h2·sin⁴λ), and
  y = k·(Q+1) / (Q + cos λ) · tan φ · (1 + Σ v_ij·sin^(2i+2) φ·sin^(2j) λ).
  P is the horizontal Panini parameter, Q the vertical one (Q > P makes the vertical scale grow towards the sides).

21 numbers in all (in `defish-hemi-model.json`), chosen by cross-validation from a family of fits: fewer terms
miss the pairs' shape, more start following the noise in the hand placement. Against the 596 pairs' splines over
the whole frame the model differs by 0.89 px rms, at most 1.85 px (6000 x 4000 frame); on a 42 MP photo the two
renders differ by 1.2 px rms, in a pattern that alternates from place to place: the splines' small wobbles from
points placed a pixel off, which the model, being smooth and exactly symmetric, doesn't have.

**The adjustment** works on that model:

* *Lens strength* s divides the lens' focal length: r_s(θ) = r(θ)/s, so a stronger lens fits a wider view on the
  frame. k is scaled so the photo's left and right edges stay at the result's edges.
* *Horizontal* and *Vertical* amounts a blend each axis of the mapping with no change at all: the photo point
  read for a result point p is p + a·(model(p) − p), separately for x and y. 0% leaves the photo as it is.
* *Across the frame*: each amount is multiplied by its curve at the point's distance from the centre (0 at the
  centre, 1 at the corners), a smooth monotone cubic through the four points (no overshoot between them).

Over the range offered (strength 50-115%, amounts 0-150%, flat curves) the result never folds; the smallest area
scale is 0.125. The geometry is checked whenever it is made, and settings that fold are refused.

---

## Output quality

### How it was measured

There is no "true" defished version of a real photo, so the comparison was made where the truth is known.
Seven scenes were made from the sample photos (the three fisheye samples and four other street photos), at
4800 x 3200. Each was then:

* **photographed through the fisheye**: a 1200 x 800 capture where each sensor pixel averages 4 x 4 samples
  of the scene, taken through the exact inverse of the projection, in linear light. This is the photo Fisheye
  Friend would be given. A second, *soft* capture adds a 0.4 px Gaussian blur (lens and anti-alias filter).
* **photographed rectilinearly**: the 4 x 4 box average of the scene, i.e. what a camera with a rectilinear
  lens and the output's pixels would have captured. This is the ground truth.

Every method defished the capture, and the result was compared with the truth (PSNR over the pixels whose
source is at least 6 px inside the photo). Fisheye Friend's own renderer was reproduced exactly: the replica
matches its actual output for `nov11Street-1.jpg` to 79 dB, so the samples were made by the CPU "Save sharp"
path and the replica is that path.

### Results (mean of the 7 scenes, dB)

| Method | clean | soft |
|---|---|---|
| Fisheye Friend as it was ("Save sharp": radial Lanczos 3, 8-bit gamma) | 24.96 | 25.74 |
| the same with its half-pixel offset fixed | 32.70 | 34.80 |
| radial Lanczos 3, sigmoidal light, elliptical anti-aliasing | 33.05 | 35.02 |
| separable Lanczos 3, gamma light, anti-aliased | 37.17 | 36.01 |
| separable Lanczos 3, linear light, anti-aliased | 35.83 | 35.47 |
| separable Lanczos 3, sigmoidal light, anti-aliased | 37.32 | 36.15 |
| the same, anti-aliasing at 75 / 50 / 25% strength | 37.37 / 37.41 / 37.43 | 36.21 / 36.26 / 36.30 |
| the same, no anti-aliasing | 37.42 | 36.33 |
| separable Lanczos 4, sigmoidal, no anti-aliasing | 37.46 | 36.38 |
| **+ aperture correction (the new default)** | **37.54** | **36.66** |
| Lanczos 3 with the same settings | 37.52 | 36.61 |
| Bicubic (Catmull-Rom) with the same settings | 37.02 | 35.83 |
| aperture correction also smoothing where the projection shrinks | 37.55 | 36.52 |
| aperture correction at 1.5x / 2x strength | 37.51 / 37.44 | 36.77 / 36.85 |

Against the output Fisheye Friend actually produced, the new default is **+12.6 dB** (clean) / **+10.9 dB**
(soft). Most of that is a half-pixel shift in the old renderer (it measured each output pixel's source from the
pixel's corner instead of its centre), which moves the whole image but is invisible on its own. Against the old
renderer with that fixed, the new default is **+4.8 dB / +1.9 dB**, and that part is sharper detail: the
radial Lanczos kernel the old renderer (and its GPU shader) used is not interpolating, so it blurs every
pixel, while the separable kernel reproduces the photo exactly where the projection does not move it.

By region (clean, new default vs. the fixed old renderer): centre 37.7 vs 31.0, middle 37.9 vs 32.1, outer
ring 36.6 vs 32.6, edges 34.7 vs 31.6 dB. The edges stay hardest because there the projection enlarges
the photo up to 1.8x: detail the photo never recorded can't be made up by a filter.

### What each choice does

* **Separable Lanczos** instead of radial: the biggest gain (4.3 dB on crisp captures, 1.1 dB on soft ones).
* **Sigmoidal light**: filtering in a contrast-curve space reduces halos where the projection enlarges; best of
  the three light spaces, as in the de-squeeze app.
* **No anti-aliasing**: widening the filter where the projection shrinks the centre (1.19x) is the textbook
  fix for aliasing, but at so mild a shrink the photo's own pixel blur already prevents it, and widening cost
  0.1-0.2 dB at every strength tried. It stays available under Advanced.
* **Aperture correction**: a camera pixel averages the light over its area. Where the projection enlarges the
  photo, each output pixel covers less of the scene than a photo pixel did, so the photo's pixel blur is partly
  undone with a 3-tap filter whose strength follows the local scale on each axis (the least-squares match of
  the two box apertures; no tuning). Sharpening only, as smoothing where the centre shrinks cost more on soft
  photos than it gained on crisp ones. Higher strengths help soft photos but hurt crisp ones, so the
  theoretical strength (1x) is used.
* **Lanczos 4**: 0.02-0.05 dB over Lanczos 3 for about a third more time; kept as the default because quality
  comes first, with Lanczos 3 one click away.

The same seven scenes and conditions decided every default; the numbers above can be reproduced with
`Tools/Bench` (commands `exp-prep`, `exp-run`, and `score.py`).

### Sharpness

Measured on `nov11Street-1` around the centre of the frame: how much of the original's detail each output keeps
at each scale of detail, relative to the original resampled with no loss at all (1.00 = all of it).

| Detail (cycles per output pixel) | 0.1 | 0.2 | 0.3 | 0.35-0.4 | 0.4-0.45 | 0.45-0.5 |
|---|---|---|---|---|---|---|
| Old Fisheye Friend | 1.04 | 1.18 | 1.50 | 1.62 | 1.73 | 1.83 |
| New, same size as photo | 1.00 | 0.99 | 1.01 | 1.00 | 1.01 | 1.04 |
| New, Full detail | 1.00 | 1.00 | 1.00 | 0.99 | 0.95 | 0.87 |

* **The resampling itself loses nothing measurable.** At the photo's size, the new output is the original
  reduced 1.19x with an ideal filter.
* **The old renderer sharpens.** Its radial Lanczos filter boosts fine detail up to 1.8x, which is why its output
  looks crisper and more contrasty than the original photo; it also boosts grain and adds halos.
* **At the photo's size, softness in the centre is unavoidable.** Reducing the centre 1.19x drops the finest 15%
  of its detail, which an image that size cannot hold. Full detail avoids the reduction, and keeps the original's
  detail to within 1% up to 0.4 cycles per pixel.
* The last 5-13% at the very finest detail is what interpolating between pixels costs, at least half a pixel
  away from the photo's grid. Restoring it would take sharpening that changes from pixel to pixel, which brings
  moiré and grain, so it is left as is.
* Reproduce with `Tools/Bench/mtf.py`.

## Speed

The old "Save sharp" path evaluated the 596-point thin-plate splines (hundreds of logarithms) for every pixel,
which is why it took minutes. Now the projection is evaluated once per adjustment and frame shape on a fine grid
(1025 nodes across) and interpolated bicubically, with derivatives, between nodes: the largest difference from
the exact model is 0.0003 px on a 42 MP frame. The model's grid takes about 0.5 s on 2 cores (the splines'
took 4.5 s); it is made in the background once an adjustment has stayed put for a moment, and the last few
are kept. The preview uses a grid of 129 nodes (a few milliseconds). Sampling then reads the photo through
precomputed separable weights.

On the 2-core machine used for development, a 7952 x 5304 (42 MP) photo defishes in about 12 s with
Lanczos 4 (9 s with Lanczos 3), against 24 s for the old sampling loop alone, even with the new fast map.
The work is spread over all cores, so a desktop with more cores is proportionally faster. The GPU preview of
the old version is gone: its shader used the same blurring radial kernel, and the CPU path is now fast enough
to show the final result directly.

---

## Project layout

```
Assets/FisheyeFriend/Core      the library: the defish model and its adjustment (DefishModel), the point pairs
                               (HemiMapping, ThinPlateSpline), DefishGeometry, Defisher, DefishPreview, the warp
                               resampler, codecs (JPEG/PNG with metadata), batch planning
Assets/FisheyeFriend/App       the app: FisheyeFriendApp (shell), InteractiveMode, BatchMode, Projection,
                               SettingsPanel (with the curve editor), LivePreview, viewer, styles, native file
                               dialogs, prefs, config
Assets/FisheyeFriend/Editor    player-settings setup and the Fisheye Friend menu
Assets/FisheyeFriend/Tests     EditMode tests (Window > General > Test Runner)
Assets/StreamingAssets         fisheye-friend-config.json, defish-hemi-model.json (the model), defish-hemi-596.json
                               (the original point pairs)
Tools/Bench                    command-line bench: self-tests, the quality experiment, timing
Legacy (Unity 2021)            the previous project, unchanged
```

The bench's `selftest` checks, among other things: the splines reproduce the 596 points; the grid matches the
exact splines to 0.002 px; the derivatives match finite differences; the inverse mapping round-trips; the warp
is exact for the identity in every light space and bit depth; a flat field stays flat; RGB and RGBA agree;
cropped outputs never read outside the photo; portraits equal the turned landscape result; and a dot in the
photo lands in the result exactly where the view-following mapping says. Given the model file
(`selftest defish-hemi-596.json defish-hemi-model.json`) it also checks the model: the file equals the built-in
copy and round-trips; it matches the pairs' splines (rms < 1.2 px, max < 3 px); its grid error; symmetry; 0%
leaves the photo unchanged; every strength keeps the width; a stronger lens corrects more; no folds over the
offered range while a zig-zag curve is caught and refused; the curves don't overshoot; out-of-range values are
clamped; the preview reproduces the photo with no correction and turns portraits correctly; adjusted portraits
equal the turned landscape result. `model-defish photo.png out.png strength=0.9 h=1 v=0.8 shapex=1,1,0.8,0.6`
renders a photo with any adjustment (`preview=640` gives the panel's quick picture instead).
