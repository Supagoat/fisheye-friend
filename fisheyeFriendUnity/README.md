# Fisheye Friend

A Unity (6.3) desktop app that defishes photos taken with a fisheye lens: open a photo, inspect it
pixel-perfectly, defish it with the measured Hemi projection, compare it with the original, and save it as a
maximum-quality JPEG or PNG. A batch mode does the same for every photo in a folder.

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

* **Output size.** *Same as photo* (default, as before) or *Full detail*: the Hemi projection squeezes the
  centre of the frame by 1.19x (and stretches the edges up to 1.8x), so at the photo's size some of the
  centre's finest detail is averaged away. Full detail makes the result 1.19x larger each way so the centre
  keeps every pixel.
* **Edges.** *Full frame* (default, as before) keeps the whole frame; along the top and bottom the projection
  reaches up to about 5 pixels past the photo (on a 42 MP frame), and those pixels repeat its edge. *Crop to
  photo* crops 0.3% each way so every pixel comes from the photo.
* **Filter.** Lanczos 4 (default), Lanczos 3, Bicubic.
* **Advanced.** Light space (Sigmoidal default), aperture correction (on), anti-aliasing of the centre (off),
  reset to defaults.

Portrait photos are defished in the sensor's landscape frame (turned the way the camera's orientation tag
says, defished, turned back), so the projection lines up with the lens. Photos that are not 3:2 are stretched
to the projection's frame, as before, and the panel says so. EXIF data and ICC profiles are carried to the
output (and a PNG's text chunks to a PNG); the orientation tag is reset as the pixels are upright, and the
EXIF pixel dimensions are updated.

`Assets/StreamingAssets/fisheye-friend-config.json` holds the defaults and the projection file
(`defish-hemi-596.json`, Fisheye Friend's 596 point pairs, same format as before). It can be edited in a
built app too; restart after editing.

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

## Speed

The old "Save sharp" path evaluated the 596-point thin-plate splines (hundreds of logarithms) for every pixel,
which is why it took minutes. Now the splines are evaluated exactly once per session on a fine grid (1025
nodes across: about 4.5 s on 2 cores, done in the background at start-up) and interpolated bicubically,
with derivatives, between nodes: the largest difference from the exact splines is 0.002 px on a 42 MP frame.
Sampling then reads the photo through precomputed separable weights.

On the 2-core machine used for development, a 7952 x 5304 (42 MP) photo defishes in about 12 s with
Lanczos 4 (9 s with Lanczos 3), against 24 s for the old sampling loop alone, even with the new fast map.
The work is spread over all cores, so a desktop with more cores is proportionally faster. The GPU preview of
the old version is gone: its shader used the same blurring radial kernel, and the CPU path is now fast enough
to show the final result directly.

---

## Project layout

```
Assets/FisheyeFriend/Core      the library: projection (HemiMapping, ThinPlateSpline, DefishGeometry),
                               Defisher, the warp resampler, codecs (JPEG/PNG with metadata), batch planning
Assets/FisheyeFriend/App       the app: FisheyeFriendApp (shell), InteractiveMode, BatchMode, SettingsPanel,
                               viewer, styles, native file dialogs, prefs, config
Assets/FisheyeFriend/Editor    player-settings setup and the Fisheye Friend menu
Assets/FisheyeFriend/Tests     EditMode tests (Window > General > Test Runner)
Assets/StreamingAssets         fisheye-friend-config.json, defish-hemi-596.json
Tools/Bench                    command-line bench: self-tests, the quality experiment, timing
Legacy (Unity 2021)            the previous project, unchanged
```

The bench's `selftest` checks, among other things: the splines reproduce the 596 points; the grid matches the
exact splines to 0.002 px; the derivatives match finite differences; the inverse mapping round-trips; the warp
is exact for the identity in every light space and bit depth; a flat field stays flat; RGB and RGBA agree;
cropped outputs never read outside the photo; portraits equal the turned landscape result; and a dot in the
photo lands in the result exactly where the view-following mapping says.
