using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.App.IO;
using FisheyeFriend.App.UI;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>Open one photo, inspect it pixel-perfectly, defish it and save the result.</summary>
    internal sealed class InteractiveMode : IDisposable
    {
        private readonly FisheyeFriendApp _app;
        private readonly ImageViewer _viewer = new ImageViewer();
        private readonly ProgressTracker _progress = new ProgressTracker();
        private readonly int[] _px = new int[4];

        private ImageLayer _originalLayer, _resultLayer;
        private PixelImage _original, _result;
        private DefishSettings _resultSettings;
        private DefishGeometry _resultGeometry;
        private string _path = "";
        private string _pathField = "";
        private bool _showResult;
        private bool _defishOnOpen;
        private CancellationTokenSource _cts;
        private string _busy;                       // non-null while working
        private FileDialogs.SaveFormat _saveFormat;
        private string _lastSaved;

        public InteractiveMode(FisheyeFriendApp app)
        {
            _app = app;
            _viewer.SetMaxZoom(app.Config.maxZoom);
            _saveFormat = app.Config.defaultSaveFormat == "png" ? FileDialogs.SaveFormat.Png : FileDialogs.SaveFormat.Jpeg;
            _saveFormat = Prefs.GetString("SaveFormat", _saveFormat.ToString()) == "Png" ? FileDialogs.SaveFormat.Png : FileDialogs.SaveFormat.Jpeg;
            _defishOnOpen = Prefs.GetBool("DefishOnOpen", app.Config.defishOnOpen);
        }

        public bool Busy => _busy != null;
        public PixelImage Original => _original;

        private bool ResultIsCurrent => _result != null && _resultSettings != null && _resultSettings.Same(_app.Settings);

        // ============================================================================ main area

        public void DrawMain(Rect area, UiStyles ui)
        {
            float barH = ui.S(42);
            float statusH = ui.S(26);
            var toolbar = new Rect(area.x, area.y, area.width, barH);
            var view = new Rect(area.x, area.y + barH, area.width, area.height - barH - statusH);
            var status = new Rect(area.x, area.yMax - statusH, area.width, statusH);

            // Snap the view to whole device pixels so the image grid lines up with the screen.
            float ppp = UiStyles.PixelsPerPoint;
            view = new Rect(Mathf.Round(view.x * ppp) / ppp, Mathf.Round(view.y * ppp) / ppp,
                Mathf.Round(view.width * ppp) / ppp, Mathf.Round(view.height * ppp) / ppp);

            _viewer.HandleInput(view);
            _viewer.Draw(view, ui);
            if (_viewer.Layer == null) DrawEmptyHint(view, ui);
            DrawToolbar(toolbar, ui);
            DrawStatus(status, ui);
        }

        private void DrawEmptyHint(Rect view, UiStyles ui)
        {
            var style = new GUIStyle(ui.Label) { alignment = TextAnchor.MiddleCenter, wordWrap = true };
            style.normal.textColor = UiStyles.Muted;
            GUI.Label(view, "Open a JPEG or PNG shot with the fisheye lens\n(Open…, or type a path above and press Enter)", style);
        }

        private void DrawToolbar(Rect r, UiStyles ui)
        {
            GUILayout.BeginArea(r, ui.Bar);
            GUILayout.BeginHorizontal();
            GUI.enabled = !Busy;
            if (GUILayout.Button(new GUIContent("Open…", "Open a photo (Ctrl/Cmd+O)"), ui.Button, GUILayout.Width(ui.S(80))))
                _app.Defer(OpenWithDialog);
            if (Event.current.type == EventType.KeyDown && (Event.current.keyCode == KeyCode.Return || Event.current.keyCode == KeyCode.KeypadEnter) &&
                GUI.GetNameOfFocusedControl() == "pathField")
            {
                Event.current.Use();
                GUI.FocusControl(null);
                string typed = _pathField.Trim().Trim('"');
                _app.Defer(() => Open(typed));
            }
            GUI.SetNextControlName("pathField");
            _pathField = GUILayout.TextField(_pathField, ui.TextField, GUILayout.MinWidth(ui.S(120)), GUILayout.ExpandWidth(true));
            GUI.enabled = !Busy && _original != null;
            if (GUILayout.Button(new GUIContent("<", "Previous photo in this folder ([ or Page Up)"), ui.Button, GUILayout.Width(ui.S(30))))
                _app.Defer(() => Step(-1));
            if (GUILayout.Button(new GUIContent(">", "Next photo in this folder (] or Page Down)"), ui.Button, GUILayout.Width(ui.S(30))))
                _app.Defer(() => Step(+1));
            GUI.enabled = true;

            GUILayout.Space(ui.S(10));
            bool hasResult = _result != null;
            if (GUILayout.Button(new GUIContent("Original", "Show the fisheye photo (Space or Tab toggles)"), !_showResult || !hasResult ? ui.ButtonSelected : ui.Button, GUILayout.Width(ui.S(84))))
                SetShowResult(false);
            GUI.enabled = hasResult;
            if (GUILayout.Button(new GUIContent("Defished", "Show the result (Space or Tab toggles)"), _showResult && hasResult ? ui.ButtonSelected : ui.Button, GUILayout.Width(ui.S(84))))
                SetShowResult(true);
            GUI.enabled = true;

            GUILayout.Space(ui.S(10));
            if (GUILayout.Button(new GUIContent("-", "Zoom out (-)"), ui.Button, GUILayout.Width(ui.S(30)))) _viewer.ZoomOut();
            GUILayout.Label(_viewer.ZoomLabel, new GUIStyle(ui.Label) { alignment = TextAnchor.MiddleCenter }, GUILayout.Width(ui.S(92)));
            if (GUILayout.Button(new GUIContent("+", "Zoom in (+)"), ui.Button, GUILayout.Width(ui.S(30)))) _viewer.ZoomIn();
            if (GUILayout.Button(new GUIContent("Fit", "Largest exact zoom that fits (F)"), ui.Button, GUILayout.Width(ui.S(44)))) _viewer.Fit();
            if (GUILayout.Button(new GUIContent("1:1", "Actual pixels (1)"), ui.Button, GUILayout.Width(ui.S(44)))) _viewer.ActualSize();
            GUILayout.EndHorizontal();
            _app.CaptureTooltip();
            GUILayout.EndArea();
        }

        private void DrawStatus(Rect r, UiStyles ui)
        {
            ui.Fill(r, UiStyles.Hex(0x1D1F22));
            var img = _viewer.Layer != null ? _viewer.Layer.Image : null;
            string text;
            if (img == null)
            {
                text = "No image";
            }
            else
            {
                string which = _showResult && _result != null ? "Defished" : "Original";
                text = which + "  ·  " + img.Describe() + "  ·  " + Path.GetFileName(_path);
                var h = _viewer.Hover;
                if (h.HasValue)
                {
                    img.GetPixel(h.Value.x, h.Value.y, _px);
                    text += "  ·  (" + h.Value.x + ", " + h.Value.y + ")  ";
                    for (int c = 0; c < img.Channels; c++) text += (c > 0 ? " " : "") + _px[c];
                }
                if (_showResult && _result != null && !ResultIsCurrent) text += "  ·  settings changed since Defish";
            }
            var style = new GUIStyle(ui.MutedSmall) { alignment = TextAnchor.MiddleLeft, wordWrap = false };
            GUI.Label(new Rect(r.x + ui.S(10), r.y, r.width - ui.S(20), r.height), text, style);
        }

        // ============================================================================ side panel actions

        public void DrawActions(UiStyles ui)
        {
            if (Busy)
            {
                ui.ProgressBar(_progress.Value, _busy + "  " + Mathf.RoundToInt(_progress.Value * 100) + "%");
                if (_cts != null && GUILayout.Button("Cancel", ui.Button)) _cts.Cancel();
                return;
            }

            GUI.enabled = _original != null && !_app.Projection.Failed;
            string applyLabel = _result != null && ResultIsCurrent ? "Defished" : "Defish";
            if (GUILayout.Button(new GUIContent(applyLabel, "Defish at full resolution with the settings above (Ctrl/Cmd+Enter)"), ui.Primary))
                _app.Defer(() => Apply(false));
            GUI.enabled = true;
            bool auto = ui.Toggle(_defishOnOpen, new GUIContent("Defish when opened", "Defish each photo as soon as it opens (also with < and >)"));
            if (auto != _defishOnOpen)
            {
                _defishOnOpen = auto;
                Prefs.SetBool("DefishOnOpen", auto);
            }

            GUI.enabled = _original != null && !_app.Projection.Failed;
            GUILayout.Space(ui.S(6));
            GUILayout.BeginHorizontal();
            GUILayout.Label("Save as", ui.Label, GUILayout.Width(ui.S(60)));
            if (GUILayout.Button(new GUIContent("JPG", "JPEG at maximum quality: Q100, 4:4:4, optimised Huffman"), _saveFormat == FileDialogs.SaveFormat.Jpeg ? ui.ButtonSelected : ui.Button))
                SetSaveFormat(FileDialogs.SaveFormat.Jpeg);
            if (GUILayout.Button(new GUIContent("PNG", "Lossless PNG; 16-bit sources stay 16-bit"), _saveFormat == FileDialogs.SaveFormat.Png ? ui.ButtonSelected : ui.Button))
                SetSaveFormat(FileDialogs.SaveFormat.Png);
            GUILayout.EndHorizontal();
            if (GUILayout.Button(new GUIContent("Save As…", "Save the defished photo (Ctrl/Cmd+S). Defishes first if needed."), ui.Button))
                _app.Defer(SaveAs);
            GUI.enabled = true;

            if (!string.IsNullOrEmpty(_lastSaved))
                GUILayout.Label("Saved " + _lastSaved, ui.Success);
        }

        private void SetSaveFormat(FileDialogs.SaveFormat f)
        {
            _saveFormat = f;
            Prefs.SetString("SaveFormat", f.ToString());
        }

        private void SetShowResult(bool show)
        {
            bool target = show && _result != null;
            var layer = target ? _resultLayer : _originalLayer;
            if (layer == null) return;
            // Keep the same part of the scene under the view: the projection moves it, so map the view centre through it.
            double cx = 0, cy = 0;
            bool mapped = false;
            if (target != _showResult && _viewer.Layer != null && _resultGeometry != null && _resultSettings != null)
            {
                double x, y;
                _viewer.GetCentre(out x, out y);
                if (target)
                {
                    mapped = Defisher.PhotoToResult(_resultGeometry, _original, _resultSettings, x - 0.5, y - 0.5, out cx, out cy);
                }
                else
                {
                    Defisher.ResultToPhoto(_resultGeometry, _original, _resultSettings, x - 0.5, y - 0.5, out cx, out cy);
                    mapped = true;
                }
            }
            _showResult = target;
            _viewer.Show(layer, true);
            if (mapped) _viewer.SetCentre(cx + 0.5, cy + 0.5);
        }

        // ============================================================================ keyboard

        public void HandleShortcuts(Event e)
        {
            if (e.type != EventType.KeyDown) return;
            bool cmd = e.command || e.control;
            if (cmd && e.keyCode == KeyCode.O) { if (!Busy) _app.Defer(OpenWithDialog); e.Use(); return; }
            if (cmd && e.keyCode == KeyCode.S) { if (!Busy && _original != null) _app.Defer(SaveAs); e.Use(); return; }
            if (cmd && (e.keyCode == KeyCode.Return || e.keyCode == KeyCode.KeypadEnter)) { if (!Busy && _original != null) _app.Defer(() => Apply(false)); e.Use(); return; }
            if (cmd) return;
            if (GUIUtility.keyboardControl != 0) return; // typing in a text field
            switch (e.keyCode)
            {
                case KeyCode.Space:
                case KeyCode.Tab: SetShowResult(!_showResult); e.Use(); break;
                case KeyCode.Equals:
                case KeyCode.Plus:
                case KeyCode.KeypadPlus: _viewer.ZoomIn(); e.Use(); break;
                case KeyCode.Minus:
                case KeyCode.KeypadMinus: _viewer.ZoomOut(); e.Use(); break;
                case KeyCode.Alpha1:
                case KeyCode.Alpha0:
                case KeyCode.Keypad1: _viewer.ActualSize(); e.Use(); break;
                case KeyCode.F: _viewer.Fit(); e.Use(); break;
                case KeyCode.LeftArrow:
                case KeyCode.A: _viewer.PanPixels(Pan(e), 0); e.Use(); break;
                case KeyCode.RightArrow:
                case KeyCode.D: _viewer.PanPixels(-Pan(e), 0); e.Use(); break;
                case KeyCode.UpArrow:
                case KeyCode.W: _viewer.PanPixels(0, Pan(e)); e.Use(); break;
                case KeyCode.DownArrow:
                case KeyCode.S: _viewer.PanPixels(0, -Pan(e)); e.Use(); break;
                case KeyCode.LeftBracket:
                case KeyCode.PageUp: if (!Busy && _original != null) _app.Defer(() => Step(-1)); e.Use(); break;
                case KeyCode.RightBracket:
                case KeyCode.PageDown: if (!Busy && _original != null) _app.Defer(() => Step(+1)); e.Use(); break;
            }
        }

        private static float Pan(Event e)
        {
            return e.shift ? 400f : 80f;
        }

        // ============================================================================ actions

        private async void OpenWithDialog()
        {
            if (Busy) return;
            string start = !string.IsNullOrEmpty(_path) ? Path.GetDirectoryName(_path) : Prefs.GetString("OpenDir", "");
            string p = await FileDialogs.OpenImageAsync("Open fisheye photo", start);
            if (!string.IsNullOrEmpty(p)) Open(p);
        }

        /// <summary>Opens the previous (-1) or next (+1) photo in the open photo's folder, by name.</summary>
        private void Step(int dir)
        {
            if (Busy || string.IsNullOrEmpty(_path)) return;
            string folder = Path.GetDirectoryName(Path.GetFullPath(_path));
            List<string> files;
            try
            {
                files = new List<string>();
                foreach (var f in Directory.GetFiles(folder))
                    if (ImageIO.IsSupportedPath(f)) files.Add(f);
                files.Sort(StringComparer.OrdinalIgnoreCase);
            }
            catch (Exception e)
            {
                _app.ShowError("Could not list " + folder + ": " + e.Message);
                return;
            }
            int i = files.FindIndex(f => string.Equals(Path.GetFullPath(f), Path.GetFullPath(_path), StringComparison.OrdinalIgnoreCase));
            int j = i < 0 ? 0 : i + dir;
            if (j < 0 || j >= files.Count)
            {
                _app.ShowInfo(dir < 0 ? "This is the first photo in the folder." : "This is the last photo in the folder.");
                return;
            }
            Open(files[j]);
        }

        public async void Open(string path)
        {
            if (Busy || string.IsNullOrEmpty(path)) return;
            if (!File.Exists(path))
            {
                _app.ShowError("File not found: " + path);
                return;
            }
            _busy = "Loading";
            _progress.BeginStage("Loading", 0f, 1f);
            _cts = null;
            bool opened = false;
            try
            {
                var img = await ImageLoader.LoadAsync(path, CancellationToken.None);
                var layer = await ImageLayer.CreateAsync(img);
                DisposeLayers();
                _original = img;
                _originalLayer = layer;
                _result = null;
                _resultSettings = null;
                _resultGeometry = null;
                _showResult = false;
                _path = path;
                _pathField = path;
                _lastSaved = null;
                _viewer.Show(_originalLayer, false);
                Prefs.SetString("OpenDir", Path.GetDirectoryName(path));
                int o = img.Metadata != null ? img.Metadata.SourceOrientation : 0;
                _app.ShowInfo("Opened " + Path.GetFileName(path) + " (" + img.Describe() + ")" +
                              (o > Orientation.Normal ? " - " + Orientation.Describe(o) + " as the camera's orientation tag asks" : ""));
                opened = true;
                GC.Collect();
            }
            catch (Exception e)
            {
                _app.ShowError("Could not open " + Path.GetFileName(path) + ": " + e.Message);
            }
            finally
            {
                _busy = null;
            }
            if (opened && _defishOnOpen && !_app.Projection.Failed) await ApplyAsync();
        }

        private async void Apply(bool thenSave)
        {
            await ApplyAsync();
            if (thenSave && ResultIsCurrent) SaveAs();
        }

        private async Task ApplyAsync()
        {
            if (Busy || _original == null) return;
            var settings = _app.Settings.Clone();
            var src = _original;
            _cts = new CancellationTokenSource();
            var ct = _cts.Token;
            _busy = "Defishing";
            _progress.BeginStage("Preparing the projection", 0f, 0f);
            var sw = System.Diagnostics.Stopwatch.StartNew();
            try
            {
                if (!_app.Projection.Ready) _busy = "Preparing the projection";
                else if (_app.Projection.PeekGeometry(Defisher.SensorAspect(src), settings.Adjust) == null) _busy = "Preparing the defish";
                var g = await _app.Projection.GeometryAsync(Defisher.SensorAspect(src), settings.Adjust);
                ct.ThrowIfCancellationRequested();
                _busy = "Defishing";
                sw.Restart();
                _progress.BeginStage("Resampling", 0f, 0.85f);
                var result = await Task.Run(() => Defisher.Process(src, g, settings, _progress, ct), ct);
                if (src != _original) return;   // another photo was opened meanwhile
                _progress.BeginStage("Preparing view", 0.85f, 0.15f);
                var layer = await ImageLayer.CreateAsync(result);
                var previous = _resultLayer;
                bool wasShowingResult = _showResult && previous != null;
                _result = result;
                _resultLayer = layer;
                _resultSettings = settings;
                _resultGeometry = g;
                _lastSaved = null;
                if (wasShowingResult)
                {
                    // Compare results of different settings at the same spot.
                    _viewer.Show(_resultLayer, true);
                }
                else
                {
                    _showResult = false;   // so the toggle maps the view through the projection
                    SetShowResult(true);
                }
                if (previous != null) previous.Dispose();
                _app.ShowInfo("Defished to " + result.Width + " x " + result.Height + " in " +
                              sw.Elapsed.TotalSeconds.ToString("0.0", CultureInfo.InvariantCulture) + " s (" + Describe(settings) + ")");
                GC.Collect();
            }
            catch (OperationCanceledException)
            {
                _app.ShowInfo("Cancelled.");
            }
            catch (Exception e)
            {
                _app.ShowError("Defish failed: " + e.Message);
            }
            finally
            {
                _busy = null;
                _cts = null;
            }
        }

        internal static string Describe(DefishSettings s)
        {
            return SettingsPanel.FilterName(s) + ", " + (s.OutputSize == DefishOutputSize.FullDetail ? "full detail" : "photo size") +
                   (s.Edges == DefishEdges.Crop ? ", cropped to the photo" : "") + DescribeAdjustment(s.Adjust);
        }

        internal static string DescribeAdjustment(DefishAdjustment a)
        {
            if (a == null || a.IsDefault) return "";
            var ci = CultureInfo.InvariantCulture;
            string t = ", lens " + (a.Strength * 100).ToString("0", ci) + "%";
            if (Math.Abs(a.AmountX - a.AmountY) < 1e-9) t += ", amount " + (a.AmountX * 100).ToString("0", ci) + "%";
            else t += ", horizontal " + (a.AmountX * 100).ToString("0", ci) + "%, vertical " + (a.AmountY * 100).ToString("0", ci) + "%";
            if (!a.ShapesFlat) t += ", curved";
            return t;
        }

        private async void SaveAs()
        {
            if (Busy || _original == null) return;
            if (!ResultIsCurrent)
            {
                await ApplyAsync();
                if (!ResultIsCurrent) return;
            }
            string dir = Prefs.GetString("SaveDir", Path.GetDirectoryName(_path));
            string ext = _saveFormat == FileDialogs.SaveFormat.Png ? ".png" : ".jpg";
            string suffix = string.IsNullOrEmpty(_app.Config.batchFileSuffix) ? "_defished" : _app.Config.batchFileSuffix;
            string name = Path.GetFileNameWithoutExtension(_path) + suffix + ext;
            string target = await FileDialogs.SaveImageAsync("Save defished photo", dir, name, _saveFormat);
            if (string.IsNullOrEmpty(target)) return;

            // The extension the user typed wins; with none, use the selected format.
            var fmt = ImageIO.FormatFromExtension(target);
            if (fmt == ImageFileFormat.Unknown)
            {
                fmt = _saveFormat == FileDialogs.SaveFormat.Png ? ImageFileFormat.Png : ImageFileFormat.Jpeg;
                target += ImageIO.DefaultExtension(fmt);
            }
            if (string.Equals(Path.GetFullPath(target), Path.GetFullPath(_path), StringComparison.OrdinalIgnoreCase))
            {
                _app.ShowError("Choose a different file name - this would overwrite the original.");
                return;
            }

            var img = _result;
            _busy = "Saving";
            _progress.BeginStage("Encoding", 0f, 1f);
            try
            {
                await Task.Run(() =>
                {
                    byte[] bytes = ImageIO.Encode(img, fmt);
                    _progress.Report(0.9f);
                    FileUtil.WriteAllBytesSafely(target, bytes);
                });
                Prefs.SetString("SaveDir", Path.GetDirectoryName(target));
                _lastSaved = Path.GetFileName(target);
                _app.ShowInfo("Saved " + target + (fmt == ImageFileFormat.Jpeg ? " (JPEG Q100, 4:4:4)" : " (PNG)"));
            }
            catch (Exception e)
            {
                _app.ShowError("Save failed: " + e.Message);
            }
            finally
            {
                _busy = null;
            }
        }

        private void DisposeLayers()
        {
            if (_originalLayer != null) _originalLayer.Dispose();
            if (_resultLayer != null) _resultLayer.Dispose();
            _originalLayer = _resultLayer = null;
            _viewer.Show(null, false);
        }

        public void Dispose()
        {
            if (_cts != null) _cts.Cancel();
            DisposeLayers();
        }
    }

    internal static class FileUtil
    {
        /// <summary>Writes to a temporary file first, then replaces the target, so a failed write never leaves a half file.</summary>
        public static void WriteAllBytesSafely(string path, byte[] bytes)
        {
            string dir = Path.GetDirectoryName(Path.GetFullPath(path));
            if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
            string tmp = path + ".partial";
            File.WriteAllBytes(tmp, bytes);
            if (File.Exists(path)) File.Delete(path);
            File.Move(tmp, path);
        }
    }
}
