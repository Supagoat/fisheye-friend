using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.App.IO;
using FisheyeFriend.App.UI;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Batch;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Defish;
using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>
    /// Defishes every JPEG/PNG in a folder with the current settings. Each output keeps its input's format:
    /// JPEG in -> maximum-quality JPEG out, PNG in -> PNG out (bit depth preserved); EXIF and ICC profiles are kept.
    /// </summary>
    internal sealed class BatchMode : IDisposable
    {
        private readonly FisheyeFriendApp _app;
        private readonly BatchOptions _options = new BatchOptions();
        private readonly ProgressTracker _itemProgress = new ProgressTracker();
        private List<BatchItem> _items = new List<BatchItem>();
        private string _planError;
        private Vector2 _scroll;
        private CancellationTokenSource _cts;
        private bool _running;
        private int _current = -1;
        private int _toProcess, _processed;
        private string _summary;
        private readonly Stopwatch _clock = new Stopwatch();

        public BatchMode(FisheyeFriendApp app)
        {
            _app = app;
            Prefs.LoadBatch(_options, app.Config.batchFileSuffix);
        }

        public bool Busy => _running;

        private bool _planned;
        private int _listFirst, _listLast;

        public void DrawMain(Rect area, UiStyles ui)
        {
            if (!_planned)
            {
                _planned = true;
                RequestPlan();
            }
            ui.Fill(area, UiStyles.Background);
            GUILayout.BeginArea(new Rect(area.x + ui.S(16), area.y + ui.S(12), area.width - ui.S(32), area.height - ui.S(24)));
            GUI.enabled = !_running;

            GUILayout.Label("BATCH DEFISH", ui.Title);
            GUILayout.Label("Every .jpg/.jpeg/.png in the input folder is defished with the settings on the right " +
                            "and written to the output folder in its original format (JPEG at maximum quality).", ui.MutedSmall);
            GUILayout.Label("Projection: " + _app.Projection.Describe() + "  ·  " + InteractiveMode.Describe(_app.Settings), _app.Projection.Failed ? ui.Error : ui.Label);

            bool changed = false;
            changed |= FolderRow(ui, "INPUT FOLDER", ref _options.InputFolder, "Choose the folder of images to process");
            changed |= FolderRow(ui, "OUTPUT FOLDER", ref _options.OutputFolder, "Choose where the defished photos go");

            GUILayout.Space(ui.S(6));
            GUILayout.BeginHorizontal();
            bool sub = ui.Toggle(_options.IncludeSubfolders, new GUIContent("Include subfolders", "Recreates the folder structure under the output folder"));
            GUILayout.Space(ui.S(16));
            bool ow = ui.Toggle(_options.OverwriteExisting, new GUIContent("Overwrite existing", "Replace output files that already exist (otherwise they are skipped)"));
            GUILayout.Space(ui.S(16));
            GUILayout.Label(new GUIContent("File suffix", "Added to each output file name, e.g. IMG_0001_defished.jpg"), ui.Label);
            string suffix = GUILayout.TextField(_options.FileSuffix ?? "", ui.TextField, GUILayout.Width(ui.S(150)));
            GUILayout.FlexibleSpace();
            GUILayout.EndHorizontal();
            if (sub != _options.IncludeSubfolders || ow != _options.OverwriteExisting || suffix != _options.FileSuffix)
            {
                _options.IncludeSubfolders = sub;
                _options.OverwriteExisting = ow;
                _options.FileSuffix = suffix;
                changed = true;
            }
            if (changed)
            {
                Prefs.SaveBatch(_options);
                RequestPlan();
            }

            GUILayout.Space(ui.S(8));
            GUILayout.BeginHorizontal();
            if (GUILayout.Button(new GUIContent("Rescan", "List the files again"), ui.Button, GUILayout.Width(ui.S(90)))) RequestPlan();
            int pending = CountPending();
            GUI.enabled = !_running && pending > 0 && !_app.Projection.Failed;
            if (GUILayout.Button(new GUIContent(pending > 0 ? "Start batch (" + pending + ")" : "Start batch", "Defish every listed file"), ui.Primary, GUILayout.Width(ui.S(220))))
                _app.Defer(Run);
            GUI.enabled = true;
            if (_running && GUILayout.Button("Cancel", ui.Button, GUILayout.Width(ui.S(90))) && _cts != null) _cts.Cancel();
            GUI.enabled = !string.IsNullOrEmpty(_options.OutputFolder) && Directory.Exists(_options.OutputFolder);
            if (GUILayout.Button(new GUIContent("Show output folder", "Open the output folder in Finder/Explorer"), ui.Button, GUILayout.Width(ui.S(160))))
                Application.OpenURL(new Uri(Path.GetFullPath(_options.OutputFolder)).AbsoluteUri);
            GUI.enabled = true; // progress and the file list stay fully visible while running
            GUILayout.FlexibleSpace();
            GUILayout.EndHorizontal();

            if (_running)
            {
                float overall = _toProcess > 0 ? (_processed + _itemProgress.Value) / _toProcess : 0f;
                string cur = _current >= 0 && _current < _items.Count ? _items[_current].RelativePath : "";
                ui.ProgressBar(overall, (_processed + 1) + " / " + _toProcess + "   " + cur + "   " + _itemProgress.Stage);
            }
            if (!string.IsNullOrEmpty(_planError)) GUILayout.Label(_planError, ui.Error);
            if (!string.IsNullOrEmpty(_summary)) GUILayout.Label(_summary, ui.Success);

            GUILayout.Space(ui.S(6));
            DrawList(ui);
            GUI.enabled = true;
            _app.CaptureTooltip();
            GUILayout.EndArea();
        }

        private bool FolderRow(UiStyles ui, string label, ref string folder, string dialogTitle)
        {
            GUILayout.Label(label, ui.Section);
            GUILayout.BeginHorizontal();
            string f = GUILayout.TextField(folder ?? "", ui.TextField, GUILayout.ExpandWidth(true));
            bool changed = f != folder;
            folder = f;
            if (GUILayout.Button("Choose…", ui.Button, GUILayout.Width(ui.S(96))))
            {
                bool input = label == "INPUT FOLDER";
                string start = folder;
                _app.Defer(() => PickFolder(input, dialogTitle, start));
            }
            GUILayout.EndHorizontal();
            return changed;
        }

        private async void PickFolder(bool input, string title, string start)
        {
            if (string.IsNullOrEmpty(start)) start = input ? _options.OutputFolder : _options.InputFolder;
            string p = await FileDialogs.PickFolderAsync(title, start);
            if (string.IsNullOrEmpty(p)) return;
            p = p.TrimEnd('/', '\\');
            if (p.Length == 0) p = "/";
            if (input) _options.InputFolder = p;
            else _options.OutputFolder = p;
            Prefs.SaveBatch(_options);
            Plan();
        }

        private void DrawList(UiStyles ui)
        {
            if (_items.Count == 0) return;
            float h = ui.S(22);
            GUILayout.BeginHorizontal();
            GUILayout.Label("Status", ui.Section, GUILayout.Width(ui.S(90)));
            GUILayout.Label("File", ui.Section, GUILayout.Width(ui.S(320)));
            GUILayout.Label("Result", ui.Section, GUILayout.ExpandWidth(true));
            GUILayout.EndHorizontal();
            // Only lay out visible rows so very large folders stay fast. The range is decided in the
            // Layout pass and reused for the other passes of the same event so they match.
            if (Event.current.type == EventType.Layout)
            {
                _listFirst = Mathf.Max(0, (int)(_scroll.y / h) - 2);
                int visible = Mathf.CeilToInt(Screen.height / Mathf.Max(1f, h * UiStyles.PixelsPerPoint)) + 4;
                _listLast = Mathf.Min(_items.Count, _listFirst + visible);
            }
            int first = Mathf.Min(_listFirst, _items.Count);
            int last = Mathf.Min(_listLast, _items.Count);
            _scroll = GUILayout.BeginScrollView(_scroll, GUILayout.ExpandHeight(true));
            GUILayout.Space(first * h);
            for (int i = first; i < last; i++)
            {
                var it = _items[i];
                var style = (i & 1) == 0 ? ui.ListRow : ui.ListRowAlt;
                GUILayout.BeginHorizontal(style, GUILayout.Height(h));
                var st = new GUIStyle(ui.ListRow) { normal = { background = null, textColor = StatusColor(it.Status) } };
                GUILayout.Label(StatusText(it, i), st, GUILayout.Width(ui.S(90)), GUILayout.Height(h));
                GUILayout.Label(it.RelativePath, ui.ListRow, GUILayout.Width(ui.S(320)), GUILayout.Height(h));
                string res = it.Message;
                if (it.Status == BatchItemStatus.Done) res += "   " + it.Seconds.ToString("0.0", CultureInfo.InvariantCulture) + " s   » " + Path.GetFileName(it.OutputPath);
                GUILayout.Label(res, ui.ListRow, GUILayout.ExpandWidth(true), GUILayout.Height(h));
                GUILayout.EndHorizontal();
            }
            GUILayout.Space((_items.Count - last) * h);
            GUILayout.EndScrollView();
        }

        private string StatusText(BatchItem it, int index)
        {
            if (it.Status == BatchItemStatus.Processing) return "Working " + Mathf.RoundToInt(_itemProgress.Value * 100) + "%";
            return it.Status.ToString();
        }

        private static Color StatusColor(BatchItemStatus s)
        {
            switch (s)
            {
                case BatchItemStatus.Done: return UiStyles.Good;
                case BatchItemStatus.Failed: return UiStyles.Bad;
                case BatchItemStatus.Processing: return UiStyles.Accent;
                case BatchItemStatus.Skipped:
                case BatchItemStatus.Cancelled: return UiStyles.Muted;
                default: return UiStyles.Text;
            }
        }

        private int CountPending()
        {
            int n = 0;
            foreach (var it in _items)
                if (it.Status == BatchItemStatus.Pending) n++;
            return n;
        }

        private bool _planRequested;

        /// <summary>Re-scans on the next Update (never inside OnGUI: the list changes the layout).</summary>
        public void RequestPlan()
        {
            if (_planRequested) return;
            _planRequested = true;
            _app.Defer(() =>
            {
                _planRequested = false;
                Plan();
            });
        }

        public void Plan()
        {
            if (_running) return;
            _summary = null;
            if (string.IsNullOrWhiteSpace(_options.InputFolder) && string.IsNullOrWhiteSpace(_options.OutputFolder))
            {
                _items = new List<BatchItem>();
                _planError = null;
                return;
            }
            _items = BatchPlanner.Plan(_options, out _planError);
        }

        private async void Run()
        {
            if (_running) return;
            Plan();
            if (CountPending() == 0) return;
            var settings = _app.Settings.Clone();
            _cts = new CancellationTokenSource();
            var ct = _cts.Token;
            _running = true;
            _processed = 0;
            _toProcess = CountPending();
            _summary = null;
            _clock.Restart();
            int done = 0, failed = 0;
            try
            {
                if (!_app.Projection.Ready) _itemProgress.BeginStage("Preparing the projection", 0f, 0f);
                Directory.CreateDirectory(_options.OutputFolder);
                for (int i = 0; i < _items.Count; i++)
                {
                    var it = _items[i];
                    if (it.Status != BatchItemStatus.Pending) continue;
                    if (ct.IsCancellationRequested) { it.Status = BatchItemStatus.Cancelled; continue; }
                    _current = i;
                    it.Status = BatchItemStatus.Processing;
                    var sw = Stopwatch.StartNew();
                    try
                    {
                        _itemProgress.BeginStage("Loading", 0f, 0.1f);
                        var img = await ImageLoader.LoadAsync(it.InputPath, ct);
                        // the geometry belongs to the frame's shape: made once for each shape in the batch (and kept)
                        var g = await _app.Projection.GeometryAsync(Defisher.SensorAspect(img), settings.Adjust);
                        _itemProgress.BeginStage("Defishing", 0.1f, 0.6f);
                        bool shape = Defisher.ShapeDiffers(g, img);
                        var outImg = await Task.Run(() => Defisher.Process(img, g, settings, _itemProgress, ct), ct);
                        img = null;
                        _itemProgress.BeginStage("Encoding", 0.7f, 0.3f);
                        var fmt = it.Format;
                        string path = it.OutputPath;
                        await Task.Run(() =>
                        {
                            byte[] bytes = ImageIO.Encode(outImg, fmt);
                            ct.ThrowIfCancellationRequested();
                            FileUtil.WriteAllBytesSafely(path, bytes);
                        }, ct);
                        it.Status = BatchItemStatus.Done;
                        it.Message = outImg.Width + " x " + outImg.Height + (shape ? "   not 3:2 - stretched to the projection's frame" : "");
                        done++;
                    }
                    catch (OperationCanceledException)
                    {
                        it.Status = BatchItemStatus.Cancelled;
                        it.Message = "";
                    }
                    catch (Exception e)
                    {
                        it.Status = BatchItemStatus.Failed;
                        it.Message = e.Message;
                        failed++;
                    }
                    it.Seconds = sw.Elapsed.TotalSeconds;
                    _processed++;
                    _itemProgress.Reset();
                    GC.Collect(); // large buffers: keep the heap from growing across many images
                }
            }
            catch (Exception e)
            {
                _app.ShowError("Batch stopped: " + e.Message);
            }
            finally
            {
                _running = false;
                _current = -1;
                _cts = null;
                _clock.Stop();
            }
            int skipped = 0, cancelled = 0;
            foreach (var it in _items)
            {
                if (it.Status == BatchItemStatus.Skipped) skipped++;
                if (it.Status == BatchItemStatus.Cancelled) cancelled++;
            }
            _summary = done + " done, " + failed + " failed, " + skipped + " skipped" + (cancelled > 0 ? ", " + cancelled + " cancelled" : "") +
                       " in " + _clock.Elapsed.TotalSeconds.ToString("0.0", CultureInfo.InvariantCulture) + " s";
            _app.ShowInfo("Batch finished: " + _summary);
        }

        public void Dispose()
        {
            if (_cts != null) _cts.Cancel();
        }
    }
}
