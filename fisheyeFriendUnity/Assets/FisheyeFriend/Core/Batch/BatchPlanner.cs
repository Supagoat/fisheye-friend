using System;
using System.Collections.Generic;
using System.IO;
using FisheyeFriend.Core.Codecs;
using FisheyeFriend.Core.Imaging;

namespace FisheyeFriend.Core.Batch
{
    public enum BatchItemStatus
    {
        Pending,
        Processing,
        Done,
        Skipped,
        Failed,
        Cancelled,
    }

    public sealed class BatchItem
    {
        public string InputPath;
        public string OutputPath;
        /// <summary>Path relative to the input folder, for display.</summary>
        public string RelativePath;
        public ImageFileFormat Format;
        public BatchItemStatus Status = BatchItemStatus.Pending;
        public string Message = "";
        public double Seconds;
    }

    [Serializable]
    public sealed class BatchOptions
    {
        public string InputFolder = "";
        public string OutputFolder = "";
        public bool IncludeSubfolders;
        /// <summary>Appended to each file name before the extension, e.g. "_desqueezed".</summary>
        public string FileSuffix = "";
        public bool OverwriteExisting;
    }

    /// <summary>Finds the images to process and decides where each result goes. Output keeps the input's format.</summary>
    public static class BatchPlanner
    {
        public static List<BatchItem> Plan(BatchOptions o, out string error)
        {
            error = null;
            var items = new List<BatchItem>();
            if (o == null || string.IsNullOrWhiteSpace(o.InputFolder))
            {
                error = "Choose an input folder.";
                return items;
            }
            if (string.IsNullOrWhiteSpace(o.OutputFolder))
            {
                error = "Choose an output folder.";
                return items;
            }
            string input, output;
            try
            {
                input = Normalize(o.InputFolder);
                output = Normalize(o.OutputFolder);
            }
            catch (Exception e)
            {
                error = "Invalid folder path: " + e.Message;
                return items;
            }
            if (!Directory.Exists(input))
            {
                error = "Input folder does not exist: " + input;
                return items;
            }
            string suffix = o.FileSuffix ?? "";
            if (suffix.IndexOfAny(Path.GetInvalidFileNameChars()) >= 0)
            {
                error = "The file suffix contains characters that can't be used in file names.";
                return items;
            }
            bool sameFolder = string.Equals(input, output, PathComparison);
            if (sameFolder && suffix.Length == 0)
            {
                error = "Output folder is the input folder: add a file suffix so the originals are not overwritten.";
                return items;
            }

            IEnumerable<string> files;
            try
            {
                files = Directory.EnumerateFiles(input, "*",
                    o.IncludeSubfolders ? SearchOption.AllDirectories : SearchOption.TopDirectoryOnly);
            }
            catch (Exception e)
            {
                error = "Could not read the input folder: " + e.Message;
                return items;
            }

            bool outputInsideInput = !sameFolder && IsInside(output, input);
            var list = new List<string>();
            foreach (var f in files)
            {
                string name = Path.GetFileName(f);
                if (name.StartsWith(".", StringComparison.Ordinal)) continue; // hidden and macOS "._" resource forks
                if (!ImageIO.IsSupportedPath(f)) continue;
                if (outputInsideInput && IsInside(Path.GetDirectoryName(f), output)) continue; // don't re-process results
                list.Add(f);
            }
            list.Sort(StringComparer.OrdinalIgnoreCase);

            var planned = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (var f in list)
            {
                string rel = GetRelativePath(input, f);
                string relDir = Path.GetDirectoryName(rel) ?? "";
                string ext = Path.GetExtension(f);
                string outName = Path.GetFileNameWithoutExtension(f) + suffix + ext;
                string outPath = Path.Combine(output, relDir, outName);
                var item = new BatchItem
                {
                    InputPath = f,
                    OutputPath = outPath,
                    RelativePath = rel,
                    Format = ImageIO.FormatFromExtension(f),
                };
                if (string.Equals(Normalize(outPath), Normalize(f), PathComparison))
                {
                    item.Status = BatchItemStatus.Skipped;
                    item.Message = "Output would overwrite the original.";
                }
                else if (!planned.Add(outPath))
                {
                    item.Status = BatchItemStatus.Skipped;
                    item.Message = "Another file already writes to this output name.";
                }
                else if (!o.OverwriteExisting && File.Exists(outPath))
                {
                    item.Status = BatchItemStatus.Skipped;
                    item.Message = "Output already exists (turn on Overwrite to replace).";
                }
                items.Add(item);
            }
            if (items.Count == 0) error = "No .jpg, .jpeg or .png files found in " + input;
            return items;
        }

        private static StringComparison PathComparison =>
            Path.DirectorySeparatorChar == '\\' ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal;

        private static string Normalize(string p)
        {
            string full = Path.GetFullPath(p.Trim());
            return full.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
        }

        private static bool IsInside(string path, string folder)
        {
            if (path == null) return false;
            string p = Normalize(path) + Path.DirectorySeparatorChar;
            string f = Normalize(folder) + Path.DirectorySeparatorChar;
            return p.StartsWith(f, PathComparison);
        }

        private static string GetRelativePath(string folder, string file)
        {
            string f = Normalize(folder) + Path.DirectorySeparatorChar;
            string full = Path.GetFullPath(file);
            return full.StartsWith(f, PathComparison) ? full.Substring(f.Length) : Path.GetFileName(file);
        }
    }
}
