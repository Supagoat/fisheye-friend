using System;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Threading.Tasks;

namespace FisheyeFriend.App.IO
{
    /// <summary>
    /// Native open/save/folder dialogs without plugins.
    ///  * Unity Editor: EditorUtility panels.
    ///  * macOS player: AppleScript via osascript (standard Finder dialogs).
    ///  * Windows player: WinForms dialogs via PowerShell.
    ///  * Linux player: zenity (or kdialog).
    /// External dialogs run in a separate process, so the app keeps drawing while they are open.
    /// Every method returns null when the user cancels or no dialog is available; the UI also has
    /// editable path fields as a fallback.
    /// </summary>
    internal static class FileDialogs
    {
        public enum SaveFormat
        {
            Jpeg,
            Png,
        }

        public static Task<string> OpenImageAsync(string title, string startDir)
        {
#if UNITY_EDITOR
            string p = UnityEditor.EditorUtility.OpenFilePanelWithFilters(title, ExistingDir(startDir) ?? "",
                new[] { "Images", "jpg,jpeg,png,JPG,JPEG,PNG" });
            return Task.FromResult(string.IsNullOrEmpty(p) ? null : p);
#else
            return RunExternalAsync("open", title, ExistingDir(startDir), null, SaveFormat.Jpeg);
#endif
        }

        public static Task<string> PickFolderAsync(string title, string startDir)
        {
#if UNITY_EDITOR
            string p = UnityEditor.EditorUtility.OpenFolderPanel(title, ExistingDir(startDir) ?? "", "");
            return Task.FromResult(string.IsNullOrEmpty(p) ? null : p);
#else
            return RunExternalAsync("folder", title, ExistingDir(startDir), null, SaveFormat.Jpeg);
#endif
        }

        public static Task<string> SaveImageAsync(string title, string startDir, string defaultName, SaveFormat format)
        {
#if UNITY_EDITOR
            string ext = format == SaveFormat.Png ? "png" : "jpg";
            string p = UnityEditor.EditorUtility.SaveFilePanel(title, ExistingDir(startDir) ?? "",
                Path.GetFileNameWithoutExtension(defaultName), ext);
            return Task.FromResult(string.IsNullOrEmpty(p) ? null : p);
#else
            return RunExternalAsync("save", title, ExistingDir(startDir), defaultName, format);
#endif
        }

        /// <summary>True when a native dialog can be shown on this platform.</summary>
        public static bool Available
        {
            get
            {
#if UNITY_EDITOR || UNITY_STANDALONE_OSX || UNITY_STANDALONE_WIN || UNITY_STANDALONE_LINUX
                return true;
#else
                return false;
#endif
            }
        }

        private static string ExistingDir(string dir)
        {
            try
            {
                if (string.IsNullOrEmpty(dir)) return null;
                if (File.Exists(dir)) dir = Path.GetDirectoryName(dir);
                return Directory.Exists(dir) ? Path.GetFullPath(dir) : null;
            }
            catch (Exception)
            {
                return null;
            }
        }

        // ---------------------------------------------------------------------------------------------

#pragma warning disable 0162, 0414 // some branches are unreachable on a given platform
        private static Task<string> RunExternalAsync(string kind, string title, string dir, string name, SaveFormat format)
        {
            return Task.Run(() =>
            {
                try
                {
#if UNITY_STANDALONE_OSX
                    return MacDialog(kind, title, dir, name);
#elif UNITY_STANDALONE_WIN
                    return WindowsDialog(kind, title, dir, name, format);
#elif UNITY_STANDALONE_LINUX
                    return LinuxDialog(kind, title, dir, name, format);
#else
                    return (string)null;
#endif
                }
                catch (Exception e)
                {
                    UnityEngine.Debug.LogWarning("File dialog failed: " + e.Message);
                    return null;
                }
            });
        }

        private const string MacScript = @"
set k to system attribute ""FF_KIND""
set p to system attribute ""FF_PROMPT""
set d to system attribute ""FF_DIR""
set n to system attribute ""FF_NAME""
activate
try
	if k is ""open"" then
		if d is """" then
			set f to choose file with prompt p of type {""public.jpeg"", ""public.png""}
		else
			set f to choose file with prompt p of type {""public.jpeg"", ""public.png""} default location (POSIX file d as alias)
		end if
	else if k is ""folder"" then
		if d is """" then
			set f to choose folder with prompt p
		else
			set f to choose folder with prompt p default location (POSIX file d as alias)
		end if
	else
		if d is """" then
			set f to choose file name with prompt p default name n
		else
			set f to choose file name with prompt p default name n default location (POSIX file d as alias)
		end if
	end if
	return POSIX path of f
on error number -128
	return """"
end try
";

        private static string MacDialog(string kind, string title, string dir, string name)
        {
            var psi = new ProcessStartInfo("/usr/bin/osascript", "-")
            {
                UseShellExecute = false,
                RedirectStandardInput = true,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                StandardOutputEncoding = Encoding.UTF8,
            };
            // Values travel as environment variables so no escaping is needed in the script.
            psi.EnvironmentVariables["FF_KIND"] = kind;
            psi.EnvironmentVariables["FF_PROMPT"] = title ?? "";
            psi.EnvironmentVariables["FF_DIR"] = dir ?? "";
            psi.EnvironmentVariables["FF_NAME"] = name ?? "";
            return RunProcess(psi, MacScript);
        }

        private static string WindowsDialog(string kind, string title, string dir, string name, SaveFormat format)
        {
            const string script = @"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
Add-Type -AssemblyName System.Windows.Forms
$owner = New-Object System.Windows.Forms.Form -Property @{ TopMost = $true; ShowInTaskbar = $false }
$kind = $env:FF_KIND; $dir = $env:FF_DIR; $prompt = $env:FF_PROMPT
if ($kind -eq 'folder') {
  $d = New-Object System.Windows.Forms.FolderBrowserDialog
  $d.Description = $prompt
  if ($dir) { $d.SelectedPath = $dir }
  if ($d.ShowDialog($owner) -eq 'OK') { [Console]::Out.Write($d.SelectedPath) }
} elseif ($kind -eq 'open') {
  $d = New-Object System.Windows.Forms.OpenFileDialog
  $d.Title = $prompt
  $d.Filter = 'Images (*.jpg;*.jpeg;*.png)|*.jpg;*.jpeg;*.png'
  if ($dir) { $d.InitialDirectory = $dir }
  if ($d.ShowDialog($owner) -eq 'OK') { [Console]::Out.Write($d.FileName) }
} else {
  $d = New-Object System.Windows.Forms.SaveFileDialog
  $d.Title = $prompt
  $d.Filter = 'JPEG (*.jpg)|*.jpg|PNG (*.png)|*.png'; $d.FilterIndex = [int]$env:FF_FILTER
  $d.FileName = $env:FF_NAME
  if ($dir) { $d.InitialDirectory = $dir }
  if ($d.ShowDialog($owner) -eq 'OK') { [Console]::Out.Write($d.FileName) }
}
";
            string encoded = Convert.ToBase64String(Encoding.Unicode.GetBytes(script));
            var psi = new ProcessStartInfo("powershell.exe",
                "-NoProfile -NonInteractive -STA -ExecutionPolicy Bypass -EncodedCommand " + encoded)
            {
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                StandardOutputEncoding = Encoding.UTF8,
            };
            psi.EnvironmentVariables["FF_KIND"] = kind;
            psi.EnvironmentVariables["FF_PROMPT"] = title ?? "";
            psi.EnvironmentVariables["FF_DIR"] = dir ?? "";
            psi.EnvironmentVariables["FF_NAME"] = name ?? "";
            psi.EnvironmentVariables["FF_FILTER"] = format == SaveFormat.Png ? "2" : "1";
            return RunProcess(psi, null);
        }

        private static string LinuxDialog(string kind, string title, string dir, string name, SaveFormat format)
        {
            var sb = new StringBuilder();
            sb.Append("--file-selection --title=").Append(Quote(title ?? ""));
            if (kind == "folder")
            {
                sb.Append(" --directory");
                if (dir != null) sb.Append(" --filename=").Append(Quote(dir + "/"));
            }
            else if (kind == "open")
            {
                sb.Append(" --file-filter=").Append(Quote("Images | *.jpg *.jpeg *.png *.JPG *.JPEG *.PNG"));
                if (dir != null) sb.Append(" --filename=").Append(Quote(dir + "/"));
            }
            else
            {
                sb.Append(" --save --confirm-overwrite");
                string full = dir != null ? Path.Combine(dir, name ?? "") : (name ?? "");
                sb.Append(" --filename=").Append(Quote(full));
            }
            var psi = new ProcessStartInfo("zenity", sb.ToString())
            {
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                StandardOutputEncoding = Encoding.UTF8,
            };
            return RunProcess(psi, null);
        }
#pragma warning restore 0162, 0414

        private static string Quote(string s)
        {
            return "\"" + s.Replace("\\", "\\\\").Replace("\"", "\\\"") + "\"";
        }

        private static string RunProcess(ProcessStartInfo psi, string stdin)
        {
            using (var p = Process.Start(psi))
            {
                if (p == null) return null;
                if (stdin != null)
                {
                    var bytes = Encoding.UTF8.GetBytes(stdin);
                    p.StandardInput.BaseStream.Write(bytes, 0, bytes.Length);
                    p.StandardInput.Close();
                }
                // Read stderr asynchronously so a chatty process cannot block on a full pipe.
                p.ErrorDataReceived += (s, e) => { };
                p.BeginErrorReadLine();
                string output = p.StandardOutput.ReadToEnd();
                p.WaitForExit();
                output = (output ?? "").Trim('\r', '\n', ' ');
                return output.Length == 0 ? null : output;
            }
        }
    }
}
