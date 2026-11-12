using System.IO;
using UnityEditor;
using UnityEngine;

namespace FisheyeFriend.EditorTools
{
    /// <summary>
    /// Editor conveniences: applies desktop-app player settings once per project (windowed, resizable,
    /// keeps running in the background), and adds a Fisheye Friend menu.
    /// </summary>
    [InitializeOnLoad]
    internal static class ProjectSetup
    {
        private const string AppliedKey = "FisheyeFriend.PlayerSettingsApplied.v1";

        static ProjectSetup()
        {
            EditorApplication.delayCall += ApplyOnce;
        }

        private static void ApplyOnce()
        {
            if (EditorUserSettings.GetConfigValue(AppliedKey) == "1") return;
            ApplyPlayerSettings();
            EditorUserSettings.SetConfigValue(AppliedKey, "1");
        }

        [MenuItem("Fisheye Friend/Apply Recommended Player Settings")]
        private static void ApplyPlayerSettings()
        {
            PlayerSettings.fullScreenMode = FullScreenMode.Windowed;
            PlayerSettings.resizableWindow = true;
            PlayerSettings.defaultScreenWidth = 1600;
            PlayerSettings.defaultScreenHeight = 1000;
            PlayerSettings.runInBackground = true;
            PlayerSettings.macRetinaSupport = true;
            AssetDatabase.SaveAssets();
            Debug.Log("[Fisheye Friend] Player settings applied: windowed, resizable, 1600x1000, runs in background, Retina on. " +
                      "(Menu: Fisheye Friend > Apply Recommended Player Settings)");
        }

        private static string ConfigPath => Path.Combine(Application.streamingAssetsPath, "fisheye-friend-config.json");

        [MenuItem("Fisheye Friend/Open Config File")]
        private static void OpenConfig()
        {
            if (File.Exists(ConfigPath)) EditorUtility.OpenWithDefaultApp(ConfigPath);
            else Debug.LogWarning("[Fisheye Friend] Config file not found: " + ConfigPath);
        }

        [MenuItem("Fisheye Friend/Reveal Config File")]
        private static void RevealConfig()
        {
            EditorUtility.RevealInFinder(ConfigPath);
        }

        [MenuItem("Fisheye Friend/Forget Remembered Settings")]
        private static void ClearPrefs()
        {
            string[] keys =
            {
                "OutputSize", "Edges", "Kernel", "Radius", "Light", "Antialias", "Aperture", "DefishOnOpen",
                "SaveFormat", "OpenDir", "SaveDir", "Mode", "BatchIn", "BatchOut", "BatchSubfolders", "BatchOverwrite", "BatchSuffix",
            };
            foreach (var k in keys) PlayerPrefs.DeleteKey("FisheyeFriend." + k);
            PlayerPrefs.Save();
            Debug.Log("[Fisheye Friend] Remembered settings cleared; the config file defaults apply on next launch.");
        }
    }
}
