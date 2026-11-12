using System;
using System.IO;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Resampling;
using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>
    /// App configuration, read from StreamingAssets/fisheye-friend-config.json at start-up.
    /// Edit that file to change the defaults or the projection; no rebuild is needed for a built player (the file sits
    /// next to the app's data). Any field left out keeps the value below.
    /// </summary>
    [Serializable]
    public sealed class AppConfig
    {
        public const string FileName = "fisheye-friend-config.json";

        /// <summary>The projection: a point-pair file in Fisheye Friend's format, in StreamingAssets (or an absolute path).</summary>
        public string projectionFile = "defish-hemi-596.json";
        /// <summary>Name shown for the projection.</summary>
        public string projectionName = "Hemi";

        /// <summary>"SameAsPhoto" or "FullDetail" (about 1.19x larger, so the centre keeps every pixel).</summary>
        public string outputSize = "SameAsPhoto";
        /// <summary>"Extend" (keep the full frame) or "Crop" (to the area the photo covers).</summary>
        public string edges = "Extend";
        /// <summary>Lanczos lobes for the warp: 4 (best), 3, or 2 for bicubic.</summary>
        public int lanczosLobes = 4;
        /// <summary>"Sigmoidal", "Gamma" or "Linear".</summary>
        public string lightSpace = "Sigmoidal";
        public bool antialias = false;
        public bool apertureCorrection = true;

        /// <summary>Defish each photo as soon as it is opened in interactive mode.</summary>
        public bool defishOnOpen = true;

        /// <summary>"jpg" or "png" for interactive Save As.</summary>
        public string defaultSaveFormat = "jpg";
        public string batchFileSuffix = "_defished";

        /// <summary>UI size multiplier. 0 picks automatically from the screen's DPI.</summary>
        public float uiScale = 0f;

        /// <summary>Largest viewer magnification (integer).</summary>
        public int maxZoom = 64;

        public static string FilePath => Path.Combine(Application.streamingAssetsPath, FileName);

        public string ProjectionPath =>
            Path.IsPathRooted(projectionFile) ? projectionFile : Path.Combine(Application.streamingAssetsPath, projectionFile);

        /// <summary>Loads the config file. Never throws: problems are reported in <paramref name="warning"/>.</summary>
        public static AppConfig Load(out string warning)
        {
            warning = null;
            var cfg = new AppConfig();
            try
            {
                string path = FilePath;
                if (File.Exists(path)) JsonUtility.FromJsonOverwrite(File.ReadAllText(path), cfg);
                else warning = "Config file not found (" + path + "); using built-in defaults.";
            }
            catch (Exception e)
            {
                cfg = new AppConfig();
                warning = "Could not read " + FileName + ": " + e.Message + " - using built-in defaults.";
            }
            cfg.Sanitize();
            return cfg;
        }

        private void Sanitize()
        {
            if (string.IsNullOrWhiteSpace(projectionFile)) projectionFile = "defish-hemi-596.json";
            if (string.IsNullOrWhiteSpace(projectionName)) projectionName = Path.GetFileNameWithoutExtension(projectionFile);
            if (lanczosLobes < 2 || lanczosLobes > 4) lanczosLobes = 4;
            if (maxZoom < 1) maxZoom = 1;
            if (maxZoom > 256) maxZoom = 256;
            if (uiScale < 0 || uiScale > 4) uiScale = 0;
            defaultSaveFormat = (defaultSaveFormat ?? "jpg").Trim().TrimStart('.').ToLowerInvariant();
            if (defaultSaveFormat != "png") defaultSaveFormat = "jpg";
            if (batchFileSuffix == null) batchFileSuffix = "";
        }

        private static T ParseEnum<T>(string value, T fallback) where T : struct
        {
            T v;
            if (!string.IsNullOrEmpty(value) && Enum.TryParse(value.Trim(), true, out v) && Enum.IsDefined(typeof(T), v)) return v;
            return fallback;
        }

        public DefishSettings CreateDefaultSettings()
        {
            var s = new DefishSettings
            {
                OutputSize = ParseEnum(outputSize, DefishOutputSize.SameAsPhoto),
                Edges = ParseEnum(edges, DefishEdges.Extend),
                Light = ParseEnum(lightSpace, LightSpace.Sigmoidal),
                Antialias = antialias,
                ApertureCorrection = apertureCorrection,
            };
            if (lanczosLobes == 2)
            {
                s.Kernel = WarpKernel.Bicubic;
                s.Radius = 2;
            }
            else
            {
                s.Kernel = WarpKernel.TensorLanczos;
                s.Radius = lanczosLobes;
            }
            s.Validate();
            return s;
        }
    }
}
