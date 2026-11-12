using System;
using System.Globalization;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Batch;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Resampling;
using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>Remembers the user's last settings and folders between sessions (PlayerPrefs).</summary>
    internal static class Prefs
    {
        private const string P = "FisheyeFriend.";

        /// <summary>Every key this app writes (without the prefix), for "Forget Remembered Settings".</summary>
        public static readonly string[] Keys =
        {
            "OutputSize", "Edges", "Kernel", "Radius", "Light", "Antialias", "Aperture", "DefishOnOpen",
            "SaveFormat", "OpenDir", "SaveDir", "Mode", "BatchIn", "BatchOut", "BatchSubfolders", "BatchOverwrite", "BatchSuffix",
        };

        public static string GetString(string key, string fallback)
        {
            return PlayerPrefs.GetString(P + key, fallback);
        }

        public static void SetString(string key, string value)
        {
            PlayerPrefs.SetString(P + key, value ?? "");
        }

        public static bool GetBool(string key, bool fallback)
        {
            return PlayerPrefs.GetInt(P + key, fallback ? 1 : 0) != 0;
        }

        public static void SetBool(string key, bool value)
        {
            PlayerPrefs.SetInt(P + key, value ? 1 : 0);
        }

        public static double GetDouble(string key, double fallback)
        {
            double v;
            string s = PlayerPrefs.GetString(P + key, "");
            return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v) ? v : fallback;
        }

        public static void SetDouble(string key, double value)
        {
            PlayerPrefs.SetString(P + key, value.ToString("R", CultureInfo.InvariantCulture));
        }

        private static T GetEnum<T>(string key, T fallback) where T : struct
        {
            T v;
            if (Enum.TryParse(GetString(key, fallback.ToString()), out v) && Enum.IsDefined(typeof(T), v)) return v;
            return fallback;
        }

        public static DefishSettings LoadSettings(DefishSettings defaults)
        {
            var s = defaults.Clone();
            s.OutputSize = GetEnum("OutputSize", s.OutputSize);
            s.Edges = GetEnum("Edges", s.Edges);
            s.Kernel = GetEnum("Kernel", s.Kernel);
            s.Radius = GetDouble("Radius", s.Radius);
            s.Light = GetEnum("Light", s.Light);
            s.Antialias = GetBool("Antialias", s.Antialias);
            s.ApertureCorrection = GetBool("Aperture", s.ApertureCorrection);
            s.Validate();
            return s;
        }

        public static void SaveSettings(DefishSettings s)
        {
            SetString("OutputSize", s.OutputSize.ToString());
            SetString("Edges", s.Edges.ToString());
            SetString("Kernel", s.Kernel.ToString());
            SetDouble("Radius", s.Radius);
            SetString("Light", s.Light.ToString());
            SetBool("Antialias", s.Antialias);
            SetBool("Aperture", s.ApertureCorrection);
        }

        public static void LoadBatch(BatchOptions o, string defaultSuffix)
        {
            o.InputFolder = GetString("BatchIn", "");
            o.OutputFolder = GetString("BatchOut", "");
            o.IncludeSubfolders = GetBool("BatchSubfolders", false);
            o.OverwriteExisting = GetBool("BatchOverwrite", false);
            o.FileSuffix = GetString("BatchSuffix", defaultSuffix);
        }

        public static void SaveBatch(BatchOptions o)
        {
            SetString("BatchIn", o.InputFolder);
            SetString("BatchOut", o.OutputFolder);
            SetBool("BatchSubfolders", o.IncludeSubfolders);
            SetBool("BatchOverwrite", o.OverwriteExisting);
            SetString("BatchSuffix", o.FileSuffix);
        }

        public static void Flush()
        {
            PlayerPrefs.Save();
        }
    }
}
