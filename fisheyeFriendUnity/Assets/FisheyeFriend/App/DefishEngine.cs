using System;
using System.Threading;
using System.Threading.Tasks;
using FisheyeFriend.Core;
using FisheyeFriend.Core.Defish;
using FisheyeFriend.Core.Imaging;
using FisheyeFriend.Gpu;
using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>
    /// Where defishing happens: on the GPU when the user has it on and it can do the job, else on the CPU. Both give the
    /// same picture (see ConformanceTests: equal but for the odd value on a rounding boundary, and then by one level).
    /// The geometry is the same exact one either way.
    /// </summary>
    internal sealed class DefishEngine
    {
        private readonly GpuDefisher _gpu;
        private bool _enabled;

        public DefishEngine(bool enabled)
        {
            string why;
            try
            {
                _gpu = GpuDefisher.Create(out why);
            }
            catch (Exception e)
            {
                _gpu = null;
                why = e.Message;
            }
            Unavailable = _gpu == null ? why : null;
            _enabled = enabled;
        }

        /// <summary>Why the GPU can't be used on this computer, or null when it can.</summary>
        public string Unavailable { get; }
        public bool Available => _gpu != null;

        /// <summary>The user's choice (remembered); the GPU is used when it is on and available.</summary>
        public bool Enabled
        {
            get { return _enabled; }
            set
            {
                _enabled = value;
                Prefs.SetBool("UseGpu", value);
            }
        }

        public bool UsesGpu => _enabled && _gpu != null;

        public sealed class Outcome
        {
            public PixelImage Image;
            public bool OnGpu;
            /// <summary>Why the CPU did it although the GPU is on, or null.</summary>
            public string Note;
        }

        /// <summary>Defishes an upright photo. Start it on Unity's main thread (the GPU is driven from there).</summary>
        public async Task<Outcome> DefishAsync(PixelImage upright, DefishGeometry g, DefishSettings settings, ProgressTracker progress, CancellationToken ct)
        {
            string note = null;
            if (UsesGpu && !g.Folds)           // folding settings: the CPU path reports them
            {
                string why;
                if (GpuDefisher.CanProcess(upright, settings, out why))
                {
                    try
                    {
                        var r = await _gpu.ProcessAsync(upright, g, settings, progress, ct);
                        return new Outcome { Image = r, OnGpu = true };
                    }
                    catch (OperationCanceledException)
                    {
                        throw;
                    }
                    catch (Exception e)
                    {
                        Debug.LogException(e);
                        note = "the GPU failed (" + e.Message + ")";
                    }
                }
                else
                {
                    note = why;
                }
            }
            var cpu = await Task.Run(() => Defisher.Process(upright, g, settings, progress, ct), ct);
            return new Outcome { Image = cpu, OnGpu = false, Note = note };
        }

        /// <summary>A few words for the status line: where it ran, and why not on the GPU if it was on.</summary>
        public static string Where(Outcome o)
        {
            if (o.OnGpu) return "on the GPU";
            return o.Note != null ? "on the CPU, as " + o.Note : "on the CPU";
        }
    }
}
