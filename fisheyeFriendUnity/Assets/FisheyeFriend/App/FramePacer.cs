using UnityEngine;

namespace FisheyeFriend.App
{
    /// <summary>
    /// Keeps the frame rate as low as the moment allows, to save battery. The whole interface is redrawn every frame, so
    /// frames are the app's main cost when nothing is happening. vSync is turned off, else Unity ignores the target and
    /// draws at the display's refresh rate (60-120 Hz) all the time.
    ///
    ///   30 fps  while the user is using the app (any input in the last second, or a drag in progress)
    ///   15 fps  otherwise, while the window is in front or a job is running
    ///    5 fps  in the background with nothing to do
    ///
    /// Work is done on worker threads and doesn't wait on frames, so a lower rate only makes the screen update less
    /// often, never makes a job slower.
    /// </summary>
    internal sealed class FramePacer
    {
        public const int Interacting = 30, Idle = 15, Background = 5;
        private const float InteractingSeconds = 1f;

        private float _lastInput = -100f;
        private Vector2 _lastMouse = new Vector2(float.NaN, float.NaN);
        private int _rate;

        public void Start()
        {
            QualitySettings.vSyncCount = 0;
            _lastInput = Time.realtimeSinceStartup;
            Apply(Interacting);
        }

        /// <summary>Notes user input; call at the start of OnGUI.</summary>
        public void Observe(Event e)
        {
            switch (e.type)
            {
                case EventType.MouseDown:
                case EventType.MouseUp:
                case EventType.MouseDrag:
                case EventType.MouseMove:
                case EventType.ScrollWheel:
                case EventType.KeyDown:
                case EventType.KeyUp:
                    _lastInput = Time.realtimeSinceStartup;
                    break;
                case EventType.Repaint:
                    // the pointer moving over the window (not every platform sends MouseMove events)
                    if (e.mousePosition != _lastMouse)
                    {
                        if (!float.IsNaN(_lastMouse.x)) _lastInput = Time.realtimeSinceStartup;
                        _lastMouse = e.mousePosition;
                    }
                    break;
            }
        }

        /// <summary>Sets the frame rate for the moment; call from Update. <paramref name="busy"/>: a job is running.</summary>
        public void Tick(bool busy)
        {
            bool focused = Application.isFocused;
            bool interacting = focused && (Time.realtimeSinceStartup - _lastInput < InteractingSeconds || GUIUtility.hotControl != 0);
            int rate = interacting ? Interacting : focused || busy ? Idle : Background;
            if (rate != _rate) Apply(rate);
        }

        private void Apply(int rate)
        {
            _rate = rate;
            Application.targetFrameRate = rate;
        }
    }
}
