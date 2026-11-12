using System.Collections.Generic;
using UnityEngine;

namespace FisheyeFriend.App.UI
{
    /// <summary>
    /// Dark theme for the IMGUI interface. Sizes are multiplied by <see cref="Scale"/> (font sizes and
    /// metrics, not a GUI matrix) so text stays crisp on high-DPI screens.
    /// Must be built from inside OnGUI.
    /// </summary>
    internal sealed class UiStyles
    {
        public static readonly Color Background = Hex(0x17181B);
        public static readonly Color PanelColor = Hex(0x222428);
        public static readonly Color PanelBorder = Hex(0x33363C);
        public static readonly Color Text = Hex(0xE4E6EA);
        public static readonly Color Muted = Hex(0x979CA6);
        public static readonly Color Accent = Hex(0xE8A13A);      // warm amber, shared with the de-squeeze app
        public static readonly Color AccentDim = Hex(0x8A6224);
        public static readonly Color Good = Hex(0x6CC18A);
        public static readonly Color Bad = Hex(0xE5675F);
        public static readonly Color ViewerBackground = Hex(0x101113);

        public float Scale { get; private set; }

        /// <summary>
        /// Screen pixels per IMGUI unit. In players and the Game view, IMGUI works directly in screen
        /// pixels (Screen.width is in pixels, Retina included), so this is 1. All viewer maths goes
        /// through this one value so the pixel-perfect drawing stays correct if that ever changes.
        /// </summary>
        public static float PixelsPerPoint => 1f;

        public GUIStyle Title, Section, Label, Small, MutedSmall, Value, Error, Success;
        public GUIStyle Button, ButtonSelected, Primary, Tab, TabSelected, TextField, Slider, Thumb;
        public GUIStyle Panel, Bar, ListRow, ListRowAlt, CheckLabel, Tooltip;
        public Texture2D White, Checker;

        private readonly List<Texture2D> _owned = new List<Texture2D>();
        private Font _font;

        public float S(float v)
        {
            return Mathf.Round(v * Scale);
        }

        public bool IsBuiltFor(float scale)
        {
            return Title != null && Mathf.Approximately(scale, Scale);
        }

        public void Build(float scale)
        {
            Release();
            Scale = scale;
            _font = GUI.skin.font;

            White = Tex(Color.white);
            Checker = MakeChecker();

            Title = Text0(17, FontStyle.Bold, Text);
            Section = Text0(11, FontStyle.Bold, Muted);
            Section.margin = new RectOffset(0, 0, (int)S(12), (int)S(4));
            Label = Text0(13, FontStyle.Normal, Text);
            Small = Text0(11, FontStyle.Normal, Text);
            Small.wordWrap = true;
            MutedSmall = Text0(11, FontStyle.Normal, Muted);
            MutedSmall.wordWrap = true;
            Value = Text0(13, FontStyle.Bold, Accent);
            Error = Text0(12, FontStyle.Normal, Bad);
            Error.wordWrap = true;
            Success = Text0(12, FontStyle.Normal, Good);
            Success.wordWrap = true;
            CheckLabel = Text0(13, FontStyle.Normal, Text);
            CheckLabel.padding = new RectOffset((int)S(24), 0, (int)S(3), (int)S(3));

            Button = MakeButton(Hex(0x34373D), Hex(0x41454C), Hex(0x2B2D32), Text);
            ButtonSelected = MakeButton(Accent, Hex(0xF0B55A), AccentDim, Hex(0x1A1A1A));
            ButtonSelected.fontStyle = FontStyle.Bold;
            Primary = MakeButton(Accent, Hex(0xF0B55A), AccentDim, Hex(0x1A1A1A));
            Primary.fontSize = (int)S(14);
            Primary.fontStyle = FontStyle.Bold;
            Primary.fixedHeight = S(38);

            Tab = MakeButton(PanelColor, Hex(0x2C2F34), Hex(0x2C2F34), Muted);
            Tab.fontSize = (int)S(13);
            Tab.fixedHeight = S(30);
            Tab.padding = new RectOffset((int)S(16), (int)S(16), (int)S(6), (int)S(6));
            TabSelected = MakeButton(Hex(0x34373D), Hex(0x34373D), Hex(0x34373D), Text);
            TabSelected.fontSize = Tab.fontSize;
            TabSelected.fontStyle = FontStyle.Bold;
            TabSelected.fixedHeight = Tab.fixedHeight;
            TabSelected.padding = Tab.padding;

            TextField = new GUIStyle(GUI.skin.textField)
            {
                font = _font,
                fontSize = (int)S(13),
                padding = new RectOffset((int)S(6), (int)S(6), (int)S(4), (int)S(4)),
                margin = new RectOffset((int)S(2), (int)S(2), (int)S(3), (int)S(3)),
                fixedHeight = S(26),
                clipping = TextClipping.Clip,
            };
            var fieldBg = Tex(Hex(0x131416));
            var fieldFocus = Tex(Hex(0x0E0F11));
            TextField.normal.background = fieldBg;
            TextField.hover.background = fieldBg;
            TextField.focused.background = fieldFocus;
            TextField.active.background = fieldFocus;
            TextField.onNormal.background = fieldBg;
            TextField.onFocused.background = fieldFocus;
            TextField.normal.textColor = TextField.hover.textColor = TextField.focused.textColor = Text;
            TextField.border = new RectOffset(0, 0, 0, 0);

            Slider = new GUIStyle(GUI.skin.horizontalSlider)
            {
                fixedHeight = S(18),
                margin = new RectOffset((int)S(6), (int)S(6), (int)S(6), (int)S(6)),
                padding = new RectOffset(0, 0, 0, 0),
                border = new RectOffset(0, 0, 0, 0),
                overflow = new RectOffset(0, 0, (int)-S(7), (int)-S(7)),
            };
            Slider.normal.background = Tex(Hex(0x3A3D43));
            Thumb = new GUIStyle(GUI.skin.horizontalSliderThumb)
            {
                fixedWidth = S(14),
                fixedHeight = S(18),
                border = new RectOffset(0, 0, 0, 0),
                padding = new RectOffset(0, 0, 0, 0),
                margin = new RectOffset(0, 0, 0, 0),
                overflow = new RectOffset(0, 0, 0, 0),
            };
            Thumb.normal.background = Tex(Accent);
            Thumb.hover.background = Tex(Hex(0xF0B55A));
            Thumb.active.background = Tex(Hex(0xFFC870));
            Thumb.focused.background = Thumb.normal.background;

            Panel = new GUIStyle
            {
                padding = new RectOffset((int)S(14), (int)S(14), (int)S(10), (int)S(10)),
            };
            Panel.normal.background = Tex(PanelColor);
            Bar = new GUIStyle
            {
                padding = new RectOffset((int)S(10), (int)S(10), (int)S(5), (int)S(5)),
            };
            Bar.normal.background = Tex(Hex(0x1D1F22));
            ListRow = Text0(12, FontStyle.Normal, Text);
            ListRow.padding = new RectOffset((int)S(6), (int)S(6), (int)S(3), (int)S(3));
            ListRow.clipping = TextClipping.Clip;
            ListRow.margin = new RectOffset(0, 0, 0, 0);
            ListRowAlt = new GUIStyle(ListRow);
            ListRowAlt.normal.background = Tex(Hex(0x26282D));
            Tooltip = Text0(11, FontStyle.Normal, Text);
            Tooltip.wordWrap = true;
            Tooltip.padding = new RectOffset((int)S(8), (int)S(8), (int)S(6), (int)S(6));
            Tooltip.normal.background = Tex(Hex(0x3A3D43));
        }

        public void Release()
        {
            foreach (var t in _owned)
                if (t != null) Object.Destroy(t);
            _owned.Clear();
            Title = null;
        }

        private GUIStyle Text0(int size, FontStyle style, Color color)
        {
            var s = new GUIStyle(GUI.skin.label)
            {
                font = _font,
                fontSize = (int)S(size),
                fontStyle = style,
                wordWrap = false,
                richText = false,
                padding = new RectOffset((int)S(2), (int)S(2), (int)S(2), (int)S(2)),
                margin = new RectOffset((int)S(2), (int)S(2), (int)S(2), (int)S(2)),
            };
            s.normal.textColor = color;
            s.normal.background = null;
            return s;
        }

        private GUIStyle MakeButton(Color normal, Color hover, Color active, Color text)
        {
            var s = new GUIStyle(GUI.skin.button)
            {
                font = _font,
                fontSize = (int)S(13),
                fixedHeight = S(28),
                padding = new RectOffset((int)S(10), (int)S(10), (int)S(4), (int)S(4)),
                margin = new RectOffset((int)S(2), (int)S(2), (int)S(2), (int)S(2)),
                border = new RectOffset(0, 0, 0, 0),
                overflow = new RectOffset(0, 0, 0, 0),
                alignment = TextAnchor.MiddleCenter,
                clipping = TextClipping.Clip,
            };
            var n = Tex(normal);
            var h = Tex(hover);
            var a = Tex(active);
            s.normal.background = n;
            s.hover.background = h;
            s.active.background = a;
            s.focused.background = n;
            s.onNormal.background = n;
            s.onHover.background = h;
            s.onActive.background = a;
            s.normal.textColor = s.hover.textColor = s.active.textColor = s.focused.textColor = text;
            s.onNormal.textColor = s.onHover.textColor = s.onActive.textColor = text;
            return s;
        }

        private Texture2D Tex(Color c)
        {
            var t = new Texture2D(1, 1, TextureFormat.RGBA32, false) { hideFlags = HideFlags.HideAndDontSave };
            t.SetPixel(0, 0, c);
            t.Apply();
            _owned.Add(t);
            return t;
        }

        private Texture2D MakeChecker()
        {
            var t = new Texture2D(2, 2, TextureFormat.RGBA32, false)
            {
                hideFlags = HideFlags.HideAndDontSave,
                filterMode = FilterMode.Point,
                wrapMode = TextureWrapMode.Repeat,
            };
            Color a = Hex(0x3A3A3A), b = Hex(0x2A2A2A);
            t.SetPixels(new[] { a, b, b, a });
            t.Apply();
            _owned.Add(t);
            return t;
        }

        public static Color Hex(int rgb)
        {
            return new Color(((rgb >> 16) & 0xFF) / 255f, ((rgb >> 8) & 0xFF) / 255f, (rgb & 0xFF) / 255f, 1f);
        }

        // ---- small widgets ------------------------------------------------------------------------

        /// <summary>A checkbox drawn to scale (the default IMGUI toggle graphic does not scale).</summary>
        public bool Toggle(bool value, GUIContent content)
        {
            Rect r = GUILayoutUtility.GetRect(content, CheckLabel);
            float box = S(14);
            var boxRect = new Rect(r.x + S(2), r.y + (r.height - box) * 0.5f, box, box);
            if (Event.current.type == EventType.Repaint)
            {
                var old = GUI.color;
                GUI.color = value ? Accent : Hex(0x4A4E55);
                GUI.DrawTexture(boxRect, White);
                if (!value)
                {
                    GUI.color = Hex(0x1A1B1E);
                    GUI.DrawTexture(new Rect(boxRect.x + S(2), boxRect.y + S(2), box - S(4), box - S(4)), White);
                }
                else
                {
                    GUI.color = Hex(0x1A1A1A);
                    GUI.DrawTexture(new Rect(boxRect.x + S(4), boxRect.y + S(4), box - S(8), box - S(8)), White);
                }
                GUI.color = old;
            }
            if (GUI.Button(r, content, CheckLabel)) value = !value;
            return value;
        }

        /// <summary>Horizontal progress bar.</summary>
        public void ProgressBar(float value, string text)
        {
            Rect r = GUILayoutUtility.GetRect(S(10), S(20), GUILayout.ExpandWidth(true));
            if (Event.current.type != EventType.Repaint) return;
            var old = GUI.color;
            GUI.color = Hex(0x131416);
            GUI.DrawTexture(r, White);
            GUI.color = Accent;
            GUI.DrawTexture(new Rect(r.x, r.y, r.width * Mathf.Clamp01(value), r.height), White);
            GUI.color = old;
            var style = new GUIStyle(Small) { alignment = TextAnchor.MiddleCenter };
            style.normal.textColor = Color.white;
            GUI.Label(r, text, style);
        }

        public void Separator()
        {
            Rect r = GUILayoutUtility.GetRect(S(10), S(9), GUILayout.ExpandWidth(true));
            if (Event.current.type != EventType.Repaint) return;
            var old = GUI.color;
            GUI.color = PanelBorder;
            GUI.DrawTexture(new Rect(r.x, r.y + S(4), r.width, Mathf.Max(1, S(1))), White);
            GUI.color = old;
        }

        public void Fill(Rect r, Color c)
        {
            if (Event.current.type != EventType.Repaint) return;
            var old = GUI.color;
            GUI.color = c;
            GUI.DrawTexture(r, White);
            GUI.color = old;
        }
    }
}
