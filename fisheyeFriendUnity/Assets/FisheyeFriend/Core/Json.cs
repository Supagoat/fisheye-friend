using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace FisheyeFriend.Core
{
    /// <summary>
    /// A minimal JSON reader: objects become <see cref="Dictionary{TKey,TValue}"/> of string to object, arrays
    /// <see cref="List{T}"/> of object, numbers double, plus string, bool and null. Throws <see cref="FormatException"/>.
    /// </summary>
    public static class Json
    {
        public static object Parse(string text)
        {
            var r = new Reader(text ?? "");
            var v = r.Value();
            r.Space();
            if (!r.End) throw r.Error("unexpected text after the end");
            return v;
        }

        public static object Get(object o, string key)
        {
            var d = o as Dictionary<string, object>;
            object v;
            return d != null && d.TryGetValue(key, out v) ? v : null;
        }

        public static double Number(object o, string key, double fallback)
        {
            var v = Get(o, key);
            return v is double ? (double)v : fallback;
        }

        public static string String(object o, string key, string fallback)
        {
            return Get(o, key) as string ?? fallback;
        }

        private sealed class Reader
        {
            private readonly string _s;
            private int _i, _depth;

            public Reader(string s)
            {
                _s = s;
            }

            public bool End => _i >= _s.Length;

            public void Space()
            {
                while (_i < _s.Length && char.IsWhiteSpace(_s[_i])) _i++;
            }

            public FormatException Error(string what)
            {
                int line = 1;
                for (int k = 0; k < Math.Min(_i, _s.Length); k++) if (_s[k] == '\n') line++;
                return new FormatException(what + " on line " + line);
            }

            public object Value()
            {
                Space();
                if (End) throw Error("unexpected end");
                char c = _s[_i];
                if (c == '{') return Obj();
                if (c == '[') return Arr();
                if (c == '"') return Str();
                if (Match("true")) return true;
                if (Match("false")) return false;
                if (Match("null")) return null;
                if (c == '-' || (c >= '0' && c <= '9')) return Num();
                throw Error("unexpected '" + c + "'");
            }

            private bool Match(string w)
            {
                if (string.CompareOrdinal(_s, _i, w, 0, w.Length) != 0) return false;
                _i += w.Length;
                return true;
            }

            private Dictionary<string, object> Obj()
            {
                if (++_depth > 64) throw Error("nested too deeply");
                var o = new Dictionary<string, object>();
                _i++;
                Space();
                if (!End && _s[_i] == '}') { _i++; _depth--; return o; }
                while (true)
                {
                    Space();
                    if (End || _s[_i] != '"') throw Error("expected a property name");
                    string k = Str();
                    Space();
                    if (End || _s[_i] != ':') throw Error("expected ':'");
                    _i++;
                    o[k] = Value();
                    Space();
                    if (End) throw Error("unexpected end");
                    if (_s[_i] == ',') { _i++; continue; }
                    if (_s[_i] == '}') { _i++; break; }
                    throw Error("expected ',' or '}'");
                }
                _depth--;
                return o;
            }

            private List<object> Arr()
            {
                if (++_depth > 64) throw Error("nested too deeply");
                var a = new List<object>();
                _i++;
                Space();
                if (!End && _s[_i] == ']') { _i++; _depth--; return a; }
                while (true)
                {
                    a.Add(Value());
                    Space();
                    if (End) throw Error("unexpected end");
                    if (_s[_i] == ',') { _i++; continue; }
                    if (_s[_i] == ']') { _i++; break; }
                    throw Error("expected ',' or ']'");
                }
                _depth--;
                return a;
            }

            private string Str()
            {
                var sb = new StringBuilder();
                _i++;
                while (true)
                {
                    if (End) throw Error("unterminated string");
                    char c = _s[_i++];
                    if (c == '"') break;
                    if (c != '\\') { sb.Append(c); continue; }
                    if (End) throw Error("unterminated string");
                    char e = _s[_i++];
                    switch (e)
                    {
                        case '"': sb.Append('"'); break;
                        case '\\': sb.Append('\\'); break;
                        case '/': sb.Append('/'); break;
                        case 'b': sb.Append('\b'); break;
                        case 'f': sb.Append('\f'); break;
                        case 'n': sb.Append('\n'); break;
                        case 'r': sb.Append('\r'); break;
                        case 't': sb.Append('\t'); break;
                        case 'u':
                            int code;
                            if (_i + 4 > _s.Length || !int.TryParse(_s.Substring(_i, 4), NumberStyles.HexNumber, CultureInfo.InvariantCulture, out code))
                                throw Error("bad \\u escape");
                            sb.Append((char)code);
                            _i += 4;
                            break;
                        default: throw Error("bad escape");
                    }
                }
                return sb.ToString();
            }

            private double Num()
            {
                int start = _i;
                if (_s[_i] == '-') _i++;
                while (_i < _s.Length && "0123456789.eE+-".IndexOf(_s[_i]) >= 0) _i++;
                double v;
                if (!double.TryParse(_s.Substring(start, _i - start), NumberStyles.Float, CultureInfo.InvariantCulture, out v))
                    throw Error("bad number");
                return v;
            }
        }
    }
}
