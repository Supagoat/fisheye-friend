"""
Golden reference data for the defish model, computed independently of the C# code.

For each case (frame shape and adjustment) and each point (u, v) of the defished frame (normalised, 0..1 across its
width and height), the point (s, t) of the fisheye photo it shows (normalised the same way), in double precision.

The model is the one in defish-hemi-model.json:
  lens:        angle from the axis  theta(r) = sum theta[i] r^(2i+1), r = distance from the photo's centre in photo widths
  projection:  x = k (P+1) sin(lon) / (P + cos lon) (1 + sum h[i] sin^(2i+2) lon)
               y = k (Q+1) / (Q + cos lon) tan(lat) (1 + sum v[i][j] sin^(2i+2) lat sin^(2j) lon)
  adjustment:  strength s divides the lens' focal length (r_s(theta) = r(theta) / s) and k is rescaled so the photo's
               edge midpoints stay at the frame's edges; the amounts a blend the displacement with none,
               p + a (model(p) - p), each multiplied by its shape curve at the distance from the centre (0 = centre,
               1 = corner): a monotone cubic through four points (centre, 1/3, 2/3, corner), flat at the centre,
               harmonic-mean tangents inside (0 where the slope changes sign), one-sided at the corner, and straight on
               past the corners (outside the frame).

The equations are solved by bisection to the last bit (not by Newton's method as in the C# code), so this is a
separate implementation to test against.

Usage: python3 golden.py ../../Assets/StreamingAssets/defish-hemi-model.json defish-golden.json
"""
import json
import math
import sys


def load(path):
    m = json.load(open(path))
    return {
        "theta": m["lens"]["theta"],
        "k": m["projection"]["scale"],
        "P": m["projection"]["panini"],
        "Q": m["projection"]["paniniVertical"],
        "h": m["projection"].get("horizontal", []),
        "v": m["projection"].get("vertical", []),
    }


def lens_angle(m, r):
    return sum(c * r ** (2 * i + 1) for i, c in enumerate(m["theta"]))


def lens_radius(m, th):
    # theta(r) increases over the range used; bracket and bisect
    lo, hi = 0.0, 1.0
    while lens_angle(m, hi) < th:
        hi *= 2
    return bisect(lambda r: lens_angle(m, r) - th, lo, hi)


def proj_x(m, lon):
    s = math.sin(lon)
    u = s * s
    return (m["P"] + 1) * s / (m["P"] + math.cos(lon)) * (1 + sum(c * u ** (i + 1) for i, c in enumerate(m["h"])))


def proj_y(m, lon, lat):
    u = math.sin(lon) ** 2
    w = math.sin(lat) ** 2
    corr = 1 + sum(m["v"][i][j] * w ** (i + 1) * u ** j for i in range(len(m["v"])) for j in range(len(m["v"][i])))
    return (m["Q"] + 1) / (m["Q"] + math.cos(lon)) * math.tan(lat) * corr


def bisect(f, lo, hi):
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        fm = f(mid)
        if fm == 0:
            return mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def shape_at(y, rho):
    n = len(y)
    h = 1.0 / (n - 1)

    def tangent(k):
        if k == 0:
            return 0.0
        if k == n - 1:
            return (y[k] - y[k - 1]) / h
        d0 = (y[k] - y[k - 1]) / h
        d1 = (y[k + 1] - y[k]) / h
        if d0 * d1 <= 0:
            return 0.0
        return 2 / (1 / d0 + 1 / d1)

    if rho <= 0:
        return y[0]
    if rho >= 1:
        return y[-1] + tangent(n - 1) * (rho - 1)     # past the corners (outside the frame): straight on
    k = min(n - 2, int(rho / h))
    t = (rho - k * h) / h
    m0, m1 = tangent(k), tangent(k + 1)
    t2, t3 = t * t, t * t * t
    return (2 * t3 - 3 * t2 + 1) * y[k] + (t3 - 2 * t2 + t) * h * m0 + (-2 * t3 + 3 * t2) * y[k + 1] + (t3 - t2) * h * m1


def map_point(m, u, v, aspect, strength=1.0, ax=1.0, ay=1.0, shape_x=(1, 1, 1, 1), shape_y=(1, 1, 1, 1)):
    x = u - 0.5
    y = (v - 0.5) * aspect
    k = m["k"] * proj_x(m, lens_angle(m, 0.5)) / proj_x(m, lens_angle(m, 0.5 * strength))
    lon = bisect(lambda a: k * proj_x(m, a) - x, -1.55, 1.55)
    lat = bisect(lambda a: k * proj_y(m, lon, a) - y, -1.5, 1.5)
    X = math.cos(lat) * math.sin(lon)
    Y = math.sin(lat)
    Z = math.cos(lat) * math.cos(lon)
    sin_th = math.hypot(X, Y)
    th = math.atan2(sin_th, Z)
    if th < 1e-12:
        sx, sy = X / (m["theta"][0] * strength), Y / (m["theta"][0] * strength)
    else:
        r = lens_radius(m, th) / strength
        sx, sy = r * X / sin_th, r * Y / sin_th
    rho = math.hypot(x, y) / math.sqrt(0.25 + 0.25 * aspect * aspect)
    a = ax * shape_at(shape_x, rho)
    b = ay * shape_at(shape_y, rho)
    sx = x + a * (sx - x)
    sy = y + b * (sy - y)
    return sx + 0.5, sy / aspect + 0.5


CASES = [
    # name, aspect (height / width), strength, amount x, amount y, shape x, shape y
    ("default 3:2", 2 / 3, 1.0, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("strength 50%", 2 / 3, 0.5, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("strength 115%", 2 / 3, 1.15, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("strength 80%, amounts 60%", 2 / 3, 0.8, 0.6, 0.6, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("no correction", 2 / 3, 1.0, 0.0, 0.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("horizontal 150%, vertical 50%", 2 / 3, 1.0, 1.5, 0.5, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("horizontal 30%, vertical 120%", 2 / 3, 1.0, 0.3, 1.2, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("curve softening the edges", 2 / 3, 1.0, 1.0, 1.0, [1, 1, 0.6, 0.3], [1, 1, 0.6, 0.3]),
    ("separate curves", 2 / 3, 0.9, 1.1, 0.9, [0.6, 0.8, 1, 1.3], [1, 0.9, 0.7, 0.5]),
    ("curve with a turn", 2 / 3, 1.05, 0.8, 0.8, [1, 1.4, 0.9, 1.2], [0.5, 1, 1, 0.8]),
    ("4:3 frame", 3 / 4, 1.0, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("16:9 frame", 9 / 16, 1.0, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("square frame", 1.0, 1.0, 1.0, 1.0, [1, 1, 1, 1], [1, 1, 1, 1]),
    ("4:3, everything adjusted", 3 / 4, 1.1, 0.8, 1.3, [1, 1.2, 0.8, 0.6], [0.9, 1, 1.1, 1.2]),
]


def points():
    pts = [(i / 12, j / 8) for j in range(9) for i in range(13)]
    # a few in between, from a fixed sequence (no random module: the same everywhere)
    a, b = 0.5, 0.5
    for _ in range(23):
        a = (a + 0.6180339887498949) % 1.0
        b = (b + 0.7548776662466927) % 1.0
        pts.append((a, b))
    return pts


def main():
    m = load(sys.argv[1])
    out = {"format": "fisheye-friend-defish-golden", "model": "defish-hemi-model.json",
           "note": "s, t = normalised photo point shown at normalised result point u, v; double precision (golden.py)",
           "cases": []}
    for name, aspect, st, ax, ay, shx, shy in CASES:
        rows = []
        for u, v in points():
            s, t = map_point(m, u, v, aspect, st, ax, ay, shx, shy)
            rows.append([u, v, s, t])
        out["cases"].append({"name": name, "aspect": aspect, "strength": st, "amountX": ax, "amountY": ay,
                             "shapeX": shx, "shapeY": shy, "points": rows})
    with open(sys.argv[2], "w") as f:
        json.dump(out, f, indent=1)
    print("%d cases, %d points each" % (len(CASES), len(points())))


if __name__ == "__main__":
    main()
