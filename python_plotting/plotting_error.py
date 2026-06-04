from mpmath import mp, mpc, mpf

mp.dps = 35  # increase if you want more than C++ max_digits10 (~15–33)

import re
from mpmath import mp, mpc, mpf

complex_pattern = re.compile(r"\(([^,]+),([^)]+)\)")

def parse_scalar(x: str):
    return mpf(x.strip())

def parse_complex(x: str):
    m = complex_pattern.match(x.strip())
    if not m:
        raise ValueError(f"Invalid complex format: {x}")
    return mpc(mpf(m.group(1)), mpf(m.group(2)))


def load_dat(filename):
    data = {
        "l_point": [],
        "l_ref": [],
        "l_tt": []
    }

    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if not line:
            i += 1
            continue

        parts = line.split()
        key = parts[0]

        if key == "l_point":
            n = int(parts[1])
            i += 1
            vals = lines[i].split()
            data["l_point"] = [parse_scalar(v) for v in vals[:n]]

        elif key == "l_ref":
            n = int(parts[1])
            i += 1
            vals = lines[i].split()
            data["l_ref"] = [parse_complex(v) for v in vals[:n]]

        elif key == "l_tt":
            n = int(parts[1])
            i += 1
            vals = lines[i].split()
            data["l_tt"] = [parse_complex(v) for v in vals[:n]]

        i += 1

    return data

import matplotlib.pyplot as plt
from mpmath import mp

def plot_tt_comparison(data, use_mpmath=True):
    """
    data: dict with keys:
        - l_point
        - l_ref (complex mpmath or Python complex)
        - l_tt  (complex mpmath or Python complex)
    """

    x = data["l_point"]
    ref = data["l_ref"]
    tt  = data["l_tt"]

    # ---------- convert to float for plotting ----------
    def cplx_to_float(z):
        return float(z.real), float(z.imag)

    ref_re = []
    ref_im = []
    tt_re = []
    tt_im = []

    diff_re = []
    diff_im = []
    diff_abs = []


    for r, t in zip(ref, tt):
        rr, ri = cplx_to_float(r)
        tr, ti = cplx_to_float(t)

        dr = abs(rr - tr)
        di = abs(ri - ti)

        diff_norm = (dr**2 + di**2) ** 0.5

        ref_re.append(rr)
        ref_im.append(ri)
        tt_re.append(tr)
        tt_im.append(ti)

        diff_re.append(dr)
        diff_im.append(di)
        diff_abs.append(diff_norm)


    # ---------- plotting ----------
    fig, axes = plt.subplots(2, 1, figsize=(12, 14), sharex=True)

    # ===== 1. reference vs TT =====
    axes[0].plot(x, ref_re, label="Ref real")
    axes[0].plot(x, tt_re, "--", label="TT real")
    axes[0].plot(x, ref_im, label="Ref imag")
    axes[0].plot(x, tt_im, "--", label="TT imag")
    axes[0].set_title("Reference vs TT")
    axes[0].legend()
    axes[0].grid(True)

    # ===== 2. absolute error =====
    axes[1].plot(x, diff_re, label="Real error", marker='x')
    axes[1].plot(x, diff_im, label="Imag error", marker='x')
    axes[1].plot(x, diff_abs, label="Abs error", marker='x')
    axes[1].set_title("Absolute Error (ref - tt)")
    axes[1].set_yscale("log")
    axes[1].legend()
    axes[1].grid(True)


    plt.tight_layout()
    plt.show()

a = load_dat("error.dat")

plot_tt_comparison(a)

print(a)