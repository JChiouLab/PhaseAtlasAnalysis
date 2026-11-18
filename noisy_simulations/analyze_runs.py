"""Classify every Fig-4 run (one npz per (A, rep)) and write a summary table.

Criterion (documented in figures/fig4/README.md):
  polarity state  s(t) = +1 if dG > +H, -1 if dG < -H, else keep the previous state
                  (hysteresis H = 0.10 on the atlas coordinate dG = n_G,L - n_G,R);
  reversal        a change of s between +1 and -1;
  dwell           time between consecutive reversals (the first dwell starts at the
                  first time |dG| > H, the last one ends at T);
  epoch label     a dwell longer than T_LONG = 3 x the deterministic half-period
                  (median A = 0 settled dwell) is a monopolar epoch M, otherwise an
                  oscillation epoch O; consecutive equal labels are merged.
  oscillating                : epochs == O
  damped_to_monopolar        : epochs == M, or O then M with M starting before T_SETTLE
                               (the initial transient rings down into one pole, stays)
  oscillating_then_monopolar : O then M, M starting after T_SETTLE (one late switch)
  switching                  : any other sequence (>= 2 switches, or M then O)
  other_never_polarized      : |dG| never exceeds H
  failed                     : the solver refused (fixed-point non-convergence)
Usage: analyze_runs.py <runs_dir> <out_prefix> [<A0_runs_dir>]
"""
import ast
import glob
import json
import os
import sys

import numpy as np

H = 0.10
T_SETTLE = 1000.0


def load(fn):
    z = np.load(fn, allow_pickle=False)
    meta = ast.literal_eval(str(z["meta"]))
    d = {k: z[k] for k in z.files if k != "meta"}
    return meta, d


def events(t, dG):
    s = 0
    st = np.zeros_like(dG)
    rev = []
    first = None
    for i, g in enumerate(dG):
        ns = 1 if g > H else (-1 if g < -H else s)
        if ns != 0 and first is None:
            first = t[i]
        if s != 0 and ns != s:
            rev.append(t[i])
        s = ns
        st[i] = s
    return first, np.array(rev), st


def dwells(first, rev, T):
    if first is None:
        return []
    b = np.r_[first, rev, T]
    return [(b[k], b[k + 1]) for k in range(len(b) - 1)]


def classify(t, dG, t_long):
    T = t[-1]
    first, rev, st = events(t, dG)
    dw = dwells(first, rev, T)
    if first is None:
        return "other_never_polarized", dict(n_rev=0, n_rev_settled=0, n_long=0, long_dwells=[],
                                             epochs="", n_switches=0, first_polar=None, last_rev=None)
    lab = ["M" if (b - a) > t_long else "O" for (a, b) in dw]
    seq = []                      # collapsed epoch sequence with start times
    for (a, b), l in zip(dw, lab):
        if not seq or seq[-1][0] != l:
            seq.append([l, a, b])
        else:
            seq[-1][2] = b
    s = "".join(e[0] for e in seq)
    if s == "O":
        c = "oscillating"
    elif s == "M" or (s == "OM" and seq[1][1] < T_SETTLE):
        c = "damped_to_monopolar"
    elif s == "OM":
        c = "oscillating_then_monopolar"
    else:
        c = "switching"
    rev_s = rev[rev >= T_SETTLE]
    long_dw = [(a, b) for (a, b), l in zip(dw, lab) if l == "M"]
    return c, dict(n_rev=int(len(rev)), n_rev_settled=int(len(rev_s)), n_long=len(long_dw),
                   long_dwells=[(float(a), float(b)) for a, b in long_dw],
                   epochs="-".join("%s[%.0f to %.0f]" % (e[0], e[1], e[2]) for e in seq),
                   n_switches=len(seq) - 1,
                   first_polar=float(first), last_rev=float(rev[-1]) if len(rev) else None)


def main():
    runs, prefix = sys.argv[1], sys.argv[2]
    a0dir = sys.argv[3] if len(sys.argv) > 3 else runs
    # deterministic half-period from the A = 0 runs (settled dwells)
    hp = []
    for fn in sorted(glob.glob(os.path.join(a0dir, "A0_rep*.npz"))):
        meta, d = load(fn)
        if meta["status"] != "ok":
            continue
        first, rev, _ = events(d["t"], d["dG"])
        r = rev[rev >= T_SETTLE]
        hp += list(np.diff(r))
    half = float(np.median(hp)) if hp else 275.0
    t_long = 3.0 * half
    rows = []
    for fn in sorted(glob.glob(os.path.join(runs, "A*_rep*.npz"))):
        meta, d = load(fn)
        row = dict(file=os.path.basename(fn), A=meta["A"], rep=meta["rep"], mode=meta["mode"],
                   T=meta["T"], dt=meta["dt"], wall_s=round(meta["wall_s"], 1),
                   noise_seed=meta["noise_seed"], ic_seed=meta["ic_seed"], status=meta["status"])
        if meta["status"] != "ok":
            row["class"] = "failed"
        else:
            c, info = classify(d["t"], d["dG"], t_long)
            row["class"] = c
            row.update(info)
            row["M_GEF_range"] = [float(d["M_GEF"].min()), float(d["M_GEF"].max())]
            row["M_PAK_range"] = [float(d["M_PAK"].min()), float(d["M_PAK"].max())]
            row["M_PAK_end"] = float(d["M_PAK"][-1])
            row["min_conc"] = float(d["min_conc"].min())
            row["max_abs_dG_settled"] = float(np.abs(d["dG"][d["t"] >= T_SETTLE]).max())
        rows.append(row)
    rows.sort(key=lambda r: (r["A"], r["rep"]))
    out = dict(H=H, T_SETTLE=T_SETTLE, half_period_A0=half, t_long=t_long, runs=rows)
    with open(prefix + "_classification.json", "w") as f:
        json.dump(out, f, indent=1)
    classes = sorted({r["class"] for r in rows})
    As = sorted({r["A"] for r in rows})
    with open(prefix + "_classification.csv", "w") as f:
        f.write("A,rep,class,n_switches,n_rev,n_rev_settled,n_long,first_polar,last_rev,M_PAK_end,min_conc,wall_s,epochs,status\n")
        for r in rows:
            f.write(",".join(str(r.get(k, "")) for k in
                             ["A", "rep", "class", "n_switches", "n_rev", "n_rev_settled", "n_long", "first_polar",
                              "last_rev", "M_PAK_end", "min_conc", "wall_s"]) + ',"%s","%s"\n' % (r.get("epochs", ""), r["status"]))
    print("half-period(A=0) %.1f  t_long %.1f" % (half, t_long))
    print("A".ljust(8) + "".join(c[:22].ljust(24) for c in classes))
    for A in As:
        print(str(A).ljust(8) + "".join(str(sum(1 for r in rows if r["A"] == A and r["class"] == c)).ljust(24)
                                        for c in classes))


if __name__ == "__main__":
    main()
