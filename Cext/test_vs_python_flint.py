#!/usr/bin/env python3
"""Compare ./acb_theta_cli against python-flint's acb_mat.theta -- the exact
call CHIMP's Theta.magma/PythonFlint.m makes -- and against finite differences
for the first derivatives.  Run from this directory after `make`.

    pip install python-flint numpy mpmath
    python3 test_vs_python_flint.py
"""
import re
import subprocess
import sys
from decimal import Decimal

import mpmath
import numpy as np
from flint import acb, acb_mat, ctx

BIN = "./acb_theta_cli"


def rand_tau(g, rng):
    X = rng.uniform(-0.5, 0.5, (g, g))
    X = (X + X.T) / 2
    A = rng.uniform(-1, 1, (g, g))
    Y = A @ A.T + np.eye(g)
    return [[complex(X[i][j], Y[i][j]) for j in range(g)] for i in range(g)]


def run_cli(g, digits, ord_, z, tau):
    bits = int(digits * 3.3219280948873626) + 64
    lines = [f"{g} {bits} {ord_}"]
    for i in range(g):
        lines += [str(Decimal(z[i].real)), str(Decimal(z[i].imag))]
    for i in range(g):
        for j in range(g):
            lines += [str(Decimal(tau[i][j].real)), str(Decimal(tau[i][j].imag))]
    open("/tmp/acb_theta_in.txt", "w").write("\n".join(lines) + "\n")
    r = subprocess.run([BIN, "/tmp/acb_theta_in.txt", "/tmp/acb_theta_out.txt"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"{BIN} failed: {r.stderr}")
    out = open("/tmp/acb_theta_out.txt").read()
    D = int(re.search(r"ComplexField\((\d+)\)", out).group(1))
    rows = []
    for m in re.finditer(r"\[ ComplexField\(\d+\) \|(.*?)\]\s*,?\s*\n", out):
        pairs = re.findall(r"\[(-?[0-9][^,\]]*),\s*(-?[0-9][^\]]*)\]", m.group(1))
        rows.append([(mpmath.mpf(a.split("p")[0]), mpmath.mpf(b.split("p")[0])) for a, b in pairs])
    return D, rows


def chimp(g, digits, z, tau):
    """What CHIMP's call_python_flint does (exact double inputs here)."""
    ctx.dps = int(digits * 1.1 + 10)
    T = acb_mat([[acb(tau[i][j].real, tau[i][j].imag) for j in range(g)] for i in range(g)])
    Z = acb_mat([[acb(zz.real, zz.imag)] for zz in z])
    return list(T.theta(Z))


def to_mpf(x, digits):
    return mpmath.mpf(x.mid().str(digits + 10, radius=False))


def main():
    mpmath.mp.dps = 130
    rng = np.random.default_rng(7)
    ok = True
    for g in (2, 4):
        tau = rand_tau(g, rng)
        for z in ([0j] * g, [complex(rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)) for _ in range(g)]):
            digits = 100
            D, rows = run_cli(g, digits, 1, z, tau)
            ref = chimp(g, digits, z, tau)
            assert len(rows) == 4 ** g and all(len(r) == g + 1 for r in rows)
            err = max(max(abs(rows[k][0][0] - to_mpf(ref[k].real, digits)),
                          abs(rows[k][0][1] - to_mpf(ref[k].imag, digits))) for k in range(4 ** g))
            # first derivatives vs central differences of python-flint values
            h = mpmath.mpf(10) ** -30
            derr = mpmath.mpf(0)
            for i in range(g):
                # python-flint takes acb inputs, so we can use a tiny exact step
                ctx.dps = int(digits * 1.1 + 10)
                T = acb_mat([[acb(tau[a][b].real, tau[a][b].imag) for b in range(g)] for a in range(g)])
                zp = [acb(zz.real, zz.imag) for zz in z]
                zm = list(zp)
                step = acb(1) / acb(10) ** 30
                zp[i] = zp[i] + step
                zm[i] = zm[i] - step
                fp = list(T.theta(acb_mat([[w] for w in zp])))
                fm = list(T.theta(acb_mat([[w] for w in zm])))
                for k in range(4 ** g):
                    fd = (fp[k] - fm[k]) / (2 * step)
                    derr = max(derr, abs(rows[k][1 + i][0] - to_mpf(fd.real, digits)),
                               abs(rows[k][1 + i][1] - to_mpf(fd.imag, digits)))
            zdesc = "0" if not any(z) else "generic"
            good = err < mpmath.mpf(10) ** -(digits - 5) and derr < mpmath.mpf(10) ** -35
            ok &= good
            print(f"g={g} z={zdesc}: D={D} certified digits, "
                  f"max|theta - python-flint| = {mpmath.nstr(err, 3)}, "
                  f"max|dtheta - finite diff| = {mpmath.nstr(derr, 3)}  {'ok' if good else 'FAIL'}")
    print("PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
