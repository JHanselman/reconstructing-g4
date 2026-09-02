# Cext — external C solver for `ThetaFlint`

`acb_theta_cli.c` is a ~300-line front end to FLINT's
[`acb_theta`](https://flintlib.org/doc/acb_theta.html) module.  It evaluates
all `2^(2g)` theta functions `theta_{a,b}(z, tau)` — and, on request, all their
partial derivatives in `z` up to a given total order — in rigorous ball
arithmetic, and prints the result as a Magma expression.  `ThetaFlint` in
`magma/FlintWrapper.m` writes a small input file, runs this binary and `eval`s
the output, exactly as `Belyi/Cext/powser_arnoldi` does for `PowerSeriesBasis`.

This replaces two earlier routes to the same FLINT code: CHIMP's
`Theta.magma/PythonFlint.m` (which builds a venv, pip-installs `python-flint`
and pipes a Python script) and the old `FlintWrapper.m` that shelled out to
FLINT's bundled example program (`~/flint/build/examples/acb_theta`, which
requires a source build of FLINT with examples).  Neither Python nor a FLINT
source tree is needed any more — only `libflint`.

## Build

Needs FLINT >= 3.1 (the `acb_theta` module was added in 3.1; 3.2 or later is
preferred and is what Homebrew ships).

* macOS: `brew install flint mpfr gmp`, then `make mac`.
* Linux with a recent distro FLINT: `make` (needs `libflint-dev` >= 3.1, `libmpfr-dev`,
  `libgmp-dev`).  Note that Ubuntu 24.04's `libflint-dev` is 3.0.1, which is too old.
* Servers without admin access / old system FLINT: `sh build_deps.sh` builds
  GMP, MPFR and FLINT into `$HOME/.local/acb_theta`, then `make server`
  static-links against them.

Sanity check: `make test` runs `./acb_theta_cli --selftest`, which checks
`theta_00(0, i) = pi^(1/4) / Gamma(3/4)` in genus 1 and compares first and
second derivatives in genus 2 against finite differences.

`test_vs_python_flint.py` (optional; needs `pip install python-flint numpy mpmath`)
compares the binary against `acb_mat.theta`, which is what CHIMP's
`PythonFlint.m` calls, to ~100 digits in genus 2 and 4.

## Use from Magma

```
AttachSpec("magma/spec");
// the binary is found automatically at <repo>/Cext/acb_theta_cli; otherwise
//   export THETA_FLINT_BIN=/path/to/acb_theta_cli
tau := ...;                                  // symmetric g x g, Im(tau) > 0
CC := BaseRing(tau); g := Nrows(tau);
th := ThetaFlint(ZeroMatrix(CC, g, 1), tau); // theta constants, th[ab+1] = [theta_ab(0,tau)]
th := ThetaFlint(ZeroMatrix(CC, g, 1), tau : ord := 1);
// th[ab+1] = [theta_ab, d/dz_1 theta_ab, ..., d/dz_g theta_ab]
ThetaFlint(char, z, tau : ord := 1);         // one characteristic, same layout
SetVerbose("ThetaFlint", 1);                 // timings; 2 also shows the command
```

Characteristics are numbered `ab = 2^g * a + b`, bits most significant first
(the FLINT and CHIMP convention); `ThetaCharacteristicToInteger` /
`ThetaCharacteristicFromInteger` convert.  Results are cached per
`(z, tau, ord, precision)`; `CacheClearThetaFlint()` empties the cache.

The certified precision of the output is decided by FLINT: every returned
value lives in `ComplexField(D)` where `D` is the number of decimal digits
FLINT guarantees for all entries, given that the inputs are only known to
Magma's working precision (they are passed as balls of relative radius
`10^(1 - prec)`).  In practice `D` is about `prec - 3`.

## File formats (for debugging)

Input: a header `g prec_bits ord`, then `2g` numbers `Re z_i`, `Im z_i`, then
`2g^2` numbers for the entries of `tau` row by row (`Re`, `Im`).  A number is
either a decimal `MID` or a ball `MID +/- RAD`, in any form `arb_set_str`
accepts; Magma's upper-case `E` exponents are fine.  The file is tokenised on
whitespace after removing backslash-newline continuations, so it does not
matter where Magma's line wrapping (`SetColumns`) puts newlines or `\`.

Output: `[ [ ComplexField(D) | [re, im], ... ], ... ]`, one inner sequence per
characteristic in the order above; inside, the theta value first, then the
`g` first partials in `z`, then higher partials in FLINT's tuple order
(by total order, then reverse-lexicographic).  Unlike FLINT's raw jets these
are true partial derivatives (the `1/k!` normalisation is undone).

Exit status: 0 success; 1 usage or parse error; 2 if fewer than 5 digits could
be certified (usually: `Im(tau)` not positive definite, or precision far too low).
