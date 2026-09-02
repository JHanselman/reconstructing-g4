/***
 *  Reconstructing genus 4 curves
 *
 *  Test: ThetaFlint (Cext/acb_theta_cli, pure FLINT) against
 *    (a) CHIMP's python-flint route (Theta.magma/PythonFlint.m), run through
 *        the same python + script CHIMP uses, if that venv (or a python3 with
 *        python-flint) is available;
 *    (b) the pure-Magma Theta() in magma/theta.m at low precision;
 *    (c) parity: odd characteristics vanish at z = 0, even ones do not
 *        (this pins down the characteristic numbering);
 *    (d) FLINT's derivatives against central finite differences, and
 *        ThetaDerivativesNumerical / TritangentPlanes consistency.
 *
 *  Run from the repository root:
 *      magma tests/theta_flint_vs_chimp.m
 *  Optionally:  CHIMP=/path/to/CHIMP magma tests/theta_flint_vs_chimp.m
 *  (default: ../CHIMP next to this repository).
 */

AttachSpec("magma/spec");
SetVerbose("ThetaFlint", 1);

QQ := Rationals();
failures := 0;
procedure check(~failures, cond, msg)
  if cond then
    printf "  ok   %o\n", msg;
  else
    printf "  FAIL %o\n", msg;
    failures +:= 1;
  end if;
end procedure;

// ---------------------------------------------------------------------------
// test data: a genus 2 and a genus 4 small period matrix from hyperelliptic
// curves, via Magma's Riemann surface machinery
// ---------------------------------------------------------------------------
function SmallPeriodMatrixOf(f, prec)
  S := RiemannSurface(f, 2 : Precision := prec);
  return SmallPeriodMatrix(S);
end function;

prec := 60;
CC<i> := ComplexField(prec);
R<x> := PolynomialRing(QQ);
tau2 := SmallPeriodMatrixOf(x^5 - x^4 + 3*x^3 - 2*x + 1, prec);
// In genus 4 a hyperelliptic Jacobian has 10 vanishing even theta constants,
// which would confuse test (c), so use a generic point of the Siegel upper
// half space instead: tau4 = X + i*(A*A^T + 1) with X symmetric.
X4 := Matrix(QQ, 4, 4, [ 3/10, -1/5, 1/8, 2/7,  -1/5, 1/6, -2/9, 1/3,  1/8, -2/9, -1/4, 1/10,  2/7, 1/3, 1/10, 2/5 ]);
A4 := Matrix(QQ, 4, 4, [ 1/2, -1/3, 1/5, 2/3,  1/4, 3/4, -1/2, 1/6,  -2/5, 1/7, 1/3, -1/2,  1/3, 1/2, 1/4, -1/3 ]);
tau4 := ChangeRing(X4, CC) + i*ChangeRing(A4*Transpose(A4) + 1, CC);

// ---------------------------------------------------------------------------
// (a) CHIMP's python-flint values
// ---------------------------------------------------------------------------
chimp_root := GetEnv("CHIMP");
if chimp_root eq "" then chimp_root := "../CHIMP"; end if;
python := "";
// CHIMP creates its venv lazily (first call of its ThetaFlint), so it may not
// exist yet; to build it without Magma:
//   python3 -m venv --without-pip CHIMP/Theta.magma/flint && python3 -m pip install --pre \
//     python-flint==0.7.0a4 --target CHIMP/Theta.magma/flint/lib/python3.X/site-packages
// or simply `pip3 install python-flint` so that plain python3 works.
for cand in [chimp_root cat "/Theta.magma/flint/bin/python", "python3", "python"] do
  try
    _ := Pipe(Sprintf("%o -c 'import flint' 2>/dev/null", cand), "");
    python := cand; break;
  catch e
    _ := true;
  end try;
end for;

// verbatim from CHIMP/Theta.magma/PythonFlint.m, except that ReplaceCharacter
// (a CHIMP intrinsic) is replaced by a local equivalent
ReplaceChar := func<s, c, d | &cat[ ch eq c select d else ch : ch in Eltseq(s) ]>;
to_arb := func<elt | ReplaceChar(Sprintf("arb('%o +/- %.*o')", elt, 3, Max(Abs(elt)*10^(1-Precision(Parent(elt))), 10^(1-Precision(Parent(elt)))) ), "E", "e")>;
to_acb := func<elt | Sprintf("acb(%o, %o)", to_arb(Real(elt)), to_arb(Imaginary(elt)))>;
to_acb_list := func<elt | Sprintf("[%o]", Join([to_acb(x) : x in elt], ", "))> ;
to_acb_matrix := func<elt | Sprintf("acb_mat([%o])", Join([to_acb_list(Eltseq(x)) : x in Rows(elt)], ", ")) >;
chimp_script := "
from flint import acb_mat, acb, arb, ctx

log_10_2 = arb.const_log2()/arb.const_log10()
def arb_to_magma(x):
    digits = (max(x.rel_accuracy_bits(), x.rel_one_accuracy_bits())*log_10_2).floor().unique_fmpz()
    digits_str = f'p{digits}'
    return x.str(digits, more=True, radius=False) + digits_str, digits
def acb_to_magma(x):
    real_str, real_digits = arb_to_magma(x.real)
    imag_str, imag_digits = arb_to_magma(x.imag)
    return f'[{real_str}, {imag_str}]', min(real_digits, imag_digits)
def acb_entries_to_magma(m):
    pairs = [acb_to_magma(x) for x in m]
    digits = min(x[1] for x in pairs)
    return f'[ ComplexField({digits}) | ' + ', '.join([x[0] for x in pairs])  + ']'

ctx.dps = %o
tau = %o
z = %o
theta = tau.theta(z)
print(acb_entries_to_magma(theta))
";
function ChimpThetaFlint(z, tau)
  digits := Precision(BaseRing(tau));
  g := Nrows(tau);
  acb_z := IsZero(z) select Sprintf("acb_mat([[0] for _ in range(%o)])", g) else to_acb_matrix(z);
  cmd := Sprintf(chimp_script, Ceiling(digits*1.1 + 10), to_acb_matrix(tau), acb_z);
  return eval Pipe(python, cmd);
end function;

print "(a) against CHIMP's python-flint route";
if python eq "" then
  printf "  skipped: no python with python-flint found (looked for %o/Theta.magma/flint/bin/python and python3;\n           set CHIMP=/path/to/CHIMP, or `pip3 install python-flint`)\n", chimp_root;
else
  printf "  using %o\n", python;
  for tau in [* tau2, tau4 *] do
    g := Nrows(tau); CC := BaseRing(tau);
    for z in [ZeroMatrix(CC, g, 1), Matrix(CC, g, 1, [CC | (k/10) + (k/7)*CC.1 : k in [1..g]])] do
      ours := ThetaFlint(z, tau);
      theirs := ChimpThetaFlint(z, tau);
      check(~failures, #ours eq 4^g and #theirs eq 4^g, Sprintf("g = %o: %o values", g, 4^g));
      D := Min(Precision(Universe(ours[1])), Precision(Universe(theirs)));
      err := Max([Abs(ours[i][1] - theirs[i]) : i in [1..4^g]]);
      check(~failures, err lt 10^(5-D), Sprintf("g = %o, z %o: max difference %o at %o digits", g, IsZero(z) select "= 0" else "generic", RealField(3)!err, D));
    end for;
  end for;
end if;

// ---------------------------------------------------------------------------
// (b) against the pure-Magma theta series in magma/theta.m (slow: low precision)
// ---------------------------------------------------------------------------
print "(b) against the pure-Magma Theta() series";
lowprec := 20;
for tau in [* tau2, tau4 *] do
  g := Nrows(tau);
  CClow := ComplexField(lowprec);
  taulow := ChangeRing(tau, CClow);
  ours := ThetaFlint(ZeroMatrix(CClow, g, 1), taulow);
  // a handful of characteristics, including the zero one and some odd ones
  ns := [0, 1, 2^g, 2^g + 1, 4^g - 1, 3*2^g + 5 mod 2^g];
  err := RealField(lowprec)!0;
  try
    for n in ns do
      char := IndexToTChar(n, g);   // [a, b] with entries in {0,1}
      magma_val := Theta([CClow | 0 : i in [1..g]], taulow : char := char);
      err := Max(err, Abs(magma_val - ours[n+1][1]));
    end for;
    check(~failures, err lt 10^(6-lowprec), Sprintf("g = %o: max difference %o over chars %o", g, RealField(3)!err, ns));
  catch e
    // theta.m's Theta() needs Real/Imaginary on matrices, which this repo
    // currently gets from CHIMP; skip rather than fail when CHIMP is absent
    printf "  skipped (Theta() raised: %o)\n", e`Object;
  end try;
end for;

// ---------------------------------------------------------------------------
// (c) parity: theta_ab(0, tau) = 0 iff (a,b) odd (i.e. a.b odd)
// ---------------------------------------------------------------------------
print "(c) parity of theta constants (pins down the characteristic numbering)";
for tau in [* tau2, tau4 *] do
  g := Nrows(tau); CC := BaseRing(tau);
  th := ThetaConstantsFlint(tau);
  D := Precision(Universe(th));
  odd := []; even := [];
  for n in [0..4^g-1] do
    ab := IndexToTChar(n, g);
    parity := (&+[ab[1][i]*ab[2][i] : i in [1..g]]) mod 2;
    if parity eq 1 then Append(~odd, n); else Append(~even, n); end if;
  end for;
  check(~failures, #odd eq 2^(g-1)*(2^g-1), Sprintf("g = %o: %o odd characteristics", g, #odd));
  check(~failures, forall{n : n in odd | Abs(th[n+1]) lt 10^(5-D)}, Sprintf("g = %o: all odd theta constants vanish", g));
  check(~failures, forall{n : n in even | Abs(th[n+1]) gt 10^(-10)}, Sprintf("g = %o: all even theta constants are nonzero", g));
  // and the matrix form of a characteristic agrees with the integer form
  c := ThetaCharacteristicFromInteger(odd[3], g);
  check(~failures, ThetaCharacteristicToInteger(c) eq odd[3] and ThetaFlint(c, ZeroMatrix(CC, g, 1), tau)[1] eq th[odd[3]+1],
        Sprintf("g = %o: ThetaFlint(char, z, tau) and integer numbering agree", g));
  // ComputeThetas moves theta_0 to the end
  ct := ComputeThetas(tau);
  check(~failures, #ct eq 4^g and ct[4^g] eq th[1] and ct[1] eq th[2], Sprintf("g = %o: ComputeThetas(tau) layout", g));
end for;

// ---------------------------------------------------------------------------
// (d) derivatives
// ---------------------------------------------------------------------------
print "(d) derivatives against finite differences";
for tau in [* tau2, tau4 *] do
  g := Nrows(tau); CC := BaseRing(tau);
  z0 := Matrix(CC, g, 1, [CC | (k/10) + (k/7)*CC.1 : k in [1..g]]);
  jets := ThetaFlint(z0, tau : ord := 2);
  check(~failures, #jets[1] eq 1 + g + g*(g+1) div 2, Sprintf("g = %o: %o jet entries per characteristic at ord = 2", g, #jets[1]));
  h := CC!10^(-Ceiling(prec/3));   // 1e-20 at prec 60: FD error ~ h^2 + eps/h ~ 1e-37
  err := RealField(prec)!0;
  for i in [1..g] do
    zp := z0; zp[i,1] +:= h;
    zm := z0; zm[i,1] -:= h;
    fp := ThetaFlint(zp, tau); fm := ThetaFlint(zm, tau);
    for n in [0..4^g-1] do
      fd := (fp[n+1][1] - fm[n+1][1]) / (2*h);
      err := Max(err, Abs(fd - jets[n+1][1+i]));
    end for;
  end for;
  check(~failures, err lt 10^(-Ceiling(prec/3) - 6), Sprintf("g = %o: first derivatives, max difference %o", g, RealField(3)!err));
  // second derivative d^2/dz_1^2 is the (g+2)-nd entry (tuple order: (2,0,..) first among order 2)
  h2 := CC!10^(-Ceiling(prec/4));  // 1e-15: FD error ~ h2^2 + eps/h2^2 ~ 1e-27
  zp := z0; zp[1,1] +:= h2;
  zm := z0; zm[1,1] -:= h2;
  fp := ThetaFlint(zp, tau); fm := ThetaFlint(zm, tau);
  err2 := Max([Abs((fp[n+1][1] - 2*jets[n+1][1] + fm[n+1][1])/h2^2 - jets[n+1][g+2]) : n in [0..4^g-1]]);
  // FD error here is ~ theta'''' h2^2/12 + eps/h2^2, i.e. ~1e-25 at prec 60
  check(~failures, err2 lt 10^(-2*Ceiling(prec/4) + 9), Sprintf("g = %o: d^2/dz_1^2 (un-normalised), max difference %o", g, RealField(3)!err2));
end for;

// ThetaDerivativesNumerical and TritangentPlanes use ord := 1 in genus 4
print "(d') ThetaDerivativesNumerical / TritangentPlanes";
g := 4; CC := BaseRing(tau4);
oddchars := OddThetaCharacteristics(4);
grad := ThetaDerivativesNumerical(oddchars[1], tau4);
jets := ThetaFlint(ZeroMatrix(CC, g, 1), tau4 : ord := 1);
n := ThetaCharacteristicToInteger(Matrix(QQ, 8, 1, &cat oddchars[1])/2);
check(~failures, #grad eq 4 and forall{i : i in [1..4] | grad[i] eq jets[n+1][1+i]}, "ThetaDerivativesNumerical returns FLINT's gradient at z = 0");
Pi := HorizontalJoin(IdentityMatrix(CC, 4), tau4);
planes := TritangentPlanes(Pi, [Matrix(QQ, 8, 1, &cat c)/2 : c in oddchars[1..3]]);
check(~failures, #planes eq 3 and forall{i : i in [1..4] | Abs(planes[1][i] - grad[i]) lt 10^(10-prec)}, "TritangentPlanes(Pi = [1 | tau]) agrees with ThetaDerivativesNumerical");

printf "\n%o failure(s)\n", failures;
if failures gt 0 then
  error "theta_flint_vs_chimp: some checks failed";
end if;
quit;
