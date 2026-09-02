/***
 *  Reconstructing genus 4 curves
 *
 *  End-to-end check of the FLINT theta path on examples/paper-examples/
 *  Example-Modular.m: reconstruct the genus 4 curve from the (hard-coded,
 *  300-digit) big period matrix of the modular abelian fourfold attached to
 *  the newform 778.2.a.a, then verify the result against the LMFDB Hecke
 *  eigenvalues in 778.2.a.a.m:  L_p(C, T) = prod_sigma (1 - sigma(a_p) T + p T^2)
 *  for the four embeddings sigma of the Hecke field, at good primes p.
 *
 *  Run from examples/paper-examples:
 *      cd examples/paper-examples && magma ../../tests/run_example_modular.m
 *  (Example-Modular.m itself hard-codes ~/reconstructing-g4 and ~/Genus-4 in
 *  its AttachSpec lines, so the period matrix is read out of the file rather
 *  than loading it.)
 */

AttachSpec("../../magma/spec");
SetVerbose("Reconstruction", 1);

// --- the period matrix, read from Example-Modular.m --------------------------
CC<I> := ComplexField(300);
src := Read("Example-Modular.m");
i0 := Index(src, "Pi :=");
i1 := i0 + Index(src[i0..#src], "];");
Pi := eval src[i0 + 5 .. i1 - 1];
assert Nrows(Pi) eq 4 and Ncols(Pi) eq 8;

// --- reconstruction over Q, exactly as in the example -------------------------
t0 := Cputime();
Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
Eqs := RationalReconstructCurveG4(HorizontalJoin(Pi2, Pi1));
printf "\nRationalReconstructCurveG4 took %o s\n", Cputime(t0);
print "Recovered equations:";
print Eqs;

// --- the second half of the example: reconstruction over CC at 75 digits -----
CC1 := ComplexField(150);
Pi_ := ChangeRing(Pi, CC1);
Pi1_, Pi2_ := SplitBigPeriodMatrix(Pi_);
tau := -Pi1_^(-1) * Pi2_;
tau := (tau + Transpose(tau))/2;
tau_red, Q := SiegelReduction(tau);
t0 := Cputime();
Eqstest := ReconstructCurveG4(ChangeRing(tau_red, ComplexField(75)) : method := "Cayley");
printf "ReconstructCurveG4 over CC (75 digits, Cayley) took %o s; %o equations\n", Cputime(t0), #Eqstest;

// --- verification against the newform 778.2.a.a ----------------------------------
QQ := Rationals();
R<x,y,z,w> := PolynomialRing(QQ, 4);
got := [ R!Evaluate(e, [x,y,z,w]) : e in Eqs ];
C := Curve(ProjectiveSpace(R), got);
printf "\nrecovered curve: genus %o (expected 4), irreducible %o\n", Genus(C), IsIrreducible(C);

load "778.2.a.a.m";
an := qexpCoeffs();              // a_n, n >= 1, in the Hecke field Kf = Q(nu)
Kf := Universe(an);
printf "Hecke field: %o\n", DefiningPolynomial(Kf);
den := LCM([ Denominator(c) : c in &cat[ Coefficients(e) : e in got ] ]);
gotZ := [ den*e : e in got ];

ok := true; tested := [];
for p in PrimesInInterval(3, 60) do
  if 778 mod p eq 0 or den mod p eq 0 then continue; end if;
  Fp := GF(p);
  P3p := ProjectiveSpace(Fp, 3); Rp := CoordinateRing(P3p);
  Sp := Scheme(P3p, [ Rp!e : e in gotZ ]);
  if Dimension(Sp) ne 1 or IsSingular(Sp) then printf "  p = %2o: bad reduction of the model, skipped\n", p; continue; end if;
  t0 := Cputime();
  Lp := LPolynomial(Curve(Sp));
  ap := an[p];
  // L_p(C,T) should be the norm from Kf[T] to Q[T] of 1 - a_p T + p T^2
  hecke := Polynomial(Kf, [1, -ap, p]);
  divisible := IsDivisibleBy(ChangeRing(Lp, Kf), hecke);
  trace_ok := Coefficient(Lp, 1) eq -Trace(ap) and Degree(Lp) eq 8;
  good := divisible and trace_ok;
  ok and:= good;
  Append(~tested, p);
  printf "  p = %2o: a_p(C) = %3o, Tr(a_p(f)) = %3o, (1 - a_p T + p T^2) | L_p: %o   %o  (%o s)\n",
    p, -Coefficient(Lp, 1), Trace(ap), divisible, good select "ok" else "MISMATCH", Cputime(t0);
end for;
printf "\nL-polynomials match the newform 778.2.a.a at %o primes: %o\n", #tested, ok;
if not ok or #tested lt 5 then
  error "run_example_modular: reconstructed curve does not match the newform";
end if;
quit;
