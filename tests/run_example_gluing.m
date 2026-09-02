/***
 *  Reconstructing genus 4 curves
 *
 *  End-to-end check of the FLINT theta path: run examples/paper-examples/
 *  Example-Gluing.m (glue two genus 2 Jacobians along their 2-torsion and
 *  reconstruct the genus 4 curve, which goes through ComputeThetas ->
 *  ThetaFlint at 200 digits and TritangentPlanes with ord := 1) and check the
 *  result coordinate-independently: Jac(C) is isogenous to Jac(X1) x Jac(X2),
 *  so #C(F_p) = #X1(F_p) + #X2(F_p) - (p+1) for every prime p of good
 *  reduction.  (The equations recorded in QMexample.m for the same pair of
 *  curves are a different model, so they are not compared directly.)
 *
 *  Run from examples/paper-examples (the example uses relative paths):
 *      cd examples/paper-examples && magma ../../tests/run_example_gluing.m
 */

t0 := Cputime();
load "Example-Gluing.m";
printf "\nExample-Gluing.m finished in %o s\n", Cputime(t0);
print "Recovered equations:";
print Eqs;

QQ := Rationals();
R<x,y,z,w> := PolynomialRing(QQ, 4);
got := [ R!Evaluate(e, [x,y,z,w]) : e in Eqs ];
P3 := ProjectiveSpace(R);
C := Curve(Scheme(P3, got));
printf "\nrecovered curve: genus %o (expected 4), irreducible %o\n", Genus(C), IsIrreducible(C);

// clear denominators so that we can reduce mod p
den := LCM([ Denominator(c) : c in &cat[ Coefficients(e) : e in got ] ]);
gotZ := [ den*e : e in got ];
// Jac(C) is Qbar-isogenous to Jac(X1) x Jac(X2); over Q the isogeny may only
// be defined over a quadratic field, so a_p(C) = chi(p) (a_p(X1) + a_p(X2))
// for a quadratic character chi.  (There is no "correct twist" of C to look
// for: a non-hyperelliptic curve has no quadratic twists, and by Torelli/Serre
// exactly one of the two Q-forms A, A^(-1) of the glued variety is a Jacobian
// -- C is the curve whose Jacobian it is.)  We check this and identify chi.
ok := true; tested := []; chi := AssociativeArray(); bad := {2, 3};
for p in PrimesInInterval(5, 200) do
  Fp := GF(p);
  // good reduction of X1, X2: f_i squarefree mod p (deg 6, so also fine at infinity)
  f1p := PolynomialRing(Fp)!f1; f2p := PolynomialRing(Fp)!f2;
  if den mod p eq 0 or not IsSquarefree(f1p) or not IsSquarefree(f2p) then Include(~bad, p); continue; end if;
  P3p := ProjectiveSpace(Fp, 3); Rp := CoordinateRing(P3p);
  Sp := Scheme(P3p, [ Rp!e : e in gotZ ]);
  if Dimension(Sp) ne 1 or IsSingular(Sp) then Include(~bad, p); continue; end if;   // bad reduction of C
  Cp := Curve(Sp);
  aC := p + 1 - #RationalPoints(Cp);
  aX := (p + 1 - #Points(HyperellipticCurve(f1p))) + (p + 1 - #Points(HyperellipticCurve(f2p)));
  good := aC eq aX or aC eq -aX;
  ok and:= good;
  Append(~tested, p);
  if aX ne 0 then chi[p] := (aC eq aX) select 1 else -1; end if;
  if p le 100 or not good then
    printf "  p = %3o: a_p(C) = %4o, a_p(X1) + a_p(X2) = %4o   %o\n", p, aC, aX, good select (aC eq aX or aX eq 0 select "ok" else "ok (chi = -1)") else "MISMATCH";
  end if;
end for;
printf "\nJac(C) ~ Jac(X1) x Jac(X2) up to the sign character, checked at %o primes: %o\n", #tested, ok;
if not ok or #tested lt 5 then
  error "run_example_gluing: reconstruction does not match the glued Jacobian";
end if;
printf "chi = -1 at: %o\n", Sort([ p : p in Keys(chi) | chi[p] eq -1 ]);
// chi is unramified outside the bad primes: find its discriminant
S := Sort(Setseq(bad));
cands := [];
for sign in [1, -1] do
  for sub in Subsets({1..#S}) do
    d := sign * &*[Integers() | S[i] : i in sub];
    if d eq 1 then continue; end if;
    if forall{ p : p in Keys(chi) | KroneckerSymbol(d, p) eq chi[p] } then Append(~cands, d); end if;
  end for;
end for;
printf "quadratic character(s) supported on %o matching the traces: %o\n", S, cands;
if #cands eq 1 then
  printf "=> the isogeny Jac(C) -> Jac(X1) x Jac(X2) is defined over Q(sqrt(%o))\n", cands[1];
end if;

// Optional: compare with the model recorded in QMexample.m for the same two
// curves using Bouchet's Genus-4 package (IsIsomorphicGenus4: covariants /
// invariants + Groebner search for the PGL_4 matrix), if it can be attached.
// Genus-4 needs Thomas Bouchet's GL-Equivalence package (TransvectantGeneral
// etc.), so attach that first.  Both are looked for next to this repo and in ~.
attached := false;
function TryAttach(specs)
  for spec in specs do
    try
      AttachSpec(spec); return true;
    catch e
      _ := true;
    end try;
  end for;
  return false;
end function;
if TryAttach(["~/github/GL-Equivalence/magma/spec", "~/GL-Equivalence/magma/spec", "../../../GL-Equivalence/magma/spec"]) then
  attached := TryAttach(["~/github/Genus-4/magma/spec", "~/Genus-4/magma/spec", "../../../Genus-4/magma/spec"]);
end if;
if attached then
  expected := [
    x^2 + 28/5*x*y + 18/5*y^2 - 224/5*z^2 + 102/5*z*w - 33/10*w^2,
    -27/226*x*y^2 - 205/1808*x*z^2 + 15/113*x*z*w - 185/3616*x*w^2 - 39/452*y^3 + y*z^2 - 271/452*y*z*w + 17/452*y*w^2
  ];
  print "\nComparing with the model recorded in QMexample.m via Bouchet's invariants...";
  invs, wgt := InvariantsGenus4Curves(got[1], got[2]);
  invs1 := InvariantsGenus4Curves(expected[1], expected[2]);
  same_invs := WPSEqual(wgt, invs, invs1);
  printf "equal invariants (same point in weighted projective space): %o\n", same_invs;
  // A generic genus 4 curve has no nontrivial twists, so equal invariants means
  // isomorphic over Q.  IsIsomorphicGenus4 would also produce the matrix, but
  // its Groebner branch currently errors on this pair, so only try it.
  try
    t1 := Cputime();
    iso := IsIsomorphicGenus4(got[1], got[2], expected[1], expected[2]);
    printf "IsIsomorphicGenus4: %o (%o s)\n", iso, Cputime(t1);
    if iso then _, M := IsIsomorphicGenus4(got[1], got[2], expected[1], expected[2]); print M; end if;
  catch e
    printf "IsIsomorphicGenus4 raised an error inside the Genus-4 package (%o); relying on the invariants above\n", e`Object;
  end try;
  if not same_invs then
    error "run_example_gluing: recovered curve differs from the recorded model";
  end if;
else
  print "\n(Genus-4 package not found; skipping IsIsomorphicGenus4 comparison with the recorded model)";
end if;
quit;
