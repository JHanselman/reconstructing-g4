/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


  AttachSpec("~/github/CHIMP/CHIMP.spec");
  AttachSpec("~/github/reconstructing-g4/magma/spec");

  SetDebugOnError(true);
  SetVerbose("User1",1);
  /* Consider the hyperelliptic curves
    X1 : y^2 = f1
    X2 : y^2 = f2
    with f1, f2 as follows:
  */

  prec := 200;
  R<x> := PolynomialRing(Rationals());
  //R<x> := PolynomialRing(RationalsExtra(prec));
  f1 := 24 * x^5 + 36 * x^4 - 4 * x^3 - 12 * x^2 + 1;
  f2 := 3 * x^5 - 68 * x^4 + 159 * x^3 + 232 * x^2 - 132 * x + 16;
  printf "Creating Riemann surfaces with precision %o\n", prec;
  RS1 := RiemannSurface(f1,2 : Precision := prec);
  RS2 := RiemannSurface(f2,2 : Precision := prec);
  K2<nu> := SplittingField(f2);
  roots_f2 := [el[1] : el in Roots(f2,K2)];
  roots_f2 := [roots_f2[el] : el in [1, 3, 2, 5, 4]];

  /* Find the maximal isotropic subspace V given by the graph of the isomorphism from Jac(X1)[2] to Jac(X2)[2] that is 
     induced from the given ordering of the Weierstrass points.
  */
  print "Computing maximal isotropic subspace";
  V:= findV(RS1, RS2, roots_f2 );

  //Compute big period matrix
  print "Computing period matrices";
  P1 := BigPeriodMatrix(RS1); P2 := BigPeriodMatrix(RS2);
  P := DiagonalJoin(P1, P2);
  Q := QFromPVFor22(P, V);
 
  //Compute genus 4 curve whose Jacobian is isomorphic to Jac(X1)
  print "Recovering equations of genus 4 curve";
  Eqs := RationalReconstructCurveG4(Q);





