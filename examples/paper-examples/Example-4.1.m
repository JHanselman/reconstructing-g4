/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


  AttachSpec("~/reconstructing-g4/magma/spec");
  /* Consider the hyperelliptic curves
    X1 : y^2 = f1
    X2 : y^2 = f2
    with f1, f2 as follows:
  */

  f1 := 24 * x^5 + 36 * x^4 - 4 * x^3 - 12 * x^2 + 1;
  f2 := 3 * x^5 - 68 * x^4 + 159 * x^3 + 232 * x^2 - 132 * x + 16;
  rootsf2 := [];

  /* Find the maximal isotropic subspace V given by the graph of the isomorphism from Jac(X1)[2] to Jac(X2)[2] that is 
     induced from the given ordering of the Weierstrass points.
  */
  V:= findV(RS1, RS2, roots_f2 );

  //Compute big period matrix
  P1 := BigPeriodMatrix(RS1); P2 := BigPeriodMatrix(RS2);
  P := DiagonalJoin(P1, P2);
  Q := QFromPVFor22(P, V);
 
  //Compute genus 4 curve whose Jacobian is isomorphic to Jac(X1)
  Eqs := RationalReconstructCurveG4(Q);





