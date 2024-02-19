/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */

intrinsic SchottkyModularForm(tau::AlgMatElt : prec := -1) -> Any
  {}
   
  if prec eq -1 then
    prec := Precision(BaseRing(Parent(tau)));
  end if;
  C := ComplexField(prec);
  char := Matrix(C, 8, 1, [0,0,0,0,0,0,0,0]);
  z := Matrix(C, 4, 1, [0,0,0,0]);

  m1 := 1/2*Matrix(C, 8, 1, [1,0,1,0,1,0,1,0]);
  m2 := 1/2*Matrix(C, 8, 1, [0,0,0,1,1,0,0,0]);
  m3 := 1/2*Matrix(C, 8, 1, [0,0,1,1,1,0,1,1]);
  n0 := Matrix(C, 8, 1, [0,0,0,0,0,0,0,0]);
  n1 := 1/2*Matrix(C, 8, 1, [0,0,0,1,1,1,1,0]);
  n2 := 1/2*Matrix(C, 8, 1, [0,0,1,1,0,0,0,1]);
  n3 := 1/2*Matrix(C, 8, 1, [0,0,1,0,1,0,1,1]);
  n4 := n1+n2;
  n5 := n1+n3;
  n6 := n2+n3;
  n7 := n1+n2+n3;
  SchottkyN := [n0,n1,n2,n3,n4,n5,n6,n7];
  M1 := [m1 + n: n in SchottkyN];
  M2 := [m2 + n: n in SchottkyN];
  M3 := [m3 + n: n in SchottkyN];
  pi1 := 1;
  pi2 := 1;
  pi3 := 1;

  /*
  function CharacteristicMatrixToPair(c)
    ZZ := Integers();
    QQ := Rationals();
    c *:= 2;
    c := [QQ!(ZZ!(GF(2)!(ZZ!el))) : el in Eltseq(c)];
    return [c[1..4], c[5..8]];
  end function;

  M1 := [CharacteristicMatrixToPair(el) : el in M1];
  M2 := [CharacteristicMatrixToPair(el) : el in M2];
  M3 := [CharacteristicMatrixToPair(el) : el in M3];

  z := Eltseq(z);
  */

  tau_prec := MatrixAlgebra(C, Nrows(tau))!tau;

  for m in M1 do
    pi1 := pi1*Theta(m, z, tau_prec);
    //pi1 := pi1*Theta(z, tau : char := m, prec := prec);
  end for;

  for m in M2 do
    pi2 := pi2*Theta(m, z, tau_prec);
    //pi2 := pi2*Theta(z, tau : char := m, prec := prec);
  end for;

  for m in M3 do
    pi3 := pi3*Theta(m, z, tau_prec);
    //pi3 := pi3*Theta(z, tau : char := m, prec := prec);
  end for;

  Schottky := pi1^2 + pi2^2 + pi3^2 - 2*(pi1*pi2 + pi2*pi3 + pi1*pi3);
  return Schottky;
end intrinsic;
