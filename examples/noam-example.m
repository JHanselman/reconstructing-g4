SetDefaultRealField(RealField(300));

  N := 778;
  M := ModularSymbols(N,2);
  H1 := CuspidalSubspace(M);  // apparently this represents H_1(X_0(N),Q)
  H1N := NewSubspace(H1);      // This is not in the H142E23 example
  time D := NewformDecomposition(H1N);
  d := [Dimension(S) : S in D];
  d;
/*
[ 8, 14, 18, 24 ] -- remember these are twice the dimension of the corresponding J_0(N) factors
*/
  S := D[1];
  time P4 := Periods(S,40000);
  P := Matrix(P4);
  P;
  time E := EndomorphismRing(Transpose(P));
  E.1;
  Factorization(CharacteristicPolynomial(E.1));

// RationalReconstructCurveG4 or ReconstructCurveG4
