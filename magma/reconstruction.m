/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */




function RealPart(A)
        return Matrix(Nrows(A), Ncols(A), [[Real(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end function;

function NormalForm(A)
	N := Nrows(A);
	CC := BaseRing(A);
	CCN := VectorSpace(CC, N);
	ReA := RealPart(A);
	ImA := RealPart(-CC.1*A);
        X := BlockMatrix([[ReA, -ImA],[-ImA, -ReA]]);
	D, T := NumericalSchurForm(X);
 	pos := [i: i in [1..2*N] | D[i,i] gt 0];
	if #pos ne N then
		error("Error when computing normal form");
	end if;
        vecs := [[CCN!Eltseq(T[i])[1..N], CCN!Eltseq(T[i])[N+1..2*N]]  : i in pos];
	ret := Matrix([vecs[i][1]+CC.1*vecs[i][2]: i in [1..N ]]  );
	return ret;
end function;

function ComputeThetas(tau)
  CC := BaseRing(tau);
  g := Nrows(tau);
  thetas :=[];
  chars_even := EvenThetaCharacteristics(g);
  for i in [1..(2^(2*g) -1)] do
    characteristic := IndexToTChar(i, g);
    if characteristic in chars_even then
        thetas[i] := Theta([CC!0 : i in [1..g]], tau : char := characteristic);  
    else
        Append(~thetas, CC!0);
    end if;
  end for;
  return thetas;
end function;

function ComputeTritangents(thetas)

  tritangentsys :=
   [[[GF(2)| 0, 1, 1, 0, 0, 1, 0, 0 ], [ GF(2)|0, 1, 1, 0, 1, 1, 0, 0 ]],
    [[GF(2)| 0, 1, 0, 0, 0, 1, 0, 0 ], [GF(2)| 0, 1, 0, 0, 1, 1, 0, 0 ]],
    [[GF(2)| 0, 1, 0, 1, 0, 1, 0, 0 ], [GF(2)| 0, 1, 0, 1, 1, 1, 0, 0 ]],
    [[GF(2)| 0, 1, 1, 1, 0, 1, 0, 0 ], [GF(2)| 0, 1, 1, 1, 1, 1, 0, 0 ]],
    [[GF(2)| 0, 1, 0, 1, 0, 1, 1, 0 ], [GF(2)| 0, 1, 0, 1, 1, 1, 1, 0 ]],
    [[GF(2)| 0, 1, 0, 0, 0, 1, 1, 0 ], [GF(2)| 0, 1, 0, 0, 1, 1, 1, 0 ]],
    [[GF(2)| 0, 1, 0, 0, 0, 1, 1, 1 ], [GF(2)| 0, 1, 0, 0, 1, 1, 1, 1 ]],
    [[GF(2)| 0, 1, 1, 1, 0, 1, 1, 1 ], [GF(2)| 0, 1, 1, 1, 1, 1, 1, 1 ]],
    [[GF(2)| 0, 1, 0, 0, 0, 1, 0, 1 ], [GF(2)| 0, 1, 0, 0, 1, 1, 0, 1 ]],
    [[GF(2)| 0, 1, 1, 0, 0, 1, 0, 1 ], [GF(2)| 0, 1, 1, 0, 1, 1, 0, 1 ]]];

  tritangentbasis := [
    [GF(2)|1, 1, 1, 0, 1, 1, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 1, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 1, 1, 0],
    [GF(2)|0, 1, 1, 0, 0, 1, 0, 0]];
  
  CC := Parent(thetas[1]);
  tritangents := [[[CC!0,0,0,0], [CC!0,0,0,0]]: t in [1..10]];
  constant := [CC!0,0,0,0];
  for k in [1..4] do
    temp_list := [ Vector(t) : t in (tritangentbasis[1..k - 1] cat [tritangentbasis[5]] cat tritangentbasis[k+1 .. 4])];
    S1, S2, A := special_fundamental_system(temp_list);
    signs := signs_in_derivative_formula(A);
    T1 := CC!1;
    T2 := CC!1;
    for s in S1 do
      ss := Eltseq(s);
      T1 *:= thetas[TCharToIndex([ss[1..4],ss[5..8]])];
    end for;

    for s in S2 do
      ss := Eltseq(s);
      T2 *:= thetas[TCharToIndex([ss[1..4],ss[5..8]])];
    end for;
     constant[k] := signs[1]*T1 + signs[2]*T2;
  end for;

  for i in [1..10] do
    for j in [1..2] do
      for k in [1..4] do
        temp_list := [ Vector(t) : t in (tritangentbasis[1..k - 1] cat [tritangentsys[i][j]] cat tritangentbasis[k+1 .. 4])];
        S1, S2, A := special_fundamental_system(temp_list);
        signs := signs_in_derivative_formula(A);
        T1 := CC!1;
        T2 := CC!1;
        for s in S1 do
          ss := Eltseq(s);
          T1 *:= thetas[TCharToIndex([ss[1..4],ss[5..8]])];
        end for;

        for s in S2 do
          ss := Eltseq(s);
          T2 *:= thetas[TCharToIndex([ss[1..4],ss[5..8]])];
        end for;

        tritangents[i][j][k] := (signs[1]*T1 + signs[2]*T2)/constant[k];
      end for;
    end for;
  end for;
  return tritangents;
end function;

function ComputeQuadric(bitangents, tritangents)
  r:= #tritangents;
  CC := BaseRing(tritangents[1][1][1]);
  RR := RealField(Precision(CC));
  CC4:=PolynomialRing(CC,4);
  x:=Matrix(4,1,[CC4.i: i in [1..4]]);
  mats1new:=[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2])): i in [1..r]];
  mats1new:=[(m +Transpose(m))/2 : m in mats1new];
  mats1newx:=Matrix(CC4, 1, r,[(Transpose(x)*ChangeRing(mats1new[i], CC4)*x)[1,1]: i in [1..r]]);

  Xnew:=Matrix([&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats1new]  );
  vi:=NumericalKernel(Xnew: Epsilon:=RR!10^(-15));
  
  CC3 := PolynomialRing(CC, 3);
  
  fs := [&+[el[i]*CC3.i : i in [1..3]] : el in bitangents[1..r]];
  mons:=MonomialsOfDegree(CC3,2);

  fsq_mat := [];
  for f in fs do
    cs := [];
    for m in mons do
      Append(~cs, MonomialCoefficient(f^2,m));
    end for;
    Append(~fsq_mat, cs);
  end for;
  fsq_mat := Matrix(fsq_mat);
  si := NumericalKernel(fsq_mat);
  sirows:=Nrows(si);
  
  N:=HorizontalJoin([  DiagonalMatrix(Eltseq(vi[i]))*fsq_mat : i in [1..Nrows(vi) ]]);
  //TODO: Check singular values to see if rank is too small. If so then compute more tritangents.


  gammaiinv:=NumericalKernel(N: Epsilon:=RR!(10^(-15)));
  D, U, V:=SingularValueDecomposition(Xnew);
  Upart:=Matrix(U[1..7]);
  phi:=Upart*DiagonalMatrix(Eltseq(gammaiinv))*fsq_mat;

  //Kernel is not deterministic
  Qpre:=Kernel(phi);
  Qpre1:=Qpre*Upart;
  Qnew:=&+[Eltseq(Basis(Qpre1)[1])[i]*mats1new[i]: i in [1..r] ];

  return Qnew;
end function;

intrinsic ReconstructCurveG4(tau::AlgMatElt) -> Crv
{}
  g := Nrows(tau);
  thetas := ComputeThetas(tau);
  print thetas;
  tritangents := ComputeTritangents(thetas);
  print tritangents;
end intrinsic;

