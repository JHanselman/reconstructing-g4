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
        thetas[i] := CC!0;
    end if;
  end for;
  thetas[2^(2*g)] := Theta([CC!0 : i in [1..g]], tau : char := [[0: i in [1..g]],[0: i in [1..g]]]);  
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

function ModuliFromTheta(thetas);
  I:=Parent(thetas[1]).1;
  a1:=I*thetas[33]*thetas[5]/(thetas[40]*thetas[12]);
  a2:=I*thetas[21]*thetas[49]/(thetas[28]*thetas[56]);
  a3:=I*thetas[7]*thetas[35]/(thetas[14]*thetas[42]);
  ap1:=I*thetas[5]*thetas[54]/(thetas[27]*thetas[40]);
  ap2:=I*thetas[49]*thetas[2]/(thetas[47]*thetas[28]);
  ap3:=I*thetas[35]*thetas[16]/(thetas[61]*thetas[14]);
  as1:=-thetas[54]*thetas[33]/(thetas[12]*thetas[27]);
  as2:=thetas[2]*thetas[21]/(thetas[56]*thetas[47]);
  as3:=thetas[16]*thetas[7]/(thetas[42]*thetas[61]);
  return [a1,a2,a3,ap1,ap2,ap3,as1,as2,as3];
end function;

function RiemannModelFromModuli(mods);
  a1:=mods[1];a2:=mods[2];a3:=mods[3];
  ap1:=mods[4];ap2:=mods[5];ap3:=mods[6];
  as1:=mods[7];as2:=mods[8];as3:=mods[9];
  F:=Parent(a1);
  P<x1,x2,x3>:=PolynomialRing(F,3);
  k:=1;kp:=1;ks:=1;
  M:=Matrix([[1,1,1],[k*a1,k*a2,k*a3],[kp*ap1,kp*ap2,kp*ap3]]);
  Mb:=Matrix([[1,1,1],[1/a1,1/a2,1/a3],[1/ap1,1/ap2,1/ap3]]);
  U:=-Mb^(-1)*M;
  u1:=U[1];
  u2:=U[2];
  u3:=U[3];
  u1:=u1[1]*x1+u1[2]*x2+u1[3]*x3;
  u2:=u2[1]*x1+u2[2]*x2+u2[3]*x3;
  u3:=u3[1]*x1+u3[2]*x2+u3[3]*x3;
  return (x1*u1+x2*u2-x3*u3)^2-4*x1*u1*x2*u2, u1, u2, u3;
end function;
/* Testing
CC := Parent(thetas[1]);
chars_even := EvenThetaCharacteristics(3);
g3thetas := [];
for i := 1 to 64 do
  s := Intseq(i mod 64,2,6);
  s := Reverse(s);   delta := [s[1..3], s[4..6]];
  if delta in chars_even then
    print delta;
    delta_new1 := [[0] cat el : el in delta];
    delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
    // formula from Lemma 1 (p. 148) of Farkas
    g3thetas[i] := Sqrt(thetas[TCharToIndex(delta_new1)]*thetas[TCharToIndex(\
delta_new2)]);
  else
   Append(~g3thetas, 0);
  end if;
end for;
g3thetas:=correct_signs(g3thetas);
mods := ModuliFromTheta(g3thetas);
quartic := RiemannModelFromModuli(mods);

R:= Parent(quartic);
W<t> := PolynomialRing(CC);

bit := ComputeBitangents(thetas);
for j in [1..28] do
  print j;
  L := &+[bit[j][i] * R.i : i in [1..3]];
  res := Resultant(quartic, L, 1);
  rr := Evaluate(res, [0,t,1]);
  print Roots(rr), "\n";
end for;

function TritangentPlane(Pi, char)
  CC := BaseRing(Pi);
  prec := Precision(CC);
  Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
  tau := Pi1^-1*Pi2;
  g := Nrows(tau);
  cs := [];
  for i := 1 to g do
     dz := [0: i in [1..g]];
    dz[i] := 1;
    Append(~cs, Theta([CC | 0: i in [1..g]], tau : char := char, dz := [dz], prec := prec));
  end for;
  cs := Eltseq(Matrix(1,g,cs)*(Pi1^-1));
  //cs := [cs[i]/cs[g] : i in [1..g]];
  return cs;
end function;
*/
function ComputeBitangents(thetas)

  CC := Parent(thetas[1]);
  chars_even := EvenThetaCharacteristics(3);
  g3thetas := [];
  for i := 1 to 64 do
    s := Intseq(i mod 64,2,6);
    s := Reverse(s);
    delta := [s[1..3], s[4..6]];
    if delta in chars_even then
      print delta;
      delta_new1 := [[0] cat el : el in delta];
      delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
    // formula from Lemma 1 (p. 148) of Farkas
      g3thetas[i] := Sqrt(thetas[TCharToIndex(delta_new1)]*thetas[TCharToIndex(delta_new2)]);
    else
      Append(~g3thetas, 0);
    end if;
  end for;
  
  g3thetas:=correct_signs(g3thetas);
  mods := ModuliFromTheta(g3thetas);

  A := Transpose(Matrix(3,3,[1/el : el in mods]));
  Ainv := A^-1;
  lambdas := Ainv*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);

  mods_mat := [[mods[i], mods[i+1], mods[i+2]] : i in [1,4,7]];
  L := DiagonalMatrix(Eltseq(lambdas));
  B := Transpose(Matrix(3,3, mods))*L;
  Binv := Inverse(B);
  ks := Binv*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
  print ks;
  bitangents := [];
  bitangents := [ [CC | 1, 0, 0], [CC | 0,1,0], [CC | 0,0,1], [CC | 1,1,1]];
  bitangents cat:= mods_mat;
  F, u0, u1, u2 := RiemannModelFromModuli(mods);
  bitangents cat:= [Coefficients(el) : el in [u0, u1, u2]];
  CC3<t0,t1,t2> := Parent(u0);
  bitangents cat:= [Coefficients(el) : el in [t0+t1+u2, t0+u1+t2, u0+t1+t2]];
  mods_mat := Transpose(Matrix(mods_mat));
// (3)
  for i := 1 to 3 do
    new := u0/mods_mat[1,i] + ks[i,1]*(mods_mat[2,i]*t1 + mods_mat[3,i]*t2);
    Append(~bitangents, Coefficients(new));
  end for;
// (4)
  for i := 1 to 3 do
    new := u1/mods_mat[2,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[3,i]*t2);
    Append(~bitangents, Coefficients(new));
  end for;
// (5)
  for i := 1 to 3 do
    new := u2/mods_mat[3,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[2,i]*t1);
    Append(~bitangents, Coefficients(new));
  end for;
// (6)
  modsinv:=Inverse(mods_mat);
  modstra:= Transpose(DiagonalMatrix([1/el: el in  Eltseq(modsinv*Matrix(CC,3,1,[1,1,1]))])*modsinv);
  Atra := Transpose(Matrix(3,3,[1/el : el in Eltseq(modstra)]));
  Ainvtra := Atra^-1;
  lambdastra := Ainvtra*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);
  Ltra := DiagonalMatrix(Eltseq(lambdastra));
  Btra := Transpose(modstra)*Ltra;
  Binvtra := Inverse(Btra);
  kstra := Binvtra*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
  k:=kstra[1,1];
  kp:=kstra[2,1];

  M:=Matrix([[1,1,1],Eltseq(k*modstra[1]),Eltseq(kp*modstra[2])]);
  Mb:=Matrix([[1,1,1],Eltseq(Transpose(Atra)[1]),Eltseq(Transpose(Atra)[2])]);
  U:=-Mb^(-1)*M;
  u1tra:=U[1]*Inverse(modstra);
  u2tra:=U[2]*Inverse(modstra);
  u3tra:=U[3]*Inverse(modstra);

  bitangents cat:= [Eltseq(el) : el in [u1tra, u2tra, u3tra]];
// (7)
  for i := 1 to 3 do
    new := u0/(mods_mat[1,i]*(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + u1/(mods_mat[2,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + u2/(mods_mat[3,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
    Append(~bitangents, Coefficients(new));
  end for;

  return bitangents;

end function;

function ComputeSquareRootOnP1xP1(detqdualonsegre)

  CC := BaseRing(Parent(detqdualonsegre));
  CC4:=PolynomialRing(CC,4);
  P1P1<x1,x2,y1,y2> := PolynomialRing(CC,4);
  x:=[x1,x2];
  y:=[y1,y2];
  f := Evaluate(detqdualonsegre, [x1*y1, x2*y2, x1*y2, x2*y1]);

  //xy[a,b] = x1^a * x2^(3-a) * y1^b * y2^(3-b) (with 0<= a,b ,+3) 
  xy:= AssociativeArray();
  P3mons:= AssociativeArray();

  //P3mons defines pushforward of P1 x P1 -> P3 where P3mons[[a,b]] is the image of xy[a,b]
  P3mons[[3,3]] := CC4.1^3;
  P3mons[[3,0]] := CC4.3^3;
  P3mons[[0,3]]:= CC4.4^3;
  P3mons[[0,0]] := CC4.2^3;
  P3mons[[3,1]] := CC4.3^2*CC4.1;
  P3mons[[3,2]] := CC4.1^2*CC4.3;
  P3mons[[0,1]] := CC4.2^2*CC4.4;
  P3mons[[0,2]] := CC4.4^2*CC4.2;
  P3mons[[1,3]] := CC4.4^2*CC4.1;
  P3mons[[2,3]] := CC4.1^2*CC4.4;
  P3mons[[1,0]] := CC4.2^2*CC4.3;
  P3mons[[2,0]] := CC4.3^2*CC4.2;
  P3mons[[2,1]] :=CC4.3^2 * CC4.4;
  P3mons[[2,2]] := CC4.1^2 * CC4.2;
  P3mons[[1,1]] := CC4.2^2 * CC4.1;
  P3mons[[1,2]] := CC4.4^2 * CC4.3;


  //Compute the coefficients of the monomials in the sqaureroot of f.
  xy[[3,3]] := Sqrt(MonomialCoefficient(f, x1^6*y1^6));
  xy[[3,0]] := Sqrt(MonomialCoefficient(f, x1^6*y2^6));
  xy[[0,3]] := Sqrt(MonomialCoefficient(f, x2^6*y1^6));
  xy[[0,0]] := Sqrt(MonomialCoefficient(f, x2^6*y2^6));

  max:= 0;
  start := [];
  for i in Keys(xy) do
    absval := Abs(xy[i]);
    if absval gt max then
      max := absval;
      start := i;
    end if;
  end for;

  hstep := 1;
  vstep := 1;

  if start[1] gt 0 then
    hstep := -1;
  end if;

  if start[2] gt 0 then
    vstep := -1;
  end if;

  for i in [0..3] do
    for j in [0..3] do
      if i eq 0 and j eq 0 then
        continue;
      end if;
      n1 := start[1] + hstep * i; n2 := 3 - n1;
      n3 := start[2] + vstep * j; n4 := 3 - n3;
      mon := MonomialCoefficient(f, x1^(start[1] +n1)*x2^(3 - start[1] + n2) * y1^(start[2] + n3)*y2^(3 - start[2] + n4));
      rect := &cat[[[[start[1] + hstep * mu,start[2] + vstep * nu], [start[1] + hstep * (i - mu),start[2] + vstep * (j - nu)]]  : mu in [0..i] | (mu ne 0 or nu ne 0) and (mu ne i or nu ne j)] : nu in [0..j]];
      print "i =", i, "j=", j, "coords ", n1, n3, "\n";
      print rect;
      print x1^(start[1] +n1)*x2^(3 - start[1] + n2) * y1^(start[2] + n3)*y2^(3 - start[2] + n4), "\n";
    
      subtractsum := &+([xy[tup[1]] * xy[tup[2]] : tup in rect] cat [CC!0]);
      xy[[n1,n3]] := ((mon - subtractsum  )/xy[start])/2;
    end for;
  end for;

  sqrt:= P1P1!0;
  SegreCubic := CC4!0;

  for i in [0..3] do
    for j in [0..3] do
      sqrt +:= xy[[i,j]] *x1^i*x2^(3-i)*y1^j*y2^(3-j);
      SegreCubic +:= xy[[i,j]] * P3mons[[i,j]];
    end for;
  end for;
  return SegreCubic;
end function;

function ComputeCurve(bitangents, tritangents)
  r:= #tritangents;
  CC := Parent(tritangents[1][1][1]);
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
  dualelt:=mats1newx*ChangeRing(Transpose(Upart), CC4);
  
  D, U, V:=SingularValueDecomposition(phi);
  phiext:=HorizontalJoin(phi, Matrix(7,1, Eltseq(Conjugate(U)[7])));
  phiTinv:=ChangeRing(Transpose(phiext)^(-1), CC4);
  phiL:=dualelt*phiTinv;
  qdual:=ZeroMatrix(CC4, 3,3);
  count:=1;
  for i in [1..3] do
    for j in [i..3] do
      qdual[i,j]:=phiL[1,count];
	qdual[j, i]:=phiL[1, count];
	count+:=1;	
    end for;
  end for;
  detqdual:=Determinant(qdual);
  
  S:= NormalForm(Qnew);

  //Compute the matrix S2 whose inverse will transform the normal form into one where all coefficients are 1.

  S2 := S * Qnew * Transpose(S);
  for i in [1..4] do
    S2[i,i] := Sqrt(S2[i,i]);
  end for;

  I:=CC.1;
  //Coordinate Transformation that maps x^2 + y^2 +z^2 +w^2 to xy - zw.
  IdToSegre := Matrix(CC, 4, 4, [[0,1/2,0,1/2*I],[0,1/2,0,-1/2*I],[1/2,0,-1/2*I,0], [-1/2, 0,-1/2*I,0]]);

  //Complete Transformation
  QtoSegre := IdToSegre * (S2)^(-1) * S; 


  v:=Matrix(CC4, [[CC4.1, CC4.2, CC4.3, CC4.4]]);

  //Apply coordinate transformation to detqdual to map the quadric to the Segre quadric
  detqdualonsegre := Evaluate(detqdual, Eltseq((v * ChangeRing(QtoSegre, CC4))[1]));
  
  //Pull back to P1 x P1 and take the square root there
  SegreCubic := ComputeSquareRootOnP1xP1(detqdualonsegre);

  //Reverse the coordinate transformation
  cubic := Evaluate(SegreCubic, Eltseq((v * ChangeRing(QtoSegre^(-1), CC4))[1]));
  quadric:=(v*ChangeRing(quadric, CC4) *Transpose(v))[1,1];
  return quadric, cubic;


end function;

intrinsic ReconstructCurveG4(tau::AlgMatElt)->SeqEnum
{}
  g := Nrows(tau);
  thetas := ComputeThetas(tau);
  tritangents := ComputeTritangents(thetas);
  bitangents := ComputeBitangents(thetas);
  quadric, cubic := ComputeCurve([bitangents[i]: i in [10, 23, 4, 20,  17, 9, 12,  1, 5, 11]], tritangents);
  return [quadric, cubic];
end intrinsic;

