/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


import "rosenhain.m": FindDelta, ComputeGamma, EtaFunction, TakaseQuotient;

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

intrinsic ComputeThetas(tau::AlgMatElt) -> SeqEnum
{}
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
end intrinsic;

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
     if IsVerbose("User1",1) then
        print "precision loss in constant=", Abs(constant[k]/T1);
  end if;

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
	if IsVerbose("User1",1) then
              print "precision loss in numerator of JDF=", Abs((signs[1]*T1 + signs[2]*T2)/T1);
        end if;

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
intrinsic ComputeBitangents(thetas::SeqEnum) -> SeqEnum
  {}

  CC := Parent(thetas[1]);
  chars_even := EvenThetaCharacteristics(3);
  g3thetas := [];
  for i := 1 to 64 do
    s := Intseq(i mod 64,2,6);
    s := Reverse(s);
    delta := [s[1..3], s[4..6]];
    if delta in chars_even then
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
  mods_mat := [[mods[i], mods[i+1], mods[i+2]] : i in [1,4,7]];
  ks := Matrix(CC, 3, 1, [1,1,1]);
  bitangents := [];
  bitangents := [ [CC | 1, 0, 0], [CC | 0,1,0], [CC | 0,0,1], [CC | 1,1,1]];
  bitangents cat:= mods_mat;
  F, u0, u1, u2 := RiemannModelFromModuli(mods);
  bitangents cat:= [Coefficients(el) : el in [u0, u1, u2]];
  CC3<t0,t1,t2> := Parent(u0);
  bitangents cat:= [Coefficients(el) : el in [t0+t1+u2, t0+u1+t2, u0+t1+t2]];
  if IsVerbose("User1",1) then
     print "precision loss in bitangent comp=", [[Abs(bitangents[8][i]/bitangents[7][i]), Abs(bitangents[9][i]/bitangents[6][i]),Abs(bitangents[10][i]/bitangents[5][i])]: i in [1..3]];
  end if;
  mods_mat := Transpose(Matrix(mods_mat));
// (3)
  for i := 1 to 3 do
    new := u0/mods_mat[1,i] + ks[i,1]*(mods_mat[2,i]*t1 + mods_mat[3,i]*t2);
    Append(~bitangents, Coefficients(new));
    if IsVerbose("User1",1) then
      print "precision loss in bitangent comp=",  [[Abs(mods_mat[1,i]*bitangents[10+i][nu]/bitangents[5][nu])]: nu in [1..3]];
    end if;
  end for;
// (4)
  for i := 1 to 3 do
    new := u1/mods_mat[2,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[3,i]*t2);
    Append(~bitangents, Coefficients(new));
    if IsVerbose("User1",1) then
       print "precision loss in bitangent comp=",  [[Abs(mods_mat[2,i]*bitangents[13+i][nu]/bitangents[6][nu])]: nu in [1..3]];
    end if;

  end for;
// (5)
  for i := 1 to 3 do
    new := u2/mods_mat[3,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[2,i]*t1);
    Append(~bitangents, Coefficients(new));
    if IsVerbose("User1",1) then
      print "precision loss in bitangent comp=", [[Abs(mods_mat[3,i]*bitangents[16+i][nu]/bitangents[7][nu])]: nu in [1..3]];
    end if;

  end for;
// (6)
  modsinv:=Inverse(mods_mat);
  D := DiagonalMatrix([1/el: el in Eltseq(modsinv*Matrix(CC,3,1,[1,1,1]))]);
  modstra:= Transpose(D*modsinv);
  Atra := Transpose(Matrix(3,3,[1/el : el in Eltseq(modstra)]));
  lambdastra := Solution(Transpose(Atra), Vector([BaseRing(Parent(Atra)) | -1,-1,-1]));
  Ltra := DiagonalMatrix(Eltseq(lambdastra));
  Btra := Transpose(modstra)*Ltra;
  kstra := Solution(Transpose(Btra), Vector([BaseRing(Parent(Btra)) | -1,-1,-1]));
  printf "kstra = %o\n", kstra;
  k:=kstra[1];
  kp:=kstra[2];

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
  /*
  if IsVerbose("User1",1) then
      DA:=SingularValueDecomposition(A);
      DB:=SingularValueDecomposition(B);
      Dmods_mat:=SingularValueDecomposition(mods_mat);
      Dmodstra:=SingularValueDecomposition(modstra);
      DAtra:=SingularValueDecomposition(A);
      DBtra:=SingularValueDecomposition(Btra);
      DMb:=SingularValueDecomposition(Mb);
      print "\n SVs of A", Diagonal(DA),"\n SVs of B", Diagonal(DB),"\n SVs of mods_mat", Diagonal(Dmods_mat),"\n SVs of modstra", Diagonal(Dmodstra),"\n SVs of Atra", Diagonal(DAtra),"\n SVs of Btra", Diagonal(DBtra),"\n SVs of Mb", Diagonal(DMb);
  end if;
  */
  return bitangents;

end intrinsic;

function OddThetaPrymInfoHyp(thetas)

  CC := Parent(thetas[1]);
  chars_even := EvenThetaCharacteristics(3);
  g3thetassq := [];
  for i := 1 to 64 do
    s := Intseq(i mod 64,2,6);
    s := Reverse(s);
    delta := [s[1..3], s[4..6]];
    if delta in chars_even then
      delta_new1 := [[0] cat el : el in delta];
      delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
    // formula from Lemma 1 (p. 148) of Farkas
      g3thetassq[i] := thetas[TCharToIndex(delta_new1)]*thetas[TCharToIndex(delta_new2)];
    else
      Append(~g3thetassq, 0);
    end if;
  end for;
  rosens:=RosenhainInvariantsFromThetaSquares(g3thetassq, 3: SpTrafo:=false);
  wp:=[[CC|1,0,0], [CC|1,1,1]] cat [[ro^j: j in [0..2]]: ro in rosens] cat [[CC| 0,0,1]];
  
  bitangents:=[];
  for i in [1..8] do
     for j in [i+1..8] do
        Mat:=Transpose(Matrix([wp[i], wp[j]]));
        Append(~bitangents, Eltseq(NumericalKernel(Mat)));
     end for;
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
  if IsVerbose("User1",1) then
       print "\n Absolute value of leading coeff of sqrt=", absval, "\n";
  end if;
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
  print sqrt^2-f;
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
  DN:=SingularValueDecomposition(N);
  gammaiinv:=NumericalKernel(N: Epsilon:=RR!(10^(-15)));
  
  DXnew, U, V:=SingularValueDecomposition(Xnew);
  Upart:=Matrix(U[1..7]);
  phi:=Upart*DiagonalMatrix(Eltseq(gammaiinv))*fsq_mat;

  //Kernel is not deterministic
  
  Qpre:=Kernel(phi);
  Qpre1:=Qpre*Upart;
  Qnew:=&+[Eltseq(Basis(Qpre1)[1])[i]*mats1new[i]: i in [1..r] ];
  dualelt:=mats1newx*ChangeRing(Transpose(Upart), CC4);
  
  D, U, V:=SingularValueDecomposition(phi);
  phiext:=HorizontalJoin(phi, Matrix(7,1, Eltseq(Conjugate(U)[7])));
  if IsVerbose("User1", 1) then
        Dfsqmat:=SingularValueDecomposition(fsq_mat);
        Dphiext:=SingularValueDecomposition(Transpose(phiext));

        print "SVs of fsq_mat", Diagonal(Dfsqmat),"\n SVs of Xnew=", Diagonal(DXnew), "\n SVs of N=", Diagonal(DN), "\n SVs of phi=", Diagonal(D), "\n SVs of phiext=", Diagonal(Dphiext);
  end if;


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
  quadric:=(v*ChangeRing(Qnew, CC4) *Transpose(v))[1,1];
  return quadric, cubic;

end function;

function liftSpN(A, n)
    ZZ:=Integers();
    p:=Characteristic(BaseRing(A));
    assert(BaseRing(A) eq GF(p));
    N:=Nrows(A);
    assert(IsDivisibleBy(N,2));
    g:=N div 2;
    id:=IdentityMatrix(ZZ, g);
    zer:=ZeroMatrix(ZZ, g,g);
    J:=BlockMatrix([[zer, id],[-id, zer]]);
    JGFp:=ChangeRing(J, GF(p));
    Alift:=ChangeRing(A, ZZ);
    for i in [2..n] do
	Mat:=(J-Transpose(Alift)*J*Alift) div p^(i-1);
	Math:=Matrix(GF(p), [[(mu ge nu select Mat[mu,nu] else 0): nu in [1..N]]: mu in [1..N]]);
        X1:=(Transpose(A)*JGFp)^(-1)*Math;
        Alift+:=p^(i-1)*ChangeRing( X1, ZZ);
    end for;
    return Alift;
end function;






function MapSympl(vec1)
	//Finds a symplectic transformation mapping vec1 to vec2
	vec2:=Vector([GF(2)|0,1,1,1,0,1,0,1]);
	ZZ:=Integers();
        id:=IdentityMatrix(ZZ, 4);
        zer:=ZeroMatrix(ZZ, 4,4);
        idGF2:=IdentityMatrix(GF(2), 4);
        zerGF2:=ZeroMatrix(GF(2), 4,4);
	J1:=BlockMatrix([[zerGF2, idGF2],[zerGF2, zerGF2]]);
        J:=BlockMatrix([[zerGF2, idGF2],[idGF2, zerGF2]]);

        if vec1 eq 0 then
		//Maybe transpose/inverse of this matrix?

		return Matrix(ZZ, 8,8, [[1, 0, 0, 0, 0, 0, 0, 0]
                                       [0, 1, 1, 1, 0, -1, 1, 1],
                                       [0, 0, 1, 0, 0, 0, -1, 0],
                                       [0, 0, 1, 0, 0, -1, 1, 1],
                                       [0, 0, 0, 0, 1, 0, 0, 0],
                                       [0, -1, -1, -1, 0, 2, -1, -1],
                                       [0, 0, 0, -1, 0, 1, -1, -1],
                                       [0, 0, -1, -1,0, -1, 1, 0]]);
	end if;
	V:=QuadraticSpace(J1);
	sol1, kern1:= Solution(J*Transpose(Matrix(vec1)), Vector([GF(2)!1]));
        sol2, kern2:= Solution(J*Transpose(Matrix(vec2)), Vector([GF(2)!1]));
	for k in kern1 do
		so:=sol1+k;
		if (so*J1*Transpose(Matrix(so)))[1] eq 0 then
	   	   sol1:= so;
		   break;
		end if;
	end for;
	for k in kern2 do
                so:=sol2+k;
                if (so*J1*Transpose(Matrix(so)))[1] eq 0 then
                   sol2:= so;
                   break;
                end if;
        end for;
        ortho1:=OrthogonalComplement(V, sub<V|[vec1, sol1]>);
        ortho2:=OrthogonalComplement(V, sub<V|[vec2, sol2]>);
        bool, Trafo := IsIsometric(ortho2, ortho1);
	Mat1:=Matrix([vec1, sol1] cat[Trafo(b): b in Basis(ortho2)]  );
	Mat2:=Matrix([vec2, sol2] cat Basis(ortho2) );
        return liftSpN(Mat1^(-1) * Mat2, 3);
end function;


procedure TransformThetanulls(~thetas, S)
	CC:=Parent(thetas[1]);
        I:=CC.1;
	g:=Nrows(S) div 2;
        A:=Submatrix(Transpose(S), [1..g], [1..g]);
        B:=Submatrix(Transpose(S), [1..g], [g+1..2*g]);
        C:=Submatrix(Transpose(S), [g+1..2*g], [1..g]);
        D:=Submatrix(Transpose(S), [g+1..2*g], [g+1..2*g]);
        even:=EvenThetaCharacteristics(g);
        thetatrafo:=[CC!0: i in [1..2^(2*g)]];
        for eve in even do
                diag1:=Diagonal(B*Transpose(A));
                diag2:=Diagonal(D*Transpose(C));
                trafo:=Eltseq(Matrix([eve[1] cat eve[2]])*S+ Matrix(1,2*g, diag1 cat diag2));
                trafo:=[trafo[1..g], trafo[g+1..2*g]];
                exp:=&+[(trafo[1][i]-diag1[i]  ) * (trafo[2][i] +diag2[i])-eve[1][i]*eve[2][i]: i in [1..g]] div 2;
                carry:=&+[ trafo[1][i]*(trafo[2][i] div 2): i in [1..g]];
                thetatrafo[TCharToIndex(trafo)]:= (I)^exp * (-1)^carry * thetas[TCharToIndex(eve)];
        end for;
        thetas:=thetatrafo;
end procedure;

function ComputeSquareRootOnCone(f6)
  R<z,x,y,w> := Parent(f6);
  CC<I> := BaseRing(R);
  RQ<Z,X,Y,W> := PolynomialRing(Rationals(),4);
  S, mp := quo< RQ | X*Y - Z^2 >;
  
  mons6 := [mp(el) : el in MonomialsOfDegree(RQ,6) | not IsDivisibleBy(el, Z^2)];
  mons3 := [mp(el) : el in MonomialsOfDegree(RQ,3) | not IsDivisibleBy(el, Z^2)];
  monsred := [el : el in MonomialsOfDegree(R,6) | IsDivisibleBy(el, z^2)];
  
  //TODO: Make reduction function for cleaner code
  //Reduce f6 mod xy-z^2
  for m in monsred do
    coeff := MonomialCoefficient(f6, m);
    n := Degree(m, z);
    k := n div 2;
    f6 := f6 - coeff*m + coeff * (m div z^(2*k)) * (x*y)^k;
  end for;
  
  prods := AssociativeArray();
  for i := 1 to #mons3 do
    m1 := mons3[i];
    for j := i to #mons3 do
      m2 := mons3[j];
      if IsDefined(prods, m1*m2) then
        prods[m1*m2] := Append(prods[m1*m2], [m1,m2]);
      else
        prods[m1*m2] := [[m1,m2]];
      end if;
    end for;
  end for;

  known := [mp(W^3)];
  h:= hom< RQ -> R | z, x, y, w>;
  f3 := Sqrt(CC!MonomialCoefficient(f6,h(known[1]^2)))*h(known[1]); // use @@ ?
  while #known ne #mons3 do
    for k->v in prods do
      all_pairs := &cat v;
      unknown := [el : el in all_pairs | not el in known];
      if #unknown eq 1 then
        partner := [el : el in [pair : pair in v | unknown[1] in pair][1] | el ne unknown[1]];
        if partner[1] eq known[1] then
          Append(~known, unknown[1]);
          //printf "now %o known\n", unknown[1];
          //printf "now #known = %o\n", #known;
          Remove(~prods, k);
          c := CC!MonomialCoefficient(f6,h(k));
          for pair in [pair : pair in v | not unknown[1] in pair] do
            if pair[1] eq pair[2] then
              c -:= MonomialCoefficient(f3,h(pair[1]))*MonomialCoefficient(f3,h(pair[2]));
            else
              c -:= 2*MonomialCoefficient(f3,h(pair[1]))*MonomialCoefficient(f3,h(pair[2]));
            end if;
          end for;
          c /:= 2*MonomialCoefficient(f3,h(known[1]));
          f3 +:= c*h(unknown[1]);
        end if;
      end if;
    end for;
  end while;
  return f3;
end function;

function ComputeCurveVanTheta0(thetas, v)
	g:=4;
	CC:=Parent(thetas[1]);
	S:=MapSympl(Vector([GF(2)! el : el in v[1] cat v[2]]));
        TransformThetanulls(~thetas, S);
        tritangents := ComputeTritangents(thetas);
        bitangents := OddThetaPrymInfoHyp(thetas);
	/*TODO: Find the bijection between the characteristics and the odd theta lines. The lines are ordered by the natural ordering coming from the Weierstrass points.
The list of characteristics is:
	  [
    (1 1 1 0 0 1),
    (0 0 1 0 0 1),
    (0 0 1 0 1 1),
    (0 1 0 0 1 1),
    (0 1 0 0 1 0),
    (0 1 1 0 1 0),
    (0 1 1 1 0 1),
    (0 0 1 1 0 1),
    (0 0 1 1 1 1),
    (0 1 0 1 1 1),
    (0 1 0 1 1 0),
    (0 1 1 1 1 0),
    (0 1 1 0 0 1),
    (1 1 1 1 1 1),
    (1 0 0 1 1 1),
    (1 0 0 1 1 0),
    (1 0 1 1 1 0),
    (1 0 1 0 0 1),
    (1 0 0 1 0 1),
    (1 0 0 1 0 0),
    (1 0 1 1 0 0),
    (1 0 1 0 1 1),
    (1 1 1 1 0 0),
    (1 1 0 1 0 0),
    (1 1 0 0 1 1),
    (1 1 0 1 0 1),
    (1 1 0 0 1 0),
    (1 1 1 0 1 0)
]


[[ 1, 1, 0,  1, 0, 0 ], 
    [  1, 0, 0,  1, 0, 0 ],
    [  1, 0, 1,  1, 0, 0 ],
    [  1, 1, 1,  1, 0, 0 ],
    [  1, 0, 1,  1, 1, 0 ],
    [  1, 0, 0,  1, 1, 0 ],
    [  1, 0, 0,  1, 1, 1 ],
    [  1, 1, 1,  1, 1, 1 ],
    [  1, 0, 0,  1, 0, 1 ],
    [  1, 1, 0,  1, 0, 1 ]];

    So should be that:
*/
        bitangents:=[bitangents[i]: i in  [24, 20, 21, 23, 17, 16, 15, 14, 19, 26 ]];

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
       DN:=SingularValueDecomposition(N);
       gammaiinv:=NumericalKernel(N: Epsilon:=RR!(10^(-15)));

       DXnew, U, V:=SingularValueDecomposition(Xnew);
       Upart:=Matrix(U[1..7]);
       phi:=Upart*DiagonalMatrix(Eltseq(gammaiinv))*fsq_mat;

  //Kernel is not deterministic

       Qpre:=Kernel(phi);
       Qpre1:=Qpre*Upart;
       Qnew:=&+[Eltseq(Basis(Qpre1)[1])[i]*mats1new[i]: i in [1..r] ];
       dualelt:=mats1newx*ChangeRing(Transpose(Upart), CC4);
       D, U, V:=SingularValueDecomposition(phi);
       phiext:=HorizontalJoin(phi, Matrix(7,1, Eltseq(Conjugate(U)[7])));
       if IsVerbose("User1", 1) then
          Dfsqmat:=SingularValueDecomposition(fsq_mat); 
          Dphiext:=SingularValueDecomposition(Transpose(phiext));

          print "SVs of fsq_mat", Diagonal(Dfsqmat),"\n SVs of Xnew=", Diagonal(DXnew), "\n SVs of N=", Diagonal(DN), "\n SVs of phi=", Diagonal(D), "\n SVs of phiext=", Diagonal(Dphiext);
        end if;


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
        for i in [1..3] do
          S2[i,i] := Sqrt(S2[i,i]);
        end for;
        S2[4,4] := 1;

        I:=CC.1;
  //Coordinate Transformation that maps x^2 + y^2 +z^2 to xy - z^2.
        DiagToCone := Matrix(CC, 4, 4, [[1,0,0,0],[0,1/2*I,-1/2,0],[0,1/2*I,1/2,0], [0,0,0,1]]);

  //Complete Transformation
        QtoCone := DiagToCone * S2^(-1) * S; 


        v:=Matrix(CC4, [[CC4.1, CC4.2, CC4.3, CC4.4]]);

  //Apply coordinate transformation to detqdual to map the quadric to a cone
  detqdualoncone := Evaluate(detqdual, Eltseq((v * ChangeRing(QtoCone, CC4))[1]));
  
  //Pull back to P1 x P1 and take the square root there
  ConeCubic := ComputeSquareRootOnCone(detqdualoncone);
 
  //Reverse the coordinate transformation
  cubic := Evaluate(ConeCubic, Eltseq((v * ChangeRing(QtoCone^(-1), CC4))[1]));
  quadric:=(v*ChangeRing(Qnew, CC4) *Transpose(v))[1,1];
  return [quadric, cubic];

end function;

function ComputeCurveGeneric(thetas)
  tritangents := ComputeTritangents(thetas);
  bitangents := ComputeBitangents(thetas);
  if IsVerbose("User1", 1) then
    print "\n tritangents:", tritangents;
    print "\n bitangents:", bitangents;
  end if;
  quadric, cubic := ComputeCurve([bitangents[i]: i in [10, 23, 4, 20,  17, 9, 12,  1, 5, 11]], tritangents);
  return [quadric, cubic];
end function;

function ComputeCurveHypEll(thetas, v0s)
  gamma := ComputeGamma(v0s);
  eta := EtaFunction(gamma);
  rosens := [ TakaseQuotient([theta^2 : theta in thetas], eta, 1, l, 2) : l in [3..9] ];
  return rosens;
end function;

intrinsic ReconstructCurveG4(tau::AlgMatElt)->SeqEnum
{}
  thetas := ComputeThetas(tau);
  return ReconstructCurveG4(thetas);
end intrinsic;

intrinsic ReconstructCurveG4(thetas::SeqEnum)->SeqEnum
{}
  v0s := FindDelta(thetas);
  NrOfZeros := #v0s;
  if NrOfZeros eq 0 then
    return ComputeCurveGeneric(thetas);
  end if;
  
   if NrOfZeros eq 1 then
    return ComputeCurveVanTheta0(thetas, v0s[1]);
  end if;
  
  if NrOfZeros eq 10 then
    return ComputeCurveHypEll(thetas, v0s);
  end if;
  
  error("Wrong number of even theta characteristics are zero. ");
end intrinsic;


intrinsic RationalReconstructCurveG4(Pi::AlgMatElt)->SeqEnum
{}
  QQ := Rationals();
  Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
  tau := Pi1^-1*Pi2;
  thetas := ComputeThetas(tau);
  quadric, cubic := Explode(ReconstructCurveG4(thetas));
  CC4 := Parent(quadric);
  CC := BaseRing(CC4);
  X:=Matrix(CC4, 4,1, [CC4.i: i in [1..4]]);

  tritangentbasis := [
    [GF(2)|1, 1, 1, 0, 1, 1, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 1, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 1, 1, 0],
    [GF(2)|0, 1, 1, 0, 0, 1, 0, 0]];

  TTB:=[];
  
  for c in tritangentbasis do
    chara := [Integers()!v : v in Eltseq(c)];
    chara := [chara[1..4], chara[5..8]];
    Append(~TTB, TritangentPlane(Pi, chara));
  end for;

  TtoS := Matrix(TTB[1..4]);
  D := DiagonalMatrix(Eltseq(Vector(TTB[5]) * (TtoS)^-1));
  M := TtoS^-1 * D^-1;

  quadric:=Evaluate(quadric, Eltseq(ChangeRing(Inverse(M),CC4)*X));
  quadric_C:=quadric/LeadingCoefficient(quadric);
  cubic_C:=Evaluate(cubic, Eltseq(ChangeRing(Inverse(M), CC4)*X));
  
  R<x,y,z,w> := PolynomialRing(QQ, 4);
  mons2 := MonomialsOfDegree(R,2);
  mons3 := MonomialsOfDegree(R,3);
  
  h := hom< R -> CC4 | x, y, z, w>;
  
  quadric_Q:= R!0;
  
  for m in mons2 do
    coeff_C := MonomialCoefficient(quadric_C, h(m));
    coeff_Q := BestApproximation(Real(coeff_C), 1000000);
    quadric_Q +:= coeff_Q * m;
  end for;
  
  V := VectorSpace(QQ, #mons3);
  U := sub< V |[Vector([MonomialCoefficient(quadric_Q *R.i, m) : m in mons3]) : i in [1..4]]>; 
  B := ExtendBasis(U, V);
  Uc := B[5..#mons3];
  
  MB := Matrix(B);
  MUc := Matrix(Uc);
  
  v := Vector([MonomialCoefficient(cubic_C, h(m)) : m in mons3]); 
  w := Vector(Eltseq((v * ChangeRing(MB^(-1), CC)))[5..#mons3]);
  w := (w * ChangeRing(MUc, CC));
  
  scalar, i := Maximum([Abs(k) : k in Eltseq(w)]);
  
  w := w/w[i];
  
  cubic_C := &+[w[i] * h(mons3[i]) : i in [1..#mons3]];
  cubic_Q:= R!0;
  
  for m in mons3 do
    coeff_C := MonomialCoefficient(cubic_C, h(m));
    coeff_Q := BestApproximation(Real(coeff_C), 1000000);
    cubic_Q +:= coeff_Q * m;
  end for;
  
  return [quadric_Q, cubic_Q];
end intrinsic;



