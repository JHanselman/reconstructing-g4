AttachSpec("../../magma/spec");
prec := 30;
SetDefaultRealFieldPrecision(prec);
QQ:=Rationals();
R<x,y,z,w> := PolynomialRing(QQ,4);
Q := x^2 + 2*y^2 + z^2 - w^2;
F := 3*x^3 + y^3 - 7*z^3 + w^3;
C := Curve(Proj(R),[F,Q]);
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
assert IsIrreducible(p);
S := RiemannSurface(p : Precision := prec);
g := Genus(S);
Pi_big := BigPeriodMatrix(S);
Pi1, Pi2 := SplitBigPeriodMatrix(Pi_big);
tau := SmallPeriodMatrix(S);
CC:=BaseRing(tau);
RR:=RealField(prec);
Quadric:= Matrix(CC, 4,4, [[1,0,0,4/7],[0,4/7,0,0], [0,0, 4/7, 0], [4/7, 0,0,2/7]]);

ZZ := Integers();

T:=Time();
v := Vector([GF(2) | 0, 0, 0, 0, 1, 0, 0, 0 ]);
steiner := [];
steinerrie:=[];
//1st 7 form an aronhold set.
aronhold := [
  [[1,1,1],[1,1,1]],[[0,0,1],[0,1,1]],[[0,1,1],[0,0,1]],[[1,0,1],[1,0,0]],[[1,0,0],[1,0,1]],[[1,1,0],[0,1,0]],[[0,1,0],[1,1,0]], [[0,1,0],[0,1,0]], [[1,0,0],[1,1,0]], [[1,1,0],[1,0,0]]
  ];
rie := [
  [[0,0,0],[0,0,1]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,1]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]], [[0,0,1],[1,1,0]],[[0,0,0],[1,1,0]],[[0,0,0],[0,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[1,0,0]]
  ];
r:=10;



for aro in aronhold[1..r] do

    cc1 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!0] cat [GF(2)!el: el in  aro[2]];
    cc2 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!1] cat [GF(2)!el: el in  aro[2]];
    c1:=[[0] cat aro[1] , [0] cat aro[2]];
    c2:=[[0] cat aro[1] , [1] cat aro[2]];

    Append(~steiner, [c1,c2]);


end for;


function parity_char(chars)
  zer:=ZeroMatrix(GF(2), 4,4);
  id:=IdentityMatrix(GF(2),4);
  J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
  return (chars*J1*Transpose(chars))[1,1];
end function;

function check_azygetic(chars)
  nchar := #chars;
  i:=1;
  for j in [i+1..nchar] do
  for k in [j+1..nchar] do
    chari:=Matrix(GF(2), 1, 8, chars[i]);
    charj:=Matrix(GF(2), 1, 8, chars[j]);
    chark:=Matrix(GF(2), 1, 8, chars[k]);
    if parity_char(chari)+parity_char(charj)+parity_char(chark)+parity_char(chari+charj+chark) eq 0 then
       return false;
    end if;
  end for;
  end for;
  return true;
end function;

function DotProductSeq(v1,v2)
  assert #v1 eq #v2;
  return &+[v1[i]*v2[i] : i in [1..#v1]];
end function;


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
    [GF(2)|0, 1, 1, 0, 0, 1, 0, 0]
];

/* Code for comparison
TT2 := [];
for ca in tritangentsys do
  pair := [];
  for c in ca do
    chara := [ZZ!v : v in Eltseq(c)];
    chara := [chara[1..4], chara[5..8]];
    Append(~pair, TritangentPlane(Pi_big, chara));
  end for;
  Append(~TT2, pair);
end for;

TTB := [];
for c in tritangentbasis do
  chara := [ZZ!v : v in Eltseq(c)];
  chara := [chara[1..4], chara[5..8]];
  Append(~TTB, TritangentPlane(Pi_big, chara));
end for;

TtoS := Matrix(TTB[1..4]);
D := DiagonalMatrix(Eltseq(Vector(TTB[5]) * (TtoS)^-1));
M := TtoS^-1 * D^-1;

*/

V:=VectorSpace(GF(2), 8);
Laz:=[];
for v in V do
  if (check_azygetic(tritangentbasis cat [Eltseq(v)]) and (DotProductSeq(Eltseq(v)[1..4], Eltseq(v)[5..8]) eq 1)) then
    Append(~Laz, v);
  end if;
end for;

function Tritangents(tau)
  CC := BaseRing(tau);
  theta_list := AssociativeArray();
  tritangents := [[[CC!0,0,0,0], [CC!0,0,0,0]]: t in [1..10]];
  constant := [CC!0,0,0,0];
  for k in [1..4] do
    temp_list := [ Vector(t) : t in (tritangentbasis[1..k - 1] cat [tritangentbasis[5]] cat tritangentbasis[k+1 .. 4])];
    S1, S2, A := special_fundamental_system(temp_list);
    signs := signs_in_derivative_formula(A);
    T1 := CC!1;
    T2 := CC!1;
    for c in S1 do
      test := IsDefined(theta_list, c);
      if not test then
        chara := [ZZ!v : v in Eltseq(c)];
        chara := [chara[1..4], chara[5..8]];
        theta_list[c] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
      end if;
      T1 *:= theta_list[c];
    end for;

    for c in S2 do
      test := IsDefined(theta_list, c);
      if not test then
        chara := [ZZ!v : v in Eltseq(c)];
        chara := [chara[1..4], chara[5..8]];
        theta_list[c] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
      end if;
      T2 *:= theta_list[c];
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
        for c in S1 do
          test := IsDefined(theta_list, c);
          if not test then
            chara := [ZZ!v : v in Eltseq(c)];
            chara := [chara[1..4], chara[5..8]];
            theta_list[c] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
          end if;
          T1 *:= theta_list[c];
        end for;

        for c in S2 do
          test := IsDefined(theta_list, c);
          if not test then
            chara := [ZZ!v : v in Eltseq(c)];
            chara := [chara[1..4], chara[5..8]];
            theta_list[c] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
          end if;
          T2 *:= theta_list[c];
        end for;

        tritangents[i][j][k] := (signs[1]*T1 + signs[2]*T2)/constant[k];

      end for;
    end for;
  end for;

  return tritangents;
end function;

function JacobiNullwerte(tau, odd_thetas)
  g := Nrows(tau);
  CC<I> := BaseRing(Parent(tau));
  prec := Precision(CC);
  pi := Pi(CC);
  dets := [];
  for k in [1..g] do
    derivs := [];
    temp_list := [ Vector(t) : t in (odd_thetas[1..k-1] cat [odd_thetas[g+1]] cat odd_thetas[k+1..g])];
    for c in temp_list do
      chara := [ZZ!v : v in Eltseq(c)];
      chara := [chara[1..g], chara[g+1..(2*g)]];
      for i := 1 to g do
        dz := [0 : j in [1..g]];
        dz[i] := 1;
        Append(~derivs, Theta([CC!0 : j in [1..g]], tau : char := chara, dz := [dz], prec := prec));
      end for;
    end for;
    M := Matrix(g,g,derivs);
    Append(~dets, pi^(-g)*Determinant(M));
  end for;
  return dets;
end function;

