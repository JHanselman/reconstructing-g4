AttachSpec("/home/hanselma/gluing-g4/magma/spec");
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

T:=Time();
v := Vector([GF(2) | 0, 0, 0, 0, 1, 0, 0, 0 ]);
steiner := [];
steinerrie:=[];
//1st 7 form an aronhold set.
aronhold:=[[[1,1,1],[1,1,1]],[[0,0,1],[0,1,1]],[[0,1,1],[0,0,1]],[[1,0,1],[1,0,0]],[[1,0,0],[1,0,1]],[[1,1,0],[0,1,0]],[[0,1,0],[1,1,0]], [[0,1,0],[0,1,0]], [[1,0,0],[1,1,0]], [[1,1,0],[1,0,0]]];
rie:=[[[0,0,0],[0,0,1]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,1]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]], [[0,0,1],[1,1,0]],[[0,0,0],[1,1,0]],[[0,0,0],[0,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[1,0,0]]];
r:=10;



for aro in aronhold[1..r] do

    cc1 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!0] cat [GF(2)!el: el in  aro[2]];
    cc2 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!1] cat [GF(2)!el: el in  aro[2]];
    c1:=[[0] cat aro[1] , [0] cat aro[2]];
    c2:=[[0] cat aro[1] , [1] cat aro[2]];

    Append(~steiner, [c1,c2]);


end for;

//Copied from Agostini for quick test

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


V:=VectorSpace(GF(2), 8);
Laz:=[];
for v in V do                                                
  if (check_azygetic(tritangentbasis cat [Eltseq(v)]) and (DotProductSeq(Eltseq(v)[1..4], Eltseq(v)[5..8]) eq 1)) then
    Append(~Laz, v);                                                     
  end if;
end for;


function ThetasNeeded()

  char_list := [];
  for i in [1..4] do
    
  end for;
end function;

function TritangentPlane2(tau, tritangentbasis, char)
  CC := BaseRing(Pi);
  cs := [];
  for i := 1 to 4 do
    dz := [0,0,0,0];
    dz[i] := 1;
    Append(~cs, Theta([CC | 0,0,0,0], tau : char := char, dz := [dz], prec := prec));
  end for;
  cs := Eltseq(Matrix(1,4,cs)*(Pi1^-1));
  return cs;
end function;

tritangents := [];
//for char in chars[1..4] do
//for char in chars do
for s in steiner do
  new := [];
  for char in s do
      print(char);
    Append(~new, TritangentPlane(Pi_big, char));
    //Append(~new, TritangentPlane(HorizontalJoin(IdentityMatrix(CC, 4), tau), char));
  end for;
  Append(~tritangents, new);
end for;

Time(T);
T:=Time();
chars_even := EvenThetaCharacteristics(3);
/*eps := chars_even[1];
eta_sqs := AssociativeArray(); // now just etas
for delta in chars_even do
  print delta;
  delta_new1 := [[0] cat el : el in delta];
  delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
  // formula from Lemma 1 (p. 14) of Farkas
  eta_sqs[delta] := Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2);  // TODO: fix signs here
end for;
*/

thetas := [];
for i := 1 to 64 do
    s := Intseq(i mod 64,2,6);
    s := Reverse(s);
    delta := [s[1..3], s[4..6]];
    if delta in chars_even then
        print delta;
        delta_new1 := [[0] cat el : el in delta];
        delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
      // formula from Lemma 1 (p. 148) of Farkas
        thetas[i] := Sqrt(Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2));
    else
        Append(~thetas, 0);
    end if;
end for;
thetas:=correct_signs(thetas);

Time(T);
T:=Time();

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

mods := ModuliFromTheta(thetas);

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

A := Matrix(3,3,[1/el : el in mods]);
Ainv:=A^-1;
Ainv*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);
lambdas := $1;

mods_mat := [[mods[i], mods[i+1], mods[i+2]] : i in [1,4,7]];
L:=DiagonalMatrix(Eltseq(lambdas));
B := Matrix(3,3,mods)*L;
Binv := Inverse(B);
ks:=Binv*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
bitangents := [];
   bitangents := [ [CC | 1, 0, 0], [CC | 0,1,0], [CC | 0,0,1], [CC | 1,1,1]];
bitangents cat:= mods_mat;
F, u0, u1, u2 := RiemannModelFromModuli(mods);
bitangents cat:= [Coefficients(el) : el in [u0, u1, u2]];
CC3<t0,t1,t2> := Parent(u0);
bitangents cat:= [Coefficients(el) : el in [t0+t1+u2, t0+u1+t2, u0+t1+t2]];

// (3)
for i := 1 to 3 do
  new := u0/mods_mat[1,1] + ks[i,1]*(mods_mat[2,i]*t1 + mods_mat[3,i]*t2);
  Append(~bitangents, Coefficients(new));
end for;
// (4)
for i := 1 to 3 do
  new := u1/mods_mat[2,1] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[3,i]*t2);
  Append(~bitangents, Coefficients(new));
end for;
// (5)
for i := 1 to 3 do
  new := u2/mods_mat[3,1] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[2,i]*t1);
  Append(~bitangents, Coefficients(new));
end for;
// (6)
for i := 1 to 3 do
  new := u0/(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i]) + u1/(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i]) + u2/(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]);
  Append(~bitangents, Coefficients(new));
end for;
for i := 1 to 3 do
  new := u0/(mods_mat[1,i]*(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + u1/(mods_mat[2,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + u2/(mods_mat[3,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
  Append(~bitangents, Coefficients(new));
end for;




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
/*sum_i gamma_i s_1i l_i l_i'= Q
sum_i gamma_i s_ji l_i l_i'= a_j Q

Compute gammai with Q:*/
mats1:=[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2]) * si[1][i]): i in [1..r]] cat [ ZeroMatrix(CC, 4,4): j in [1..sirows-1]];
mats1:=[(m +Transpose(m))/2 : m in mats1];
W:= Matrix( [&cat[[CC!Quadric[i,j]: j in [i..4]]: i in [1..4]] cat [CC!0: j in [1..10*(sirows-1)]  ]] );
mats2:=[[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2]) * si[j][i]): i in [1..r]]  cat [ Quadric: nu in [1..sirows-1]] : j in [2..sirows]];
mats2:=[[(m +Transpose(m))/2 : m in mats2[j]] : j in [1..sirows-1]];
mats:=[mats1] cat mats2;
X:=HorizontalJoin([Matrix(  [&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats[j]]  ): j in [1..sirows]]);
//gammai:=NumericalSolution(X, W: Epsilon:=RR!10^(-15));
//gammai:=Eltseq(gammai1/gammai1[1,1]);


//


//Compute gammai without Q:
//
//
CC4:=PolynomialRing(CC,4);
x:=Matrix(4,1,[CC4.i: i in [1..4]]);
mats1new:=[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2])): i in [1..r]];
mats1new:=[(m +Transpose(m))/2 : m in mats1new];
mats1newx:=Matrix(CC4, 1,r,[(Transpose(x)*ChangeRing(mats1new[i], CC4)*x)[1,1]: i in [1..r]]);

Xnew:=Matrix([&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats1new]  );
vi:=NumericalKernel(Xnew: Epsilon:=RR!10^(-15));
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
Time(T);



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


//For testing the correctness:
CC1<t>:=PolynomialRing(CC);
Ccan, map:=CanonicalImage(S);
Cplane:=Domain(map);
for i in [-10..-1] cat [1..10] do
	f1:=Evaluate(DefiningEquation(Cplane), [1/CC!i+CC.1, t]);
	ys:=[roo[1]: roo in Roots(f1)];
	coord:=[[Evaluate(DefiningEquations(map)[nu], [1/CC!i+CC.1, y]) : nu in [1..4]  ]: y in ys];
	print [Abs(Evaluate(detqdual, coo )): coo in coord];
end for;

S:= NormalForm(Qnew);

//COmpute the matrix S2 whose inverse will transform the normal form into one where all coefficients are 1.

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

//Apply coordinate transformation to detqdual
detqdualonsegre := Evaluate(detqdual, Eltseq((v * ChangeRing(QtoSegre, CC4))[1]));

//Homogeneous
P1P1<x1,x2,y1,y2> := PolynomialRing(CC,4);
x:=[x1,x2];
y:=[y1,y2];
f := Evaluate(detqdualonsegre, [x1*y1, x2*y2, x1*y2, x2*y1]);


xy:= AssociativeArray();
P3mons:= AssociativeArray();

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



v:=Matrix(CC4, [[CC4.1, CC4.2, CC4.3, CC4.4]]);

//Apply coordinate transformation to detqdual
cubic := Evaluate(SegreCubic, Eltseq((v * ChangeRing(QtoSegre^(-1), CC4))[1]));

//For testing the correctness:
for i in [-10..-1] cat [1..10] do
	f1:=Evaluate(DefiningEquation(Cplane), [1/CC!i+CC.1, t]);
	ys:=[roo[1]: roo in Roots(f1)];
	coord:=[[Evaluate(DefiningEquations(map)[nu], [1/CC!i+CC.1, y]) : nu in [1..4]  ]: y in ys];
	print [Abs(Evaluate(cubic, coo )): coo in coord];
end for;



