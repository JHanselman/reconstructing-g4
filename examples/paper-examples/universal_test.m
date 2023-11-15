
function QuarticFromAronhold(bits)
	R := Universe(bits[1]);
	F := FieldOfFractions(R);
	S := Matrix(F, bits[1..3]);
	Sinv := (S)^(-1);
	D0 := DiagonalMatrix(Eltseq(Matrix(F, 1,3, bits[4])* Sinv  ));
	mods_mat := Matrix([Matrix(F, 1, 3, bits[i]) * Sinv * (D0)^(-1) : i in [5..7]]);
	mods := Eltseq(mods_mat);
	A := Transpose(Matrix(3,3,[1/el : el in mods]));
	Ainv := A^-1;
	lambdas := Ainv*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);
	L := DiagonalMatrix(Eltseq(lambdas));
	B := Transpose(mods_mat) * L;
	Binv := B^(-1);
	ks := Binv*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
	bitangents := [];
	bitangents := [ [R | 1, 0, 0], [R | 0,1,0], [R | 0,0,1], [R | 1,1,1]];
	bitangents cat:= [mods[1..3], mods[4..6], mods[7..9]];

	a1 := mods[1];a2:=mods[2];a3:=mods[3];
	ap1 := mods[4];ap2:=mods[5];ap3:=mods[6];
	as1 := mods[7];as2:=mods[8];as3:=mods[9];
	P<x1,x2,x3> := PolynomialRing(F,3);
	M := Matrix([[F|1,1,1],[ks[1,1]*a1,ks[1,1]*a2,ks[1,1]*a3],[ks[2,1]*ap1,ks[2,1]*ap2,ks[2,1]*ap3]]);
	Mb := Matrix([[F|1,1,1],[1/a1,1/a2,1/a3],[1/ap1,1/ap2,1/ap3]]);
	U := -Mb^(-1)*M;
	u0 := U[1];
	u1 := U[2];
	u2 := U[3];
	u0 := u0[1]*x1+u0[2]*x2+u0[3]*x3;
	u1 := u1[1]*x1+u1[2]*x2+u1[3]*x3;
	u2 := u2[1]*x1+u2[2]*x2+u2[3]*x3;
	Quart := (x1*u0+x2*u1-x3*u2)^2-4*x1*u0*x2*u1;
	X := Eltseq(ChangeRing(D0 * S, P) * Matrix(P, 3, 1, [x1, x2, x3]));
	Quart := Evaluate(Quart, X);

	bitangents cat:= [[MonomialCoefficient(el, P.i): i in [1..3]] : el in [u0, u1, u2]];
  	
	R3<t0,t1,t2> := P;
	bitangents cat:= [[MonomialCoefficient(el, P.i): i in [1..3]] : el in [t0+t1+u2, t0+u1+t2, u0+t1+t2]];
  	mods_mat := Transpose(Matrix(mods_mat));
// (3)
  for i := 1 to 3 do
    new := u0/mods_mat[1,i] + ks[i,1]*(mods_mat[2,i]*t1 + mods_mat[3,i]*t2);
    Append(~bitangents, [MonomialCoefficient(new, P.i): i in [1..3]]);
  end for;
// (4)
  for i := 1 to 3 do
    new := u1/mods_mat[2,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[3,i]*t2);
    Append(~bitangents, [MonomialCoefficient(new, P.i): i in [1..3]]);

  end for;
// (5)
  for i := 1 to 3 do
    new := u2/mods_mat[3,i] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[2,i]*t1);
    Append(~bitangents, [MonomialCoefficient(new, P.i): i in [1..3]]);

  end for;


  // (6)
  for i := 1 to 3 do
    new := t0/((1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + t1/((1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + t2/((1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
    Append(~bitangents, [MonomialCoefficient(new, P.i): i in [1..3]]);
  end for;
  
  
  // (7)
  for i := 1 to 3 do
    new := u0/(mods_mat[1,i]*(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + u1/(mods_mat[2,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + u2/(mods_mat[3,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
    Append(~bitangents, [MonomialCoefficient(new, P.i): i in [1..3]]);
  end for;
	bitangents := [Eltseq(Matrix(1,3, bit)* D0 * S): bit in bitangents];
	return Quart, bitangents;
end function;
/*
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



den:=[Lcm([Denominator(bi) : bi in bitangents[i]]): i in [1..28]];
bitangents:=[[R!(bitangents[i][j]*den[i]): j in [1..3]]: i in [1..28]];

Paux<q200, q201,q202,q211,q212, q300,q301,q302,q311,q312,q322>:=PolynomialRing(P1,11);

q2:=Matrix(Paux, [[q200,q201, q202], [q201, q211, q212], [q202, q212,0]]  );
q3:=Matrix(Paux, [[q300,q301, q302], [q301, q311, q312], [q302, q312,q322]]);
X:=Matrix(Paux, 3,1, [x1,x2,x3]);
LHs:=(Transpose(X)*ChangeRing(q0, Paux)*X*Transpose(X)*q3*X-Transpose(X)*ChangeRing(q1, Paux)*X*Transpose(X)*q2*X)[1,1];
mons:=MonomialsOfDegree(P1, 4);
Mat:=Matrix([[MonomialCoefficient(MonomialCoefficient(LHs, Paux.i), m): i in [1..11]]  : m in mons]);
SubMat:=Submatrix(Mat,1,1,11,11);
Cof:=Matrix([[Cofactor(SubMat, i, j): i in [1..11]]: j in [1..11]]);
det:=Determinant(SubMat);
RHs:=Matrix([[MonomialCoefficient(Quart, m)]  : m in mons]);
SubRHs:=Submatrix(RHs,1,1,11,1);
qcoeff:=Cof*SubRHs;
gcdqcoeff:= Gcd(Eltseq(qcoeff));
qcoeff:= qcoeff div gcdqcoeff;
det := det div gcdqcoeff;
I:=Ideal(Eltseq(Mat*qcoeff-det*RHs));
print [Degree(elt) : elt in Eltseq(Mat*qcoeff-det*RHs)];
print Gcd( Eltseq(Mat*qcoeff-det*RHs));
S:=quo<R|I>;
SetOutputFile("GroebnerBasis.txt ");
SetVerbose("Groebner", 4);
T:=Time();GroebnerBasis(I); T1:=Time(T);
print T1;
q200v:=qcoeff[1,1];
q201v:=qcoeff[2,1];
q202v:=qcoeff[3,1];
q211v:=qcoeff[4,1];
q212v:=qcoeff[5,1];
q300v:=qcoeff[6,1];
q301v:=qcoeff[7,1];
q302v:=qcoeff[8,1];
q311v:=qcoeff[9,1];
q312v:=qcoeff[10,1];
q322v:=qcoeff[11,1];
q2:=Matrix(R, [[q200v,q201v, q202v], [q201v, q211v, q212v], [q202v, q212v,0]]  );
q3:=Matrix(R, [[q300v,q301v, q302v], [q301v, q311v, q312v], [q302v, q312v,q322v]]);
X:=ChangeRing(X, P1);

for i in [1..1] do
break;
	bit:=Matrix(R, 3,1, bitangents[i]);
	Lhs:=[bit*Matrix(RSpace(R,3).j): j in [1..3]];
	Lhs:=[l+Transpose(l): l in Lhs];
	qs:=[q0, q1, q2, q3];
	Mat:=Matrix([[Eltseq(l)[i]: i in [1,2,3,5,6,9]]: l in Lhs cat qs]  );
        H:=[(-1)^i*Minor(Mat, [1..i-1] cat [i+1..7], [1..6]): i in [4..7] ];
	bitres:=&+[H[i]*qs[i]: i in [1..4]];
	//bitres:=Transpose(X)*ChangeRing(bitres, P1)*X;
	//bitX:=&+[bitangents[i]*P1.i: i in [1..3]];
	//res:= bitres div bitX;
	//res:=Matrix(R, 3,1, [MonomialCoefficient(res, P1.i): i in [1..3]]);

end for;
*/

/*C:=HorizontalJoin([ChangeRing(q0, P1)*ChangeRing(X, P1), ChangeRing(q1, P1)*ChangeRing(X, P1),ChangeRing(q2, P1)*ChangeRing(X, P1),ChangeRing(q3, P1)*ChangeRing(X, P1) ]);
ci:=[(-1)^i * Minor(C, [1..3], [1..i-1] cat [i+1..4]): i in [1..4]];*/
