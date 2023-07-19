/*
R<y0,y1,y2,y3,y4> := PolynomialRing(QQ,5);
f := (y0*y1 + y0*y2 + y1*y2 - y3*y4)^2 - 4*y0*y1*y2*(y0 + y1 + y2 + y3 + y4);
*/

// Thm 10.3.18, p. 543

function EtaFunction0(g);
ZZ:=Integers();
zer:=ZeroMatrix(ZZ, g,g);
zer1:=ZeroMatrix(ZZ, g,1);
id:=IdentityMatrix(ZZ,g);
        triang:=zer;
        for i in [1..g] do
                for j in [i..g] do
                        triang[i,j]:=1;
                end for;
        end for;
        M:=VerticalJoin(HorizontalJoin(id, zer1), HorizontalJoin(triang, zer1));
        N:=VerticalJoin(HorizontalJoin(id, zer1), HorizontalJoin(zer1, triang));
        return &cat[[Transpose(Matrix(Transpose(N)[i])), Transpose(Matrix(Transpose(M)[i])) ]: i in [1..g+1]];
end function;




/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Integers(), [[ 0 : i in [1..2*g] ]]));
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;


F<nu> := CyclotomicField(5);
R<x> := PolynomialRing(F);
f := &*[x-i: i in [0..5]];
roots := [el[1] : el in Roots(f)];
U:={1,3,5};
eta := EtaFunction0(2);

theta4 := [F!0: i in [1..16]];
for i in [1..5] do
	for j in [i+1..5] do
	T := {i,j, 6};	
	Tcomp := {1..6} diff T;
        S := T sdiff U;
	sign := (-1)^#(S meet U);
	cha := EtaValue(eta, Setseq(S));
	theta4[TCharToIndex(cha)] := sign * &*[ &*[(roots[nu]-roots[mu])^(-1): nu in T]: mu in Tcomp];
	end for;
end for;
	
M := ZeroMatrix(Rationals(), 5, 16);
M[1, 8]:=1;
M[1, 12]:=-1;
M[1, 15]:=1;
M[1, 9]:=-1;
M[2, 16]:=1;
M[2, 1]:=-1;
M[2, 6]:=-1;
M[2, 12]:=-1;
M[3, 6]:=1;
M[3, 2]:=-1;
M[3, 15]:=-1;
M[3, 3]:=1;
M[4, 4]:=1;
M[4, 16]:=-1;
M[4, 9]:=1;
M[4, 3]:=-1;
M[4, 3]:=1;
M[5, 4]:=1;
M[5, 8]:=-1;
M[5, 1]:=1;
M[5, 2]:=-1;
Ech := EchelonForm(M);
pivots := [Min([i: i in [1..16]| Ech[j,i] ne 0  ]): j in [1..5]];
even := [ 1, 2, 3, 4, 6, 8, 9, 12, 15, 16 ];
R16 := PolynomialRing(Rationals( ), 16);
t8 := [R16.i^2: i in even];
t16 := [R16.i^4: i in even];
equ := t8^2-4*t16;
for i in [1..5] do
piv:= pivots[i];
subs[piv]:= -&+[Ech[i, j]*R16.j: j in [piv+1..16]];
end for;
equnew := Evaluate(equ, subs);
R4 := PolynomialRing(Rationals(),5);
subsfin := [R4!0: i in [1..16]];
subsfin[6] := R4.1;
subsfin[9] := R4.2;
subsfin[12] := R4.3;
subsfin[15] := R4.4;
subsfin[16] := R4.5;
equfin:=Evaluate(equnew, subsfin);

for loo in [0..7] do
bin:=Intseq(loo,2,3);
rostest:= [ros[water]*(-1)^bin[water]: water in [1..3]];
rootstest := [0,1] cat [-ro/(ro-2): ro in rostest] cat [-1];
ftest:=x*(x-1)*&*[(x-Rationals()!ro): ro in rostest];
//print IgusaInvariantsEqual(IgusaInvariants(f), IgusaInvariants(ftest));
theta4test:= [Rationals()!0: i in [1..16]];
for i in [1..5] do
        for j in [i+1..5] do
        T := {i,j, 6};
        Tcomp := {1..6} diff T;
        S := T sdiff U;
        sign := (-1)^#(S meet U);
        cha := EtaValue(eta, Setseq(S));
        theta4test[TCharToIndex(cha)] := sign * &*[ &*[(rootstest[nu]-rootstest[mu])^(-1): nu in T]: mu in Tcomp];
        end for;
end for;
print &and[theta4test[i]*theta42[1] eq theta42[i]*theta4test[1]: i in [1..16]];
end for;

/*
 * Two gluable curves:
 * y^2= x^6 - 15*x^5 + 85*x^4 - 225*x^3 + 274*x^2 - 120*x
 * y^2=x^5 - 5/3*x^4 + 35/48*x^3 - 5/96*x^2 - 1/96*x
 take the 7th maximal isotropic subgroup
 * */


/*
SpecialPlanes(roots);
[[MonomialCoefficient(p,m) : m in MonomialsOfDegree(Parent(p),1)] : p in planes];
*/

// copypasta
/*
F<nu> := CyclotomicField(5);
R<x> := PolynomialRing(F);
f := x^6-x;
roots := [el[1] : el in Roots(f)];
AddTwoTorsionPoint(f, roots);
mp_two := $1;
AddTwoTorsionPoint(f, roots[1..2]);
ts := roots[1..2];
  t1, t2 := Explode(ts);
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
g := f;
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
TwoTorsionMatrix(G,H);
mat_two := $1;
DiagonalForm(mat_two);
DiagonalForm;
Diagonalization;
DiagonalisingMatrix;
CharacteristicPolynomial(mat_two);
Factorization($1);
mat_two;
mat_two^2;
MinimalPolynomial(mat_two);
DiagonalisingMatrix(mat_two);
JordanForm(mat_two);
M := mat_two;
D, P := JordanForm(M);
P;
Determinant(P);
P*M*(P^-1);
Binomial(6,2);
M1 := M;
ts;
ts1 := ts;
ts2 := roots[3..4];
ts := ts2;
t1, t2 := Explode(ts);
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
M2 := TwoTorsionMatrix(G,H);
M1;
M2;
M2^2;
M1*M2 - M2*M1;
P*M2*(P^-1);
IsDiagonal($1);
DiagonalisingMatrix;
Eigenspaces(M1);
Eigenspace;
Eigenvalues(M1);
D1 := D;
P1 := P;
D2, P2 := JordanForm(M2);
P2*M1*(P2^-1);
P1*M2*(P1^-1);
Diagonalization([M1, M2]);
Ds, Ps := Diagonalization([M1, M2]);
Parent(Ps[1]);
BaseRing($1);
#Ps;
P := Ps;
Parent(P);
#Ds;
ts3 := [roots[1], roots[3]];
ts := ts3;
assert #ts eq 2;
  t1, t2 := Explode(ts);
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
M3 := TwoTorsionMatrix(G,H);
M3;
M3^2;
P;
Rows(Transpose(P))[1];
P*M1*(P^-1);
IsDiagonal($1);
Rows(Transpose(P^-1))[1];
v1 := $1;
v2 := v1*Transpose(M3);
v2 := v1*Transpose(ChangeRing(M3,BaseRing(Parent(v1)));
v2 := v1*Transpose(ChangeRing(M3,BaseRing(Parent(v1))));
ts1;
ts2;
ts3;
ts := [roots[1], roots[5]];
ts3 := [roots[1], roots[5]];
assert #ts eq 2;
  t1, t2 := Explode(ts);
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
M3 := TwoTorsionMatrix(G,H);
v2 := v1*Transpose(ChangeRing(M3,BaseRing(Parent(v1))));
v1;
v2;
ts4 := [roots[3], roots[6]];
ts := ts4;
assert #ts eq 2;
  t1, t2 := Explode(ts);
  R<x> := Parent(g);
  k := BaseRing(R);
  G := (x-t1)*(x-t2);
  assert g mod G eq 0;
  H := g div G;
M4 := TwoTorsionMatrix(G,H);
v3 := v1*Transpose(ChangeRing(M4,BaseRing(Parent(v1))));
v4 := v1*Transpose(ChangeRing(M3*M4,BaseRing(Parent(v1))));
M3*M4 - M4*M3;
A := Transpose(Matrix([Eltseq(el) : el in [v1,v2,v3,v4]]));
A;
K := CalculateKummer(f);
K;
fK := DefiningEquation(K);
E := BaseRing(Parent(A));
E;
fKE := ChangeRing(fK,E);
fKE;
BaseChange(A,Parent(fKE));
BaseChange(A,MatrixAlgebra(4,Parent(fKE)));
ChangeRing(A,Parent(fKE));
A := $1;
MonomialsOfDegree(Parent(fKE),1);
SetToSequence($1);
var_mat := Matrix(4,1,$1);
A*var_mat;
new_vars := $1;
Evaluate(fKE, new_vars);
Evaluate(fKE, Eltseq(new_vars));
f_new := $1;
Monomials(f_new);
*/
