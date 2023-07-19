AttachSpec("../../../CHIMP/CHIMP.spec");
AttachSpec("../../magma/spec");
load "gluingfuncs.m";

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



// make curves to glue using formulas in Howe--Leprevost--Poonen
QQ := RationalsExtra(100);
//QQ := RationalsExtra(20);
//QQ := Rationals();
K<t,u> := RationalFunctionField(QQ,2);
bt := (16*t^3+16*t^2+6*t+1)/(8*t^2-1)^2;
ct := (16*t^3+16*t^2+6*t+1)/(2*t*(4*t+1)*(8*t^2-1));
bu := (16*u^3+16*u^2+6*u+1)/(8*u^2-1)^2;
cu := (16*u^3+16*u^2+6*u+1)/(2*u*(4*u+1)*(8*u^2-1));

Eainvs := [1-ct, -bt, -bt, 0, 0];
Fainvs := [1-cu, -bu, -bu, 0, 0];
E := EllipticCurve(Eainvs);
F := EllipticCurve(Fainvs);

f := HyperellipticPolynomials(WeierstrassModel(E));
g := HyperellipticPolynomials(WeierstrassModel(F));
alpha1, alpha2, alpha3 := Explode([Roots(el[1])[1][1] : el in Factorization(f)]);
beta1, beta2, beta3 := Explode([Roots(el[1])[1][1] : el in Factorization(g)]);
//beta2, beta3, beta1 := Explode([Roots(el[1])[1][1] : el in Factorization(g)]);

Deltaf := Discriminant(f);
Deltag := Discriminant(g);

a1 := (alpha3-alpha2)^2/(beta3-beta2) + (alpha2-alpha1)^2/(beta2-beta1) + (alpha1-alpha3)^2/(beta1-beta3);
b1 := (beta3-beta2)^2/(alpha3-alpha2) + (beta2-beta1)^2/(alpha2-alpha1) + (beta1-beta3)^2/(alpha1-alpha3);
a2 := alpha1*(beta3-beta2) + alpha2*(beta1-beta3) + alpha3*(beta2-beta1);
b2 := beta1*(alpha3-alpha2) + beta2*(alpha1-alpha3) + beta3*(alpha2-alpha1);

A := Deltag*a1/a2;
B := Deltaf*b1/b2;

R<x> := PolynomialRing(K);
h := -(A*(alpha2-alpha1)*(alpha1-alpha3)*x^2 + B*(beta2-beta1)*(beta1-beta3))*(A*(alpha3-alpha2)*(alpha2-alpha1)*x^2 + B*(beta3-beta2)*(beta2-beta1))*(A*(alpha1-alpha3)*(alpha3-alpha2)*x^2 + B*(beta1-beta3)*(beta3-beta2));

H := HyperellipticCurve(h);

Ss<s> := PolynomialRing(QQ);

K2<T> := RationalFunctionField(QQ);
R2<U,Y> := PolynomialRing(K2,2);
C0 := Curve(Spec(R2),(8*T^2-1)*(8*T^2 + 8*T + 1)*Y^2 - (8*U^2-1)*(8*U^2 + 8*U + 1));
C := ProjectiveClosure(C0);
EC, mp := EllipticCurve(C, C![T,1,1]);

P := C![T,-1,1];
Q := EvaluateByPowerSeries(mp, P);
P2 := Points((2*Q) @@ mp)[1];
/*
ev := [T,Eltseq(P2)[1]];
Evaluate(ev,1);
*/
hev := [Evaluate(el, [T,Eltseq(P2)[1]]) : el in Coefficients(h)];
R3<z> := PolynomialRing(K2);
hev := R3!hev;



/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Integers(), [[ 0 : i in [1..2*g] ]]));
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;


roots := [el[1] : el in Roots(hev)];
U:={1,3,5};
eta := EtaFunction0(2);

theta4 := [K2!0: i in [1..16]];

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
t8 := &+[R16.i^2: i in even];
t16 := &+[R16.i^4: i in even];
equ := t8^2-4*t16;

subs := Setseq(MonomialsOfDegree(R16,1)); 

for i in [1..5] do
piv:= pivots[i];
subs[piv]:= -&+[Ech[i, j]*R16.j: j in [piv+1..16]];
end for;
equnew := Evaluate(equ, subs);

FP<t1, t2> := PolynomialRing(Rationals(), 2);
//FF<t1, t2> := RationalFunctionField(Rationals(), 2);

R4 := PolynomialRing(FP,5);
subsfin := [R4!0: i in [1..16]];
subsfin[6] := R4.1;
subsfin[9] := R4.2;
subsfin[12] := R4.3;
subsfin[15] := R4.4;
subsfin[16] := R4.5;
equfin:=Evaluate(equnew, subsfin);

Point :=[ theta4[j]: j in [6,9,12,15,16]];
Den := Lcm([Denominator(t): t in Point]);
Point := [Numerator(Den * p) : p in Point];
Point1 := [Evaluate(p, t1): p in Point];
Point2 := [Evaluate(p, t2): p in Point];

TSpace := &+[Evaluate(Derivative(equfin, i), Point1) * R4.i: i in [1..5]];  
Sols := Evaluate(TSpace, Point2);
Sc:=Scheme(AffineSpace(FP), Sols);
Scomp := IrreducibleComponents(Sc);
Scomp := [ReducedSubscheme(gh) : gh in Scomp];
Ps := PointSearch(Scomp[12], 100); 


P:=Ps[1];

s1 := P[1];
s2 := P[2];


hev1 := [Evaluate(el, s1) : el in Coefficients(hev)];
hev1 := Ss!hev1;
X1 := HyperellipticCurve(hev1);

hev2 := [Evaluate(el, s2) : el in Coefficients(hev)];
hev2 := Ss!hev2;
X2 := HyperellipticCurve(hev2);

Vs := AllVs2For22();

P1 := PeriodMatrix(X1); P2 := PeriodMatrix(X2);
P := DiagonalJoin(P1, P2);

Qs := [];

for V in Vs do
  Q := QFromPVFor22(P, V);
  Append(~Qs, Q);
end for;

good := [];
for i->Q in Qs do
  print i;
  tau := SmallPeriodMatrix(Q);
  err := Abs(SchottkyModularForm(tau : prec := 20));
  print err;
  if err lt 10^-10 then
    Append(~good, i);
  end if;
end for;


V:= Vs[311];

hev1 := ChangeRing(hev1, RationalsExtra(100));
hev2 := ChangeRing(hev2, RationalsExtra(100));
X1 := HyperellipticCurve(hev1);
X2 := HyperellipticCurve(hev2);

P1 := PeriodMatrix(X1); P2 := PeriodMatrix(X2);
P := DiagonalJoin(P1, P2);

Q := QFromPVFor22(P, V);
tau := SmallPeriodMatrix(Q);

CC:= Parent(Q[1,1]);

//Change polarization convention
Q := Q* DiagonalJoin(IdentityMatrix(CC, 4), -IdentityMatrix(CC, 4)) ;

eqs := RationalReconstructCurveG4(Q);  
eqs;
