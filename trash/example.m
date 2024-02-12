AttachSpec("../../../CHIMP/CHIMP.spec");
AttachSpec("../../magma/spec");
load "gluingfuncs.m";


QQ := RationalsExtra(20);
R<x> := PolynomialRing(QQ);

X1:= HyperellipticCurve(x^6 - 15*x^5 + 85*x^4 - 225*x^3 + 274*x^2 - 120*x);  
X2:= HyperellipticCurve(x^5 - 5/3*x^4 + 35/48*x^3 - 5/96*x^2 - 1/96*x);  

Vs := AllVs2For22();

V:= Vs[7];

P1 := PeriodMatrix(X1); P2 := PeriodMatrix(X2);
P := DiagonalJoin(P1, P2);

Q := QFromPVFor22(P, V);
tau := SmallPeriodMatrix(Q);
T1 := ComputeThetas(tau/2);
T2 := ComputeThetas_old(tau);

CC:= Parent(Q[1,1]);

//Change polarization convention
Q := Q* DiagonalJoin(IdentityMatrix(CC, 4), -IdentityMatrix(CC, 4)) ;

eqs := RationalReconstructCurveG4(Q);  
eqs;

P3 := ProjectiveSpace(Parent(eqs[1]));
C := Curve(Scheme(P3, eqs));
p:= 17;
LPolynomial(C, p);

LPolynomial(Curve(Reduction(X1, p)));
LPolynomial(Curve(Reduction(X2, p)));





