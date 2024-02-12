AttachSpec("../magma/spec");
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

Q1, C1 := ReconstructCurveG4(tau);



