AttachSpec("../magma/spec");
load "../../Invariants_genus_4.m";
T:=Time();
SetDebugOnError(true);
SetVerbose("User1", 1);
ZZ:=Integers();
prec := 40;
SetDefaultRealFieldPrecision(prec);
QQ:=Rationals();
R<x,y,z,w> := PolynomialRing(QQ,4);
Q:=-10*x^2 - x*y + 8*x*z + 3*x*w - 9*y*z + y*w - 6*z^2 - 5*z*w - 5*w^2;
F:=-x^2*y - x^2*z - 5*x^2*w - 6*x*y^2 - 2*x*y*z - 9*x*y*w + 7*x*z^2 + 3*x*z*w -
    8*x*w^2 - 10*y^3 + 3*y^2*z - y^2*w - 3*y*z^2 - 7*y*w^2 - z^3 + z^2*w -
    10*z*w^2 + 3*w^3;
C := Curve(Proj(R),[F,Q]);
print C;
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
assert IsIrreducible(p);
S := RiemannSurface(p : Precision := prec);
g := Genus(S);
assert g eq 4;
Pi_big := BigPeriodMatrix(S);
Pi1, Pi2 := SplitBigPeriodMatrix(Pi_big);
tau := SmallPeriodMatrix(S);
CC:=BaseRing(tau);

Equ:=ReconstructCurveG4(tau);
quadric:=Equ[1];
cubic:=Equ[2];
CC4:=Parent(quadric);
f1 := new_cubic(quadric,cubic);
f2:=new_cubic(CC4!Q, CC4!F);

f1 := Evaluate(f1, [CC4.1*CC4.3  , CC4.2*CC4.3, CC4.1*CC4.4, CC4.2*CC4.4]);
f2 := Evaluate(f2, [CC4.1*CC4.3  , CC4.2*CC4.3, CC4.1*CC4.4, CC4.2*CC4.4]);

I1 := eval_inv(list_invariants, f1);
I2 := eval_inv(list_invariants, f2);
J1, J2 := same_wps(list_invariants, I1, I2);
for i in [1..#J1] do
        print Abs(J1[i]-J2[i]);
end for;
Time(T);

































