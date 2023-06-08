AttachSpec("../magma/spec");
SetDebugOnError(true);
SetOutputFile("reconstruction-test-invs.txt");
load "../../Invariants_genus_4.m";
ZZ:=Integers();
prec := 30;
SetDefaultRealFieldPrecision(prec);
QQ:=Rationals();
R<x,y,z,w> := PolynomialRing(QQ,4);

tritangentbasis := [
    [GF(2)|1, 1, 1, 0, 1, 1, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 1, 1, 0, 0, 0, 1, 0],
    [GF(2)|1, 0, 1, 0, 0, 1, 1, 0],
    [GF(2)|0, 1, 1, 0, 0, 1, 0, 0]];
mons2:=MonomialsOfDegree(R,2);
mons3:=MonomialsOfDegree(R,3);


for nu in [1..100] do
Q := &+[Random(-10,10)*m: m in mons2];
F := &+[Random(-10,10)*m: m in mons3];
C := Curve(Proj(R),[F,Q]);
print C;
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
if not IsIrreducible(p) then
	continue;
end if;
S := RiemannSurface(p : Precision := prec);
g := Genus(S);
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



end for;

