AttachSpec("../magma/spec");
SetDebugOnError(true);
SetVerbose("User1", 1);
ZZ:=Integers();
prec := 40;
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


Q:= 4*x^2 + 6*x*z + 8*x*w + y^2 + y*z + 3*y*w - 10*z^2 - 9*z*w - 10*w^2;
F:=5*x^3 - 3*x^2*y - 8*x^2*z + 3*x^2*w + x*y^2 + 4*x*y*z - 6*x*y*w + 4*x*z^2 +
    3*x*z*w + 7*x*w^2 - y^3 - 6*y^2*z + 10*y^2*w + 8*y*z^2 - 4*y*z*w + 3*y*w^2 -
    9*z^3 - 10*z^2*w - z*w^2 + 6*w^3;


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
CC1<t>:=PolynomialRing(CC);
X:=Matrix(CC4, 4,1, [CC4.i: i in [1..4]]);

TTB:=[];

for c in tritangentbasis do
  chara := [ZZ!v : v in Eltseq(c)];
  chara := [chara[1..4], chara[5..8]];
  Append(~TTB, TritangentPlane(Pi_big, chara));
end for;


TtoS := Matrix(TTB[1..4]);
D := DiagonalMatrix(Eltseq(Vector(TTB[5]) * (TtoS)^-1));
M := TtoS^-1 * D^-1;


Ccan, map:=CanonicalImage(S);
Cplane:=Domain(map);

for i in [-10..-1] cat [1..10] do
        f1:=Evaluate(DefiningEquation(Cplane), [1/CC!i+CC.1, t]);
        ys:=[roo[1]: roo in Roots(f1)];
        coord:=[[Evaluate(DefiningEquations(map)[nu], [1/CC!i+CC.1, y]) : nu in [1..4]  ]: y in ys];
        print [Abs(Evaluate( (Transpose(X)*ChangeRing(quadric, CC4) *X)[1,1], Eltseq(Vector(coo)*Transpose(Inverse(M))) )): coo in coord];
        print [Abs(Evaluate( (Transpose(X)*ChangeRing(cubic, CC4) *X)[1,1], Eltseq(Vector(coo)*Transpose(Inverse(M))) )): coo in coord];
end for;
