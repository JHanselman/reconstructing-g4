AttachSpec("../magma/spec");
SetDebugOnError(true);
SetOutputFile("reconstruction-test-biell.txt");
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
F:=Evaluate(F, [x+y,z,w,x+y+z+w]);
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

end for;

