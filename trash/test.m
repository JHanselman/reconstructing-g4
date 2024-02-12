AttachSpec("~/reconstructing-g4/magma/spec");
load "signs_gen2.m";

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
    return [[ 0 : i in [1..g] ], [ 0 : i in [1..g] ]];
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;



p := 37;
K := GF(p^2);
_<x> := PolynomialRing(K);
f := x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5);
C1 := HyperellipticCurve(f);
roots := Roots(f, K);
roots := [ro[1]: ro in roots];
Sort(~roots);
U:={1,3,5};
eta := EtaFunction0(2);

theta4 := [K!0: i in [1..16]];

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


theta2 := [Sqrt(th): th in theta4];
theta2 := correct_signs(theta2);
rosens := RosenhainInvariantsFromThetaSquares(theta2, 2: SpTrafo := false);
f2 := x*(x-1)*&*[x-ro: ro in rosens];
C2 := HyperellipticCurve(f2);

