AttachSpec("../../magma/spec");
load "Galois.m";

function TakaseQuotientSq(theta4, eta, k, l, m)
g := #Eltseq(eta[1]) div 2;
U := {1,3,5};
Bm := { 1..2*g+1 }; L := [ bp : bp in (Bm diff { k, l, m }) ];
V := { L[i]: i in [1..g-1] }; W := { L[i]: i in [g..2*(g-1)] };
num1 := theta4[TCharToIndex(EtaValue(eta, U sdiff (V join { k, l })))];
num2 := theta4[TCharToIndex(EtaValue(eta, U sdiff (W join { k, l })))];
den1 := theta4[TCharToIndex( EtaValue(eta, U sdiff (V join { k, m })))];
den2 := theta4[TCharToIndex(EtaValue(eta, U sdiff (W join { k, m })))];
return (num1*num2)/(den1*den2);
end function;



function AreThetasRational(rosens)
g:=2;
L:=Parent(rosens[1]);
eta:=EtaFunction0(2);
U:={1,3,5};
roots:=[0,1,-1] cat rosens;
theta4 := [L!0: i in [1..16]];
for i in [1..5] do
        for j in [i+1..5] do
        	T := {i,j, 6};
		Tcomp := {1..6} diff T;
        	S := T sdiff U;
		assert(IsEven( #S));
        	sign := (-1)^(#(S meet U));
        	cha := EtaValue(eta, Setseq(S));
        	theta4[TCharToIndex(cha)] := sign * &*[ &*[(roots[nu]-roots[mu])^(-1): nu in T]: mu in Tcomp];
        end for;
end for;
print [theta4[i]/theta4[1]: i in [1..16]];
return &and[IsPower(theta4[i]/theta4[1], 4)  : i in [1..16]];


end function;


p:=101;

for S in Subsets({2..p-2},3) do
	if AreThetasRational([GF(p) |s:s in S]) then
		R<x>:=PolynomialRing(GF(p));
		H:=HyperellipticCurve(&*[x-s: s in S]*x*(x^2-1));
		J:=Jacobian(H);
		A:=AbelianGroup(J);
		if not &and[Valuation(e,2) ne 1: e in ElementaryDivisors(A)] then 
			print S;
			break;
		end if;
	end if;
end for;

 p := [19, 29, 61, 89];
 q := [5^3, 7^3];
