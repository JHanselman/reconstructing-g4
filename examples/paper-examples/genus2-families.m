function Dihedral_inv(f)
f:=Normalize(f);
K:=BaseRing(Parent(f));
coefff := Reverse(Coefficients(f)[1..5]);
R<t> :=PolynomialRing(K);
P<[x]> := PolynomialRing(Rationals(), 5);
Pa<xa>:=PolynomialRing(P);
Poly:=&*[&*[xa-x[i]-x[j]: i in [j+1..5]]: j in [1..4]];
M := SFAElementary(Rationals());
syms:=[M!c: c in Coefficients(Poly)];
ret:= R!0;
for nu in [1..#syms] do
s:= syms[nu];
supp, coeff := Support(s);
ind := [i: i in [1..#supp]| &and[j le 5: j in supp[i]] ];
supp:= [supp[i]: i in ind]; 
coeff:= [coeff[i]: i in ind]; 
ret +:= &+[coeff[i] * &*([(-1)^j*coefff[j]: j in supp[i]] cat [K!1]): i in [1..#ind]] * z^(nu-1);
end for;
return ret;
end function;


//Mestre's family with RM by Z[(1+sqrt(5))/2]:
K<u,t>:=RationalFunctionField(Rationals(),2);
R<z>:=PolynomialRing(K);
f:=(1-z)^3+u*z*((1-z)^3+u*z^2-z^3*(1-z))-t*z^2*(z-1)^2;
L:=quo<R|f>;          
RL:=PolynomialRing(L);
Factorization(RL!f);

//Family of curves with QM by B(2, 3) from Hashimoto & Murabayashi:
FF0<s>:=RationalFunctionField(Rationals());
R<x>:=PolynomialRing(FF0);
g:=4*s^2*x^2-s^2+x^2+2;
FF<t>:=ext<FF0| g>;
P:=-2*(s+t); R:= -2*(s-t); Q:=(1+2*t^2)*(11-28*t^2+8*t^4)/(3*(1-t^2)*(1-4*t^2));
S<X>:=PolynomialRing(FF);
fQ:= X*(X^4-P*X^3+Q*X^2-R*X+1);
fQtrans:=(X-t)^6*Evaluate(fQ, (X+t)/(X-t));
S0<X>:=PolynomialRing(FF0);
fQtrans:=S0!S!fQtrans;

