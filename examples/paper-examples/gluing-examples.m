// make curves to glue using formulas in Howe--Leprevost--Poonen
QQ := RationalsExtra(20);
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

S<s> := PolynomialRing(QQ);

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
assert &and[el[2] eq 1 : el in Factorization(hev)];
//Hev := HyperellipticCurve(hev);
//
hev1 := [Evaluate(el, 1) : el in Coefficients(hev)];
hev1 := S!hev1;
H1 := HyperellipticCurve(hev1);

hev2 := [Evaluate(el, 2) : el in Coefficients(hev)];
hev2 := S!hev2;
H2 := HyperellipticCurve(hev2);

// now glue
AttachSpec("~/github/gluing/magma/spec");
AllArithmetic2GluingsCCFor22(H1, H2, BaseField(H1));

