AttachSpec("../../magma/spec");
load "Galois.m";
load "gluingfuncs.m";
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



function HyperellipticCurveFromTheta4(theta4)
g:=2;
L:=Parent(theta4[1]);
eta:=EtaFunction0(2);
ros := [ Sqrt(TakaseQuotientSq(theta4, eta, 1, l, 2)) : l in [3..2*g+1] ];
U:={1,3,5};
R<X>:=PolynomialRing(L);
for loo in [0..7] do
bin:=Intseq(loo,2,3);
rostest:= [ros[water]*(-1)^bin[water]: water in [1..3]];
rootstest := [0,1] cat [-ro/(ro-2): ro in rostest] cat [-1];
ftest:=X*(X-1)*&*[(X-L!ro): ro in rostest];
//print IgusaInvariantsEqual(IgusaInvariants(f), IgusaInvariants(ftest));
theta4test:= [L!0: i in [1..16]];
for i in [1..5] do
        for j in [i+1..5] do
        T := {i,j, 6};
        Tcomp := {1..6} diff T;
        S := T sdiff U;
        sign := (-1)^#(S meet U);
        cha := EtaValue(eta, Setseq(S));
        theta4test[TCharToIndex(cha)] := sign * &*[ &*[(rootstest[nu]-rootstest[mu])^(-1): nu in T]: mu in Tcomp];
        end for;
end for;
if &and[theta4test[i]*theta4[1] eq theta4[i]*theta4test[1]: i in [1..16]] then
return ftest, rostest;
end if;
end for;

end function;





R<x>:=PolynomialRing(Rationals());
f := 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;

if Degree(f) ne 6 then
  f := x*Reverse(f);
end if;

Ig, Poi, U0 := IgusaTwist(f);
R4:=PolynomialRing(Rationals(),5);
Ig:=R4!Ig;
TSpace := &+[Evaluate(Derivative(Ig, i), Poi) * R4.i: i in [1..5]];
P4:=ProjectiveSpace(R4);
Sch:=Scheme(P4, [Ig, TSpace]);
//points := PointSearch(Sch, 500);
L := BaseRing(U0);
for loo in [1..#points] do
loo:= 142;
Poi2:=Eltseq(points[loo]);
	TSpace2 := &+[Evaluate(Derivative(Ig, i), Poi2) * R4.i: i in [1..5]];
	if Evaluate(TSpace2, Poi) eq 0 then
		continue;
	end if;
	theta41part := Eltseq(Vector(L, Poi2)*U0);
	vecs, equfin, xi := IgusaCoordinates();
    evs := [Evaluate(x, theta41part) : x in xi];
    if #Seqset(evs) ne #evs then
      print "Found multiple gluings!";
      printf "point: %o\n", points[loo];
    end if;
	theta41 := [&+[vec[i]*theta41part[i]: i in [1..5]]: vec in vecs];
	if #[th: th in theta41| th eq 0] gt 6 then
		continue;
	end if;
	f2, ros := HyperellipticCurveFromTheta4(theta41);
	pl:=InfinitePlaces(L)[1];
	S<X,Y>:=PolynomialRing(L,2);
	equ := Y^2- Evaluate(f2, X);
    /*
	Surf:= RiemannSurface(equ, pl: Precision:= 80);
	Pi := BigPeriodMatrix(Surf);
	EndJ2:= EndomorphismAlgebra(Pi);
	print EndJ2, loo;
	print "\n";
    */
end for;

R:= Parent(f2);      
IgusaInvariants(f2, R!0: normalize:=true);
S, W:=IgusaInvariants(f2, R!0: normalize:=true);
S:=[Rationals()!s : s in S];

X1 := HyperellipticCurve(f);
X2 := HyperellipticCurveFromIgusaInvariants(S);
f2Q, g := HyperellipticPolynomials(X2);

if Degree(f2Q) ne 6 then
  f2Q := x*Reverse(f2Q);
  X2 := HyperellipticCurve(f2Q);
end if;



bool, c := IsQuadraticTwist(BaseChange(X2, L), HyperellipticCurve(f2));
twistX2:= HyperellipticCurve(c*f2);
bool2, phi := IsIsomorphic(twistX2 , BaseChange(X2, L));

roots_f2 := [phi(twistX2![r, 0])[1] : r in ([L!0,1] cat ros)] cat [phi(twistX2![1,0,0])[1]];

RS1 := RiemannSurface(f, 2: Precision:= 20); RS2 := RiemannSurface(f2Q, 2: Precision:= 20);

V:= findV(RS1, RS2, roots_f2 );

P1 := BigPeriodMatrix(RS1); P2 := BigPeriodMatrix(RS2);
Vs := AllVs2For22();
P := DiagonalJoin(P1, P2);
Qs := [];
for V in Vs do
  Q := QFromPVFor22(P, V);
  Append(~Qs, Q);
end for;

good := [];
for i->Q in Qs do
  print i;
  w1, w2 := SplitBigPeriodMatrix(Q);
  tau := w1^(-1) *w2;
  err := Abs(SchottkyModularForm(tau : prec := 20));
  print err;
  if err lt 10^-10 then
    Append(~good, i);
  end if;
end for;


V := Vs[120];

Q := QFromPVFor22(P, V);
//y^2 = 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;
//y^2 = 3*x^5 - 68*x^4 + 159*x^3 + 232*x^2 - 132*x + 16
//Output:
[
    x^2 + 28/5*x*y + 18/5*y^2 - 224/5*z^2 + 102/5*z*w - 33/10*w^2,
    -27/226*x*y^2 - 205/1808*x*z^2 + 15/113*x*z*w - 185/3616*x*w^2 - 39/452*y^3 
        + y*z^2 - 271/452*y*z*w + 17/452*y*w^2
]


//tau := SmallPeriodMatrix(Q);

CC:= Parent(Q[1,1]);

//Change polarization convention
//Q := Q* DiagonalJoin(IdentityMatrix(CC, 4), -IdentityMatrix(CC, 4)) ;

eqs := RationalReconstructCurveG4(Q);  
eqs;

P3 := ProjectiveSpace(Parent(eqs[1]));
C := Curve(Scheme(P3, eqs));

p:= 19;
Factorisation(LPolynomial(Curve(Reduction(C, p))));
Factorisation(LPolynomial(Curve(Reduction(X1, p))));
Factorisation(LPolynomial(Curve(Reduction(X2, p))));


P2<x, y> := PolynomialRing(Rationals(), 2);
plane_eq := Evaluate(Resultant(eqs[1], eqs[2], P3.1), [0,x,y,1]);
result_C:=Curve(Scheme(AffineSpace(P2),plane_eq));
result_P := BigPeriodMatrix(RiemannSurface(plane_eq:Precision:=500));
EndomorphismAlgebra(result_P);

/*
Found multiple gluings!
point: (-11/28 : -11/24 : 125/336 : -1/336 : 1)
Found multiple gluings!
point: (29/98 : -157/98 : 40/49 : 3/392 : 1)
Found multiple gluings!
point: (-1/28 : -29/56 : 5/112 : 3/112 : 1)

Poi2 := [-11/28, -11/24, 125/336, -1/336,1];  

Result:

238
1.99692468311987660851692368119E-28


552
6.94746444747365946288892131774E-30


*/

/*
[ (-1/4 : -9/20 : -9/20 : 1/40 : 1), (-17/20 : 3/20 : 3/80 : -1/80 : 1), (-9/40 
: -3/4 : 9/20 : 0 : 1), (3/46 : -7/23 : 5/92 : 7/184 : 1), (-3/14 : 9/14 : 
-15/28 : 3/56 : 1), (17/28 : 5/28 : 29/112 : 9/112 : 1), (-13/10 : 27/10 : -9/20
: 1/40 : 1), (61/144 : -3/4 : -5/48 : 1/18 : 1), (5/112 : -9/28 : 9/112 : 1/28 :
1), (1/20 : -13/40 : 1/80 : 3/80 : 1), (-7/16 : 17/12 : 5/48 : 1/24 : 1), (1/9 :
0 : -7/24 : 1/18 : 1), (7/20 : -3/8 : -21/80 : 1/16 : 1), (-5/68 : 21/68 : -6/17
: 7/136 : 1), (5/12 : -31/36 : 7/18 : 1/24 : 1), (53/28 : -75/28 : 6/7 : 5/56 : 
1), (1/76 : -25/76 : 2/19 : 5/152 : 1), (-1/128 : -21/64 : 15/128 : 1/32 : 1), 
(7/48 : -1/2 : -1/16 : 1/24 : 1), (7/30 : -1/2 : 1/5 : 1/24 : 1), (-7/60 : 1/4 :
-1/10 : 1/24 : 1), (13/8 : 15/16 : 21/32 : 5/32 : 1), (-19/4 : 21/8 : 3/32 : 
-7/32 : 1), (-17/8 : 21/16 : 21/64 : -5/64 : 1), (-9/16 : 0 : 3/16 : 0 : 1), 
(-129/80 : 17/20 : 11/80 : -1/20 : 1), (5 : -6 : -3/8 : 1/4 : 1), (89/48 : -11/4
: 15/16 : 1/12 : 1), (-1/2 : -1/2 : -1/8 : 0 : 1), (5/42 : -2/7 : 1/14 : 1/24 : 
1), (-35/26 : 33/52 : 21/104 : -1/26 : 1), (5/16 : -13/16 : -7/64 : 3/64 : 1), 
(-1/28 : -29/56 : 5/112 : 3/112 : 1), (1/2 : -3/4 : -3/16 : 1/16 : 1), (5/2 : -3
: 3/4 : 1/8 : 1), (5/18 : -1/12 : 1/8 : 1/18 : 1), (1/10 : -1/3 : -1/30 : 1/24 :
1), (148/11 : 108/11 : 57/11 : 1 : 0), (-1/6 : 1/2 : 0 : 1/24 : 1), (-1/4 : 
21/40 : -51/160 : 7/160 : 1), (7/48 : -1/3 : 5/48 : 1/24 : 1), (83/16 : -21/4 : 
15/16 : 1/4 : 1), (-5/2 : -41/20 : -37/40 : -3/20 : 1), (1/12 : -1/4 : 0 : 1/24 
: 1), (-5/48 : 1/4 : -1/16 : 1/24 : 1), (39/20 : -63/20 : 6/5 : 3/40 : 1), 
(-7/60 : 19/60 : -1/30 : 1/24 : 1), (-15/44 : -7/44 : 2/11 : 1/88 : 1), (27/2 : 
-25/2 : 7/4 : 5/8 : 1), (-13/22 : -5/44 : -1/88 : 0 : 1), (20 : -18 : -3 : 1 : 
0), (1/14 : -15/28 : -3/56 : 1/28 : 1), (1 : -2 : 1 : 0 : 0), (-1/10 : -3/10 : 
3/20 : 1/40 : 1), (-3/2 : 9/4 : -3/8 : 0 : 1), (-13/4 : 9/4 : 0 : -1/8 : 1), 
(5/12 : -17/12 : -1/6 : 1/24 : 1), (1/36 : -5/12 : -1/48 : 5/144 : 1), (89/4 : 
-39/2 : 9/4 : 1 : 0), (-15/104 : 33/52 : -27/52 : 3/52 : 1), (23/72 : -5/12 : 
-1/12 : 1/18 : 1), (-3 : 10 : 1 : 0 : 0), (-13/14 : 2/7 : 5/28 : -1/56 : 1), 
(15/16 : -51/32 : 63/64 : 3/64 : 1), (-1/6 : 1/3 : -1/6 : 1/24 : 1), (-9/10 : 
21/20 : 9/40 : 0 : 1), (-3/28 : 23/84 : -1/21 : 1/24 : 1), (-3/5 : 0 : 3/40 : 0 
: 1), (7/30 : -4/5 : -1/10 : 1/24 : 1), (13/176 : -3/22 : -15/176 : 1/22 : 1), 
(1/8 : -9/40 : 9/160 : 7/160 : 1), (1/26 : 0 : -9/52 : 5/104 : 1), (-1/80 : 3/20
: -21/80 : 1/20 : 1), (5/42 : -17/42 : -1/21 : 1/24 : 1) ]
*/

