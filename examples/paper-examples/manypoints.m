AttachSpec("../../magma/spec");
SetDebugOnError(true);
//load "Galois.m";
//load "gluingfuncs.m";
p := 19;

R<x>:=PolynomialRing(GF(p));
//f := 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;
//f :=(x^2-2)*(x-3)*(x-4)*(x-5)*(x-6);
f:= 15*x^6+4*x^5+4*x^3+9*x^2+x+9; 
/*
if Degree(f) ne 6 then
  f := x*Reverse(f);
end if;
*/

Ig, Poi, U0 := IgusaTwist(f);
R4:=PolynomialRing(GF(p),5);
Ig:=R4!Ig;
TSpace := &+[Evaluate(Derivative(Ig, i), Poi) * R4.i: i in [1..5]];
P4:=ProjectiveSpace(R4);
Sch:=Scheme(P4, [Ig, TSpace]);
print "looking for points on Igusa quartic";
points := Points(Sch);
L := BaseRing(U0);
print "looping over points found";
cur_min := 100;
cur_max := 0;
for loo in [1..#points] do
  Poi2:=Eltseq(points[loo]);
  TSpace2 := &+[Evaluate(Derivative(Ig, i), Poi2) * R4.i: i in [1..5]];
  if Evaluate(TSpace2, Poi) eq 0 then
    continue;
  end if;
  theta41part := Eltseq(Vector(L, Poi2)*U0);
  vecs, equfin, xi := IgusaCoordinates(GF(p));
  evs := [Evaluate(x, theta41part) : x in xi];
  theta41 := [&+[vec[i]*theta41part[i]: i in [1..5]]: vec in vecs];
  if #[th: th in theta41| th eq 0] gt 6 then
      continue;
  end if;
  f2, ros := HyperellipticCurveFromTheta4(theta41);
  S<X,Y>:=PolynomialRing(L,3);
  equ := Y^2- Evaluate(f2, X);


  /*
  Surf:= RiemannSurface(equ, pl: Precision:= 80);
  Pi := BigPeriodMatrix(Surf);
  EndJ2:= EndomorphismAlgebra(Pi);
  print EndJ2, loo;
  print "\n";
  */

  R:= Parent(f2);      
  IgusaInvariants(f2, R!0: normalize:=true);
  S, W:=IgusaInvariants(f2, R!0: normalize:=true);
  S:=[GF(p)!s : s in S];

  X1 := HyperellipticCurve(f);
  X2 := HyperellipticCurveFromIgusaInvariants(S);
  f2Q, g := HyperellipticPolynomials(X2);
  f2Q := EvenModel(f2Q);
  X2 := HyperellipticCurve(f2Q);
  len := #X2;
  if len lt cur_min then
       cur_min := len;
  elif len gt cur_max then
       cur_max := len;
  end if;
  if len in [33, 7] then
      print X2, Factorization(LPolynomial(X2));
  end if;
  /*
  if Degree(f2Q) ne 6 then
    f2Q := x*Reverse(f2Q);
  end if;
  */

  /*
  bool, c := IsQuadraticTwist(BaseChange(X2, L), HyperellipticCurve(f2));
  twistX2:= HyperellipticCurve(c*f2);
  bool2, phi := IsIsomorphic(twistX2 , BaseChange(X2, L));

  roots_f2 := [phi(twistX2![r, 0])[1] : r in ([L!0,1] cat ros)] cat [phi(twistX2![1,0,0])[1]];

  RS1 := RiemannSurface(f, 2: Precision:= 20); RS2 := RiemannSurface(f2Q, 2: Precision:= 20);

  V:= findV(RS1, RS2, roots_f2 );

  P1 := BigPeriodMatrix(RS1); P2 := BigPeriodMatrix(RS2);
  P := DiagonalJoin(P1, P2);
  Q := QFromPVFor22(P, V);
  */
end for;
print cur_max, cur_min;

/*
prec := 300;
F := RationalsExtra(prec);
P3<x,y,z,w> := ProjectiveSpace(F, 3);

f1 := x^2 + 28/5*x*y + 18/5*y^2 - 224/5*z^2 + 102/5*z*w - 33/10*w^2;
f2 := -27/226*x*y^2 - 205/1808*x*z^2 + 15/113*x*z*w - 185/3616*x*w^2 - 39/452*y^3 + y*z^2 - 271/452*y*z*w + 17/452*y*w^2;
X := Curve(P3, [f1, f2]);
*/

/*
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

loo = 67, 121 are wrong

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

