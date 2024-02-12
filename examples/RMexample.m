AttachSpec("~/reconstructing-g4/magma/spec");
load "Galois.m";
load "g2-rmlmfdb.m";




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
return ftest;
end if;
end for;

end function;



curves := make_data();
for loo0 in [1..2] cat [4..#curves] do
H1:= SimplifiedModel(curves[loo0]);
f:= HyperellipticPolynomials(H1);
x:= Parent(f).1;
if Degree(f) eq 5 then
if ConstantCoefficient(f) ne 0 then
f:= x*Reverse(f);
else 
continue;
end if;
end if;
Ig, Poi, U := IgusaTwist(f);
R4:=PolynomialRing(Rationals(),5);
Ig:=R4!(Ig/LeadingCoefficient(Ig));
TSpace := &+[Evaluate(Derivative(Ig, i), Poi) * R4.i: i in [1..5]];
P4:=ProjectiveSpace(R4);
Sch:=Scheme(P4, [Ig, TSpace]);
points := PointSearch(Sch, 200: Nonsingular:=false);
L := BaseRing(U);
for loo in [1..#points] do
	Poi2:=Eltseq(points[loo]);
	TSpace2 := &+[Evaluate(Derivative(Ig, i), Poi2) * R4.i: i in [1..5]];
	if Evaluate(TSpace2, Poi) eq 0 then
		continue;
	end if;
	theta41part := Eltseq(Vector(L, Poi2)*U);
	vecs := IgusaCoordinates();
	theta41 := [&+[vec[i]* theta41part[i]: i in [1..5]]: vec in vecs];
	if #[th: th in theta41| th eq 0] gt 6 then
		continue;
	end if;
	f2 := HyperellipticCurveFromTheta4(theta41);
	pl:=InfinitePlaces(L)[1];
	S<X,Y>:=PolynomialRing(L,2);
	equ := Y^2- Evaluate(f2, X);
	Surf:= RiemannSurface(equ, pl: Precision:= 80);
	Pi := BigPeriodMatrix(Surf);
	EndJ2:= EndomorphismAlgebra(Pi);
	print EndJ2, loo0, loo;
	print "\n";
end for;
end for;


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
