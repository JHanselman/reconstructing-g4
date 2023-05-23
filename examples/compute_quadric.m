AttachSpec("/home/sijsling/g4rec/gluing/magma/spec");
SetDebugOnError(true);
//prec := 40;
prec := 15;
SetDefaultRealFieldPrecision(prec);
RR := RealField(prec);
CC<I> := ComplexField(prec);
R<x,y,z,w> := PolynomialRing(QQ,4);
Q := x^2 + 2*y^2 + z^2 - w^2;
F := 3*x^3 + y^3 - 7*z^3 + w^3;
C := Curve(Proj(R),[F,Q]);
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
assert IsIrreducible(p);
S := RiemannSurface(p : Precision := prec);
g := Genus(S);
Pi_big := BigPeriodMatrix(S);
Pi1, Pi2 := SplitBigPeriodMatrix(Pi_big);
tau := SmallPeriodMatrix(S);
chars := OddThetaCharacteristics(g);

// make Steiner complex
v := Vector([GF(2) | 0, 0, 0, 0, 1, 0, 0, 0 ]);
steiner := [];
steinerrie:=[];
//First seven are an Aronhold set. 8th one is l2 + l3
aronhold:=[[[1,1,1],[1,1,1]],[[0,0,1],[0,1,1]],[[0,1,1],[0,0,1]],[[1,0,1],[1,0,0]],[[1,0,0],[1,0,1]],[[1,1,0],[0,1,0]],[[0,1,0],[1,1,0]], [[0,1,0],[0,1,0]]];
rie:=[[[0,0,0],[0,0,1]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,1]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]], [[0,0,1],[1,1,0]],[[0,0,0],[1,1,0]],[[0,0,0],[0,0,0]],[[0,0,0],[0,1,0]],   [[0,0,0],[1,0,0]]];

for aro in aronhold do

    cc1 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!0] cat [GF(2)!el: el in  aro[2]];
    cc2 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!1] cat [GF(2)!el: el in  aro[2]];
    c1:=[[0] cat aro[1] , [0] cat aro[2]];
    c2:=[[0] cat aro[1] , [1] cat aro[2]];

    Append(~steiner, [c1,c2]);


end for;

for cha in rie do

    cc1 := [GF(2)!0] cat  [GF(2)!el : el in cha[1]] cat [GF(2)!0] cat [GF(2)!el: el in  cha[2]];
    cc2 := [GF(2)!0] cat  [GF(2)!el : el in cha[1]] cat [GF(2)!1] cat [GF(2)!el: el in  cha[2]];
    c1:=[Vector(QQ, [0] cat cha[1]) , Vector(QQ, [0] cat cha[2])];
    c2:=[Vector(QQ, [0] cat cha[1]) , Vector(QQ, [1] cat cha[2])];
    Append(~steinerrie,  c1);
    Append(~steinerrie,  c2);
end for;

tritangents := [];
//for char in chars[1..4] do
//for char in chars do
for s in steiner[1..8] do
  new := [];
  for char in s do
      print(char);
    Append(~new, TritangentPlane(Pi_big, char));
    //Append(~new, TritangentPlane(HorizontalJoin(IdentityMatrix(CC, 4), tau), char));
  end for;
  Append(~tritangents, new);
end for;


Pij, thetas := ThetaBatch(tau, steinerrie, 3);
sqrt_thetas := [Sqrt(th) : th in thetas];


// d/dPij sqrt(theta_1)sqrt(theta_9) = 1/2(dtheta1/dPij sqrt(theta9))/sqrt(theta1) + 1/2(dtheta9/dPij sqrt(theta1))/sqrt(theta9)

r := [sqrt_thetas[1+2*t] * sqrt_thetas[2*t + 2] : t in [0..11]];



Q:= ZeroMatrix(CC, 4);
errors := [];
signlist := [[1,1,1], [1, 1,-1], [1,-1,1], [1,-1,-1]];
for s in signlist do
    sum := Abs(&+[s[i+1] * &*[r[4*i+j] : j in [1..4]] : i in [0..2]]);
    Append(~errors, sum);
end for; 

min, i:= Minimum(errors);
assert min lt 10^(-prec/2);

signs := signlist[i];
for i in [1..4] do
    for j in [i..4] do
       
       Q[i,j] := &+[ signs[t+1] * &+[  (Pij[2*s -1 +8*t][i,j]*sqrt_thetas[2*s +8*t]/sqrt_thetas[2*s -1 + 8*t] +  Pij[2*s +8*t][i,j] * sqrt_thetas[2*s - 1 +8*t]/sqrt_thetas[2*s + 8*t])* &*[r[u + 4*t] : u in [1..s-1] cat [s+1..4]]  :s in [1..4]] : t in [0..2] ]; 
    end for; 
end for;

Q_sym := (Q + Transpose(Q))/2;
Quadric := Transpose(Pi1^(-1)) * Q_sym * Pi1^(-1);
Q_list := Eltseq(Quadric);
_, i := Maximum([Abs(l) : l in Q_list]);
factor := Q_list[i];
Quadric := Quadric/factor;

//Compare with
DefiningEquations(CanonicalImage(S))[3];

