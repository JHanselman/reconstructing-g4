R<x0,x1,x2>:=PolynomialRing(Rationals(),3);
S<t0,t1,t2,t3>:=PolynomialRing(Rationals(),4);

//Cayley Cubic
A0:=Matrix(R, 3,3,[1,0,0,0,0,0,0,0,0]);
A1:=Matrix(R,3,3, [0,0,0,0,1,0,0,0,0]);
A2:=Matrix(R,3,3, [0,0,0,0,0,0,0,0,1]);
A3:=Matrix(R,3,3, [1,1,1,1,1,1,1,1,1]);
As:=[A0,A1,A2,A3];
x:=Matrix(R, 3,1, [x0,x1,x2]);
ci:=[(-1)^(i-1)*Determinant(HorizontalJoin([As[j]*x :  j in [1..i-1] cat [i+1..4]]))  : i in [1..4] ];
Gamma:=t1*t2*t3+t0*t2*t3+t0*t1*t3+t0*t1*t2;
N:=x0*x1*x2*(x0+x1+x2);
Aij:=&cat[[ ci[i]*ci[j] div N: j in [i+1..4]]: i in [1..4]];
tij:=&cat[[ S.i*S.j : j in [i+1..4]]: i in [1..4]];
mons:= MonomialsOfDegree(R,2);
cs:=[];
for m in mons do
	Append(~cs, [MonomialCoefficient(Aij[nu], m): nu in [1..6]]);
end for;
T:=Matrix(cs);
v:=  ChangeRing(Transpose(T)^(-1),S)*Matrix(S, 6,1,tij);
X:=ZeroMatrix(S, 3,3);
nu:=1;
for i in [1..3] do
	for j in [i..3] do
		X[i,j]:=v[nu,1];
		if i eq j then
			X[i,j] /:=2;
		end if;
		nu+:=1;
	end for;
end for;
X:=X+Transpose(X);
//Cofactor matrix recovers the symmetrization.
Cof:=Matrix([[Cofactor(X, i, j) div Gamma: i in [1..3]]: j in [1..3]]);
print Cof, As;



//Cubic with 2 A1 and one  A3
A0:=Matrix(R, 3,3,[1,0,0,0,0,0,0,0,0]);
A1:=Matrix(R,3,3, [0,0,0,0,1,0,0,0,0]);
A2:=Matrix(R,3,3, [0,0,1,0,0,1,1,1,0]);
A3:=Matrix(R,3,3, [0,0,0,0,0,0,0,0,1]);
As:=[A0,A1,A2,A3];
x:=Matrix(R, 3,1, [x0,x1,x2]);
ci:=[(-1)^(i-1)*Determinant(HorizontalJoin([As[j]*x :  j in [1..i-1] cat [i+1..4]]))  : i in [1..4] ];
tij:=[S.1*S.2,S.1*S.3,S.2*S.3,   S.3^2, S.1*S.4, S.2*S.4 ];
Aij:=[Evaluate(el, ci): el in tij];
N:=GCD(Aij);
Aij:=[el div N: el in Aij  ];

mons:= MonomialsOfDegree(R,2);
cs:=[];
for m in mons do
        Append(~cs, [MonomialCoefficient(Aij[nu], m): nu in [1..6]]);
end for;
T:=Matrix(cs);
v:=  ChangeRing(Transpose(T)^(-1),S)*Matrix(S, 6,1,tij);
X:=ZeroMatrix(S, 3,3);
nu:=1;
for i in [1..3] do
        for j in [i..3] do
                X[i,j]:=v[nu,1];
                if i eq j then
                        X[i,j] /:=2;
                end if;
                nu+:=1;
        end for;
end for;
X:=X+Transpose(X);
print Factorization(Determinant(X));
//Cofactor matrix recovers the symmetrization.
//Cof:=Matrix([[Cofactor(X, i, j) div Gamma: i in [1..3]]: j in [1..3]]);
//print Cof, As;



//Cubic with one A1 and one  A5
A0:=Matrix(R, 3,3,[1,0,0,0,0,0,0,0,0]);
A1:=Matrix(R,3,3, [0,0,1,0,1,0,1,0,0]);
A2:=Matrix(R,3,3, [0,0,0,0,0,1,0,1,0]);
A3:=Matrix(R,3,3, [0,0,0,0,0,0,0,0,1]);
As:=[A0,A1,A2,A3];
x:=Matrix(R, 3,1, [x0,x1,x2]);
ci:=[(-1)^(i-1)*Determinant(HorizontalJoin([As[j]*x :  j in [1..i-1] cat [i+1..4]]))  : i in [1..4] ];
//TODO: Find the right tij
tij:=[S.1*S.2,S.1*S.3,S.2^2,S.2*S.4-S.3^2, S.2*S.3, S.1*S.4 ];
Aij:=[Evaluate(el, ci): el in tij];
N:=GCD(Aij);
Aij:=[el div N: el in Aij  ];
mons:= MonomialsOfDegree(R,2);
cs:=[];
for m in mons do
        Append(~cs, [MonomialCoefficient(Aij[nu], m): nu in [1..6]]);
end for;
T:=Matrix(cs);
v:=  ChangeRing(Transpose(T)^(-1),S)*Matrix(S, 6,1,tij);
X:=ZeroMatrix(S, 3,3);
nu:=1;
for i in [1..3] do    
        for j in [i..3] do
                X[i,j]:=v[nu,1];
                if i eq j then
                        X[i,j] /:=2;
                end if;
                nu+:=1;
        end for;
end for;
X:=X+Transpose(X);
print Factorization(Determinant(X));
//Cofactor matrix recovers the symmetrization.
//Cof:=Matrix([[Cofactor(X, i, j) div Gamma: i in [1..3]]: j in [1..3]]);
//print Cof, As;

