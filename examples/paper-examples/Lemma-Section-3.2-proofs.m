/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */
 

p:= 37;
bitangents := [
    [GF(p)| 1, 0, 0 ],
    [GF(p)| 0, 1, 0 ],
    [GF(p)| 0, 0, 1 ],
    [GF(p)| 1, 1, 1 ],
    [GF(p)| 35, 19, 1 ],
    [GF(p)| 27, 11, 1 ],
    [GF(p)| 34, 2, 1 ],
    [GF(p)| 1, 30, 20 ],
    [GF(p)| 23, 27, 30 ],
    [GF(p)| 12, 16, 23 ],
    [GF(p)| 13, 17, 23 ],
    [GF(p)| 24, 27, 31 ],
    [GF(p)| 1, 31, 21 ],
    [GF(p)| 18, 4, 28 ],
    [GF(p)| 11, 32, 18 ],
    [GF(p)| 12, 26, 36 ],
    [GF(p)| 7, 17, 24 ],
    [GF(p)| 14, 26, 16 ],
    [GF(p)| 13, 32, 33 ],
    [GF(p)| 10, 35, 23 ],
    [GF(p)| 34, 14, 23 ],
    [GF(p)| 32, 15, 23 ],
    [GF(p)| 3, 19, 10 ],
    [GF(p)| 3, 26, 19 ],
    [GF(p)| 30, 28, 7 ],
    [GF(p)| 8, 33, 30 ],
    [GF(p)| 10, 1, 24 ],
    [GF(p)| 14, 34, 23 ]
];
tritpairs := [ MatrixAlgebra(GF(37), 4) |
  Matrix(GF(37), 4, 4, \[ 25, 13, 15, 12, 13, 32, 8, 1, 15, 8, 28, 20, 12, 1, 
20, 33 ]),

  Matrix(GF(37), 4, 4, \[ 4, 20, 1, 28, 20, 0, 27, 20, 1, 27, 35, 26, 28, 20, 
26, 15 ]),

  Matrix(GF(37), 4, 4, \[ 0, 8, 23, 17, 8, 19, 9, 15, 23, 9, 19, 30, 17, 15, 
30, 19 ]),

  Matrix(GF(37), 4, 4, \[ 14, 25, 14, 23, 25, 7, 14, 11, 14, 14, 19, 10, 23, 
11, 10, 33 ]),

  Matrix(GF(37), 4, 4, \[ 8, 0, 15, 27, 0, 19, 34, 20, 15, 34, 23, 0, 27, 20, 
0, 17 ]),

  Matrix(GF(37), 4, 4, \[ 9, 7, 7, 28, 7, 21, 2, 3, 7, 2, 7, 4, 28, 3, 4, 26 
]),

  Matrix(GF(37), 4, 4, \[ 29, 10, 28, 32, 10, 35, 16, 31, 28, 16, 7, 32, 32, 
31, 32, 5 ]),

  Matrix(GF(37), 4, 4, \[ 16, 20, 25, 20, 20, 26, 7, 11, 25, 7, 12, 10, 20, 11,
10, 36 ]),

  Matrix(GF(37), 4, 4, \[ 9, 8, 28, 3, 8, 17, 6, 7, 28, 6, 23, 26, 3, 7, 26, 32
]),

  Matrix(GF(37), 4, 4, \[ 2, 23, 9, 3, 23, 8, 10, 21, 9, 10, 15, 30, 3, 21, 30,
33 ]),

  Matrix(GF(37), 4, 4, \[ 13, 17, 14, 0, 17, 14, 25, 11, 14, 25, 29, 8, 0, 11, 
8, 14 ]),

  Matrix(GF(37), 4, 4, \[ 29, 13, 15, 18, 13, 22, 17, 24, 15, 17, 12, 36, 18, 
24, 36, 23 ]),

  Matrix(GF(37), 4, 4, \[ 7, 17, 5, 33, 17, 14, 15, 8, 5, 15, 24, 35, 33, 8, 
35, 15 ]),

  Matrix(GF(37), 4, 4, \[ 35, 7, 10, 26, 7, 28, 34, 14, 10, 34, 28, 8, 26, 14, 
8, 20 ]),

  Matrix(GF(37), 4, 4, \[ 2, 31, 10, 22, 31, 18, 7, 8, 10, 7, 35, 18, 22, 8, 
18, 28 ]),

  Matrix(GF(37), 4, 4, \[ 15, 24, 0, 33, 24, 3, 6, 21, 0, 6, 4, 1, 33, 21, 1, 
34 ]),

  Matrix(GF(37), 4, 4, \[ 24, 11, 1, 16, 11, 22, 2, 32, 1, 2, 0, 34, 16, 32, 
34, 21 ]),

  Matrix(GF(37), 4, 4, \[ 30, 4, 24, 32, 4, 3, 18, 24, 24, 18, 10, 2, 32, 24, 
2, 24 ]),

  Matrix(GF(37), 4, 4, \[ 30, 13, 24, 26, 13, 33, 19, 30, 24, 19, 0, 22, 26, 
30, 22, 24 ]),

  Matrix(GF(37), 4, 4, \[ 30, 32, 3, 33, 32, 28, 2, 8, 3, 2, 25, 2, 33, 8, 2, 
13 ]),

  Matrix(GF(37), 4, 4, \[ 18, 14, 21, 32, 14, 19, 35, 18, 21, 35, 15, 23, 32, 
18, 23, 8 ]),

  Matrix(GF(37), 4, 4, \[ 14, 12, 16, 27, 12, 25, 22, 4, 16, 22, 19, 19, 27, 4,
19, 1 ]),

  Matrix(GF(37), 4, 4, \[ 21, 31, 10, 33, 31, 2, 16, 17, 10, 16, 16, 21, 33, 
17, 21, 36 ]),

  Matrix(GF(37), 4, 4, \[ 25, 9, 31, 14, 9, 25, 14, 20, 31, 14, 8, 17, 14, 20, 
17, 32 ]),

  Matrix(GF(37), 4, 4, \[ 18, 25, 23, 20, 25, 17, 4, 36, 23, 4, 4, 5, 20, 36, 
5, 14 ]),

  Matrix(GF(37), 4, 4, \[ 29, 18, 4, 21, 18, 20, 5, 26, 4, 5, 15, 17, 21, 26, 
17, 25 ]),

  Matrix(GF(37), 4, 4, \[ 27, 28, 33, 32, 28, 11, 34, 4, 33, 34, 36, 25, 32, 4,
25, 10 ]),

  Matrix(GF(37), 4, 4, \[ 22, 9, 35, 10, 9, 19, 34, 17, 35, 34, 23, 9, 10, 17, 
9, 23 ])
];

bitangents:=[bitangents[i]: i in [10, 23, 4, 20,  17, 9, 12,  1, 5, 11]];
tritpairs10 := [tritpairs[i]: i in [10, 23, 4, 20,  17, 9, 12,  1, 5, 11]];


p:= 37;R4 := PolynomialRing(GF(p^2), 4);
mats1new := tritpairs10;
CC4<x0,x1,x2,x3> := PolynomialRing(GF(p), 4);
  x:=Matrix(4,1,[CC4.i: i in [1..4]]);
  r:= #tritpairs10;
  mats1newx:=Matrix(CC4, 1, r,[(Transpose(x)*ChangeRing(mats1new[i],
CC4)*x)[1,1]: i in [1..r]]);
  Xnew:=Matrix([&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats1new]  );
vi:= Kernel(Xnew);
CC3 := PolynomialRing(GF(p), 3);
fs := [&+[el[i]*CC3.i : i in [1..3]] : el in bitangents[1..r]];
  mons:=MonomialsOfDegree(CC3,2);
  fsq_mat := [];
  for f in fs do
    cs := [];
    for m in mons do
      Append(~cs, MonomialCoefficient(f^2,m));
    end for;
    Append(~fsq_mat, cs);
  end for;
  fsq_mat := Matrix(fsq_mat);
fsq_mat;
si:= Kernel(fsq_mat);
sirows := Dimension(si);
N:=HorizontalJoin([  DiagonalMatrix(Eltseq(Basis(vi)[i]))*fsq_mat : i in
[1..Dimension(vi) ]]);
gammaiinv := Kernel(N);
basis := ExtendBasis(vi, RSpace(GF(p),10));
Upart := Matrix(basis[4..10]);
phi:= Upart*DiagonalMatrix(Eltseq(Basis(gammaiinv)[1]))*fsq_mat;
Qpre:=Kernel(phi);
  Qpre1:=Qpre*Upart;
  Qnew:=&+[Eltseq(Basis(Qpre1)[1])[i]*mats1new[i]: i in [1..r] ];
  dualelt:=mats1newx*ChangeRing(Transpose(Upart), CC4);
basis := ExtendBasis(Image(Transpose(phi)), RSpace(GF(p),7));
  phiext:=HorizontalJoin(phi, Matrix(7,1, basis[7]));
  phiTinv:=ChangeRing(Transpose(phiext)^(-1), CC4);
  print "phiTinv";
  print phiTinv;
  phiL:=dualelt*phiTinv;
  qdual:=ZeroMatrix(CC4, 3,3);
  count:=1;
  for i in [1..3] do
    for j in [i..3] do
      qdual[i,j]:=phiL[1,count];
        qdual[j, i]:=phiL[1, count];
        count+:=1;
    end for;
  end for;
  detqdual:=Determinant(qdual);
Q := (Transpose(x)*ChangeRing(Qnew, CC4)*x)[1,1];
I:= Ideal([Q, detqdual]);

print "bitangents:", bitangents;
print "\n tritangent pairs:", tritpairs10;
print "\n Linear system for the constants lambda_i", N;
print "\n lambda_i:", gammaiinv;
print "\n Basis for V_{C eta} is given by the following lin. comb. of generators ", Upart;
print "\n Matrix for phi", Transpose(phi);
print "\n Quadric", Qnew;

print "\n dualelt", dualelt;
print "\n extension of phi", phiext;
print "\n compostion in the definition of tilde{G}", phiL;
print "\n symmetric matrix in definition of tilde{G}", qdual;
print "\n Taking the discriminant", detqdual;
print "\n Recover the curve", Radical(I);


