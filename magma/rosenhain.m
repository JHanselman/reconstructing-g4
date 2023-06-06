function IsEvenVector(v)
w := Eltseq(v); g := #w div 2;
w1 := w[1..g]; w2 := w[(g + 1)..(2*g)];
return IsZero((Integers() ! ( &+[ w1[i]*w2[i] : i in [1..g] ])) mod 2);
end function;


function ComputeGamma(v0s)
M:=Matrix(GF(2), [[1, 0, 0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 1, 0, 0, 0],
[0, 1, 0, 0, 1, 0, 0, 0],
[0, 1, 0, 0, 1, 1, 0, 0],
[0, 0, 1, 0, 1, 1, 0, 0],
[0, 0, 1, 0, 1, 1, 1, 0],
[0, 0, 0, 1, 1, 1, 1, 0],
[0, 0 ,0 ,1 ,1 ,1 ,1, 1]]);
N:=Matrix( [Vector(v0s[i][1] cat v0s[i][2])-Vector(v0s[10][1] cat v0s[10][2]) : i in [1..8]]);
zer:=ZeroMatrix(GF(2), 4,4);
id:=IdentityMatrix(GF(2),4);
J:=BlockMatrix([[zer, id],[id,zer]]);
print N*J*Transpose(N);
print M*J*Transpose(M);


return Transpose(ChangeRing( M^(-1)*(ChangeRing(N, GF(2))), Integers()));
end function;


/* Mumford's eta */
function EtaFunction0();
eta1 := [1,0,0,0,0,0,0,0]; eta2 := [1,0,0,0,1,0,0,0];
eta3 := [0,1,0,0,1,0,0,0]; eta4 := [0,1,0,0,1,1,0,0];
eta5 := [0,0,1,0,1,1,0,0]; eta6 := [0,0,1,0,1,1,1,0];
eta7 := [0,0,0,1,1,1,1,0]; eta8:= [0,0,0, 1,1,1,1,1];
eta9 := [0,0,0,0,1,1,1,1]; eta10 := [0,0,0,0,0,0,0,0];
etas := [eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10];
return [ Transpose(Matrix(Integers(), [eta])) : eta in etas ];
end function;


/* Transform Mumford's eta by gamma */
/* Note that these function are represented by tuples, but we give an evaluation later */
function EtaFunction(gamma)
return [ gamma *v : v in EtaFunction0() ];
end function;


/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Integers(), [[ 0 : i in [1..2*g] ]]));
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..4], res[5..8]];
end function;


/* Given eta, finds U in Definition 3.3 */
function UFromEtaFunction(eta)
g := #Eltseq(eta[1]) div 2;
return { i : i in [1..#eta] | not IsEvenVector(eta[i]) } join { 2*g + 2 };
end function;



/* This finds the unique even zero characteristic */
function FindDelta(thetas : prec := 20)
if prec eq 0 then
    //CC := Parent(thetas[1]); prec := Precision(CC) - 100;
    CC := Parent(thetas[1]); prec := Precision(CC) / 2;
end if;
even:=EvenThetaCharacteristics(4);
v0s:=[];
for cha in even do
    theta := thetas[TCharToIndex(cha)];
    test := (Abs(theta) lt 10^(-prec));
    if test then
        Append(~v0s, cha);
    end if;
end for;
return v0s;
end function;


function EtaInnerproduct(eta1, eta2)
eta1 := [ Integers() ! c : c in Eltseq(eta1) ];
eta2 := [ Integers() ! c : c in Eltseq(eta2) ];
g := #eta1 div 2;
return &+[ eta1[i]*eta2[g + i] : i in [1..g] ];
end function;


function LogEpsilon(U, j)
if not j in U then
    return 1;
end if;
return 0;
end function;


/* Calculates sign */
function EpsilonKLM(eta, k, l, m)
U := UFromEtaFunction(eta);
exp := EtaInnerproduct(eta[k], (eta[k] + eta[l] + eta[m])) + LogEpsilon(U, k) - 1;
return (-1)^exp;
end function;


/* Theorem 4.5 */
function TakaseQuotient(thetas_sq, eta, k, l, m)
U := UFromEtaFunction(eta);
Bm := { 1..9 }; L := [ bp : bp in (Bm diff { k, l, m }) ];
V := { L[1], L[2], L[3] }; W := { L[4], L[5], L[6] };
eps := EpsilonKLM(eta, k, l, m);
num1 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (V join { k, l })))];
num2 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (W join { k, l })))];
den1 := thetas_sq[TCharToIndex( EtaValue(eta, U sdiff (V join { k, m })))];
den2 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (W join { k, m })))];
return eps*(num1*num2)/(den1*den2);
end function;



/* Final function */
intrinsic RosenhainInvariantsFromThetaSquares(thetas_sq::.) -> .
{Computes Rosenhain invariants.}
v0s := FindDelta(thetas_sq);
if #v0s ne 10 then
    error "Not right number of even zero characteristics:", v0s;
end if;
gamma := ComputeGamma(v0s);
ZZ:=Integers();
        zer:=ZeroMatrix(ZZ, 4,4);
        id:=IdentityMatrix(ZZ,4);
        J:=BlockMatrix([[zer, id],[-id,zer]]);
print gamma*J*Transpose(gamma);
eta := EtaFunction(gamma);
rosens := [ TakaseQuotient(thetas_sq, eta, 1, l, 2) : l in [3..9] ];
return rosens;
end intrinsic;


