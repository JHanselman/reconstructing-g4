/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */

/*
Code adapted from code written by Jeroen Sijsling which can be found on: 
 https://github.com/JRSijsling/curve_reconstruction/blob/master/magma/rosenhain.m

 */


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

return Transpose(ChangeRing( M^(-1)*(ChangeRing(N, GF(2))), Integers()));
end function;


/* Mumford's eta */
function EtaFunction0(g);
ZZ:=Integers();
zer:=ZeroMatrix(ZZ, g,g);
zer1:=ZeroMatrix(ZZ, g,1);
id:=IdentityMatrix(ZZ,g);
        triang:=zer;
        for i in [1..g] do
                for j in [i..g] do
                        triang[i,j]:=1;
                end for;
        end for;
        M:=VerticalJoin(HorizontalJoin(id, zer1), HorizontalJoin(triang, zer1));
        N:=VerticalJoin(HorizontalJoin(id, zer1), HorizontalJoin(zer1, triang));
	return &cat[[Transpose(Matrix(Transpose(N)[i])), Transpose(Matrix(Transpose(M)[i])) ]: i in [1..g+1]];
end function;


/* Transform Mumford's eta by gamma */
/* Note that these function are represented by tuples, but we give an evaluation later */
function EtaFunction(gamma)
g:=Nrows(gamma)/2;
return [ gamma *v : v in EtaFunction0(g) ];
end function;


/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Integers(), [[ 0 : i in [1..2*g] ]]));
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;


/* Given eta, finds U in Definition 3.3 */
function UFromEtaFunction(eta)
g := #Eltseq(eta[1]) div 2;
return { i : i in [1..#eta] | not IsEvenVector(eta[i]) } join { 2*g + 2 };
end function;



/* This finds the even zero characteristics */
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
g := #Eltseq(eta[1]) div 2;
U := UFromEtaFunction(eta);
Bm := { 1..2*g+1 }; L := [ bp : bp in (Bm diff { k, l, m }) ];
V := { L[i]: i in [1..g-1] }; W := { L[i]: i in [g..2*(g-1)] };
eps := EpsilonKLM(eta, k, l, m);
num1 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (V join { k, l })))];
num2 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (W join { k, l })))];
den1 := thetas_sq[TCharToIndex( EtaValue(eta, U sdiff (V join { k, m })))];
den2 := thetas_sq[TCharToIndex(EtaValue(eta, U sdiff (W join { k, m })))];
return eps*(num1*num2)/(den1*den2);
end function;



/* Final function */
intrinsic RosenhainInvariantsFromThetaSquares(thetas_sq::. , g::RngIntElt : SpTrafo := true) -> .
{Computes Rosenhain invariants.}
if SpTrafo then
   if g eq 4 then
     v0s := FindDelta(thetas_sq);
     if #v0s ne 10 then
       error "Not right number of even zero characteristics:", v0s;
     end if;
     gamma := ComputeGamma(v0s);
     eta := EtaFunction(gamma);
   else
     print "only implemented for g=4";
   end if;
else
   eta:=EtaFunction0(g);
end if;


rosens := [ TakaseQuotient(thetas_sq, eta, 1, l, 2) : l in [3..2*g+1] ];
return rosens;

end intrinsic;


