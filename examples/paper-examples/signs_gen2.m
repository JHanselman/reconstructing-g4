/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


function booltoGF2(bool)
  if bool then
  return GF(2)!1;
  else return GF(2)!0;
end if;
end function;

function chartoInt(cha)
   seq:=Eltseq(cha);
   retpre:=&+[2^(4-i) * Integers()!(GF(2)!seq[i]): i in [1..4]];
   return ((retpre -1) mod 16)+1;
end function;

function preloop()
zer:=ZeroMatrix(GF(2), 2,2);
id:=IdentityMatrix(GF(2),2);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
chars:=EvenThetaCharacteristics(2);
V:=VectorSpace(GF(2), 21);
pivots:= [Index(chars, cha): cha in [[[1,0],[0,0]], [[0,1],[0,0]],[[0,0],[1,0]],[[0,0],[0,1]], [[0,0],[0,0]]]];
pivotscompl:=[i: i in [1..10]| not (i in pivots)];
W:=VectorSpace(GF(2), 4);
T:=[ W!(chars[i][1] cat chars[i][2]) : i in pivots];
Tcompl:=[ W!(chars[i][1] cat chars[i][2]) : i in pivotscompl];
print Tcompl;
bs := [w: w in W| w ne 0 ];
N:=#bs;
vec := [W!(cha[1] cat cha[2]): cha in chars];
cs := [[chars[i] : i in [1..10]|  (Matrix(b)*J*Transpose(Matrix(vec[i])))[1,1] eq (Matrix(b)*J1*Transpose(Matrix(b)))[1,1]]: b in bs  ];
rep := [[c  : c in cs[i] | (c[1] cat c[2])[Position(Eltseq(bs[i]), GF(2)!1)] eq 0] : i in [1..#bs]];
cosets := [[ [W!(rep[i][j][1] cat rep[i][j][2]), (W!(rep[i][j][1] cat rep[i][j][2])+bs[i])]: j in [1..3]   ]  : i in [1..N]];
S1:= [i  : i in [1..N] | &and[cosets[i][1][j] in T: j in [1..2] ]];
S2:= [i  : i in [1..N] | &and[cosets[i][2][j] in T: j in [1..2] ]];
S3:= [i  : i in [1..N] | &and[cosets[i][3][j] in T: j in [1..2] ]];
Si:=S1 cat S2 cat S3;

A1:=[Matrix([ [booltoGF2( &or[cosets[i][2][j] eq c: j in [1..2]])  : c in Tcompl], [booltoGF2(&or[cosets[i][3][j] eq c: j in [1..2]])       : c in Tcompl]]):  i in S1];
A2:=[Matrix([ [booltoGF2( &or[cosets[i][1][j] eq c: j in [1..2]])  : c in Tcompl], [booltoGF2(&or[cosets[i][3][j] eq c: j in [1..2]])       : c in Tcompl]]):  i in S2];
A3:=[Matrix([ [booltoGF2( &or[cosets[i][1][j] eq c: j in [1..2]])  : c in Tcompl], [booltoGF2(&or[cosets[i][2][j] eq c: j in [1..2]])       : c in Tcompl]]):  i in S3];
Ai:=A1 cat A2 cat A3;
print Ai;
As:=A1[1];
i:=2;
list:=[Si[1]];
listcoind:=[1];
while Rank(As) lt 5 do
       VJ:=VerticalJoin(As, Ai[i]);
       if Rank(VJ) gt Rank(As) then
 	 As:=VJ;
	 Append(~list, Si[i]);
	 if i lt #S1 then
	     Append(~listcoind, 1);
	 elif i lt #S1+#S2 then
	     Append(~listcoind, 2);
	 else
	     Append(~listcoind, 3);
         end if;
       end if;
       i:=i+1;
end while;
usedbs:=[bs[i]: i in list];
usedrep:=[rep[i]: i in list];
return As, usedbs, usedrep, listcoind, Tcompl;
end function;



function correct_signs(thetassq)

K := Parent(thetassq[1]);
isexact := (Type(K) ne FldCom);
X, bs, rep, coind, Tcompl:=preloop();
ZZ:=Integers();
W:=Parent(bs[1]);
W1:=RSpace(ZZ,4);
zer:=ZeroMatrix(ZZ, 2,2);
id:=IdentityMatrix(ZZ,2);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
JF2 := ChangeRing(J, GF(2));
J1F2 := ChangeRing(J1, GF(2));
N:=#bs;
cosets := [[ [W1!(rep[i][j][1] cat rep[i][j][2]),  (W1!(rep[i][j][1] cat rep[i][j][2])+W1!bs[i])]: j in [1..3]   ]  : i in [1..N]];
vec:=[];
for i in [1..N] do
        // Term of the theta relation where the sign is fixed.
        coindi:=coind[i];
        // Other two terms
        compl:=[j: j in [1..3]| j ne coindi];
	asposs := [w: w in W| Matrix(w)*JF2*Transpose(Matrix(bs[i]))+Matrix(bs[i])*J1F2*Transpose(Matrix(bs[i])) eq 0 and Matrix(w)*J1F2*Transpose(Matrix(w)) eq 1 ];
	a := W1!asposs[1];
        coeff:=[(-1)^((Matrix(a)*J*Transpose(Matrix(cosets[i][j][1])))[1,1]) *  (-1)^((Matrix(cosets[i][j][1]+W1!bs[i])*J1*Transpose(Matrix(W1!bs[i])))[1,1]): j in [1..3]];
	print coeff;
	sigsposs:=[[coeff[j]: nu in [1..4]]: j in [1..3]  ];
	/* nu=1 means no sign switch
	   nu=2 means 1st sign switches
	   nu=3 means 2nd sign switches
	   nu=4 means both signs switch
	*/
	sigsposs[compl[1]][2] *:= (-1);
	sigsposs[compl[2]][3] *:= (-1);
	sigsposs[compl[1]][4] *:= (-1);
	sigsposs[compl[2]][4] *:= (-1);
	charsInt:=[[chartoInt(cosets[i][j][xi]   )  : xi in [1..2]]: j in [1..3]];
	if isexact then
                thetarel:=[&+[sigsposs[j][nu]*  &*[thetassq[charsInt[j][xi] ]: xi in [1..2]]: j in [1..3]]   : nu in [1..4]];
		ind := [i: i in [1..4]|thetarel[i] eq 0];
		assert #ind eq 1;
		ind := ind[1];
	else
		thetarel:=[Abs(&+[sigsposs[j][nu]*  &*[thetassq[charsInt[j][xi] ]: xi in [1..2]]: j in [1..3]])   : nu in [1..4]];
		min, ind:=Min(thetarel);
		print [sigsposs[j][ind]: j in [1..3]], "\n";
	end if;
        vec cat:= Intseq(ind-1, 2,2);	
end for;
v:=Vector(GF(2), vec);
sol:=Eltseq(Solution(Transpose(X), v));
for j in [1..#Tcompl ] do
      indi:=chartoInt(Tcompl[j]);
      thetassq[indi] *:= (-1)^(ZZ!sol[j]);
end for;
return thetassq;
end function;
