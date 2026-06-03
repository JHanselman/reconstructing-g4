load "construct-Lemma-Section-3.2-proofs.m";

list := [[2,3], [1,3], [1,2], [3,4], [2,4], [1,4]] cat &cat[[[j, 4+i]: i in [1..3]]: j in [1..3]] cat [[6,7], [5,7], [5,6]] cat [[4,5], [4,6], [4,7]];
bits := [[GF(37)| 1,0,0], [0,1,0], [0,0,1], [1,1,1], [Random(GF(37)): i in [1..3]], [Random(GF(37)): i in [1..3]],[Random(GF(37)): i in [1..3]]];
quart, bitangents := QuarticFromAronhold(bits);
 R<t0,t1,t2> := PolynomialRing(GF(37),3,"glex");

bitx := [&+[b[i]*R.i: i in [1..3]]: b in bitangents];
quart := R!quart;


count := 0;
for i in [1..7] do
	for j in [i+1..7] do
		for k in [j+1..7] do
			nu := 7 + Index(list, [j,k]);
			mu := 7 + Index(list, [i,k]);
			pro := bitx[i]*bitx[j]*bitx[nu]*bitx[mu];
			I := Radical(Saturation(Ideal([quart, pro])));
			print(2 in [Degree(b): b in Basis(I)]);
			count +:= 1;
		end for;
	end for;
end for;



print count;
