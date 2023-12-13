function MinimizeReduceQuadric(Q)
	g := Nrows(Q);
	S := IdentityMatrix(Rationals(), g);
	det := Determinant(Q);
	fac := Factorization(det);
	primes := [fa: fa in fac|fa[2] gt 1];
	for pri in primes do
	 	p := pri[1];
		ei := pri[2];
		while ei gt 1 do
			Qred := ChangeRing(Q, GF(p));
			ker := KernelMatrix(Qred);
			nro := Nrows(ker);
			if nro ge ei then
				continue pri;
			end if;
			kerZ := ChangeRing(ker, Integers());
			D, U, V := SmithForm(kerZ);
			SNF := V^-1;
			S := SNF * S;
			Q := SNF * Q * Transpose(SNF);
			j := [i: i in [1..nro]| Valuation(Q[i,i],p) gt 1];
			if j ne [] then
				j := j[1];
			else
				Qsub := Submatrix(Q, 1, 1, nro, nro) div p;
				Qsub := ChangeRing(Qsub, GF(p));
				kersub := KernelMatrix(Qsub);
				kersubZ := ChangeRing(kersub, Integers());
	                        D, U, V := SmithForm(kersubZ);
	                        SNF := V^-1;
				Stra := Matrix(BlockDiagMat(<SNF, IdentityMatrix(Integers(), g-nro)>));
				S := Stra * S;
				Q := Stra * Q * Transpose(Stra);
				j := 1;
			end if;
			Dia := DiagonalMatrix([1: i in [1..j-1]] cat [1/p] cat [1: i in [1..g-j]]);
                        S := Dia * S;
                        Q := Dia * ChangeRing(Q, Rationals()) * Dia;
                        Q := ChangeRing(Q, Integers());
                        ei -:= 2;
		end while;
	end for;
	R := PolynomialRing(Rationals(), g);
	X := Matrix([[R.i: i in [1..g]]]);
	_, S1 := ReduceQuadrics([(X * ChangeRing(Q, R) * Transpose(X))[1,1]]);
	Q := S1 * Q * Transpose(S1);
	return Q, S1*S;
end function;

function MinimizeG4(Q, Gamma)
	R4 := Parent(Q);
	Qmat := ZeroMatrix(Rationals(), 4);
	for i in [1..4] do
		for j in [i..4] do
			Qmat[i,j] := MonomialCoefficient(Q, R4.i*R4.j);
		end for;
	end for;
	Qmat := Qmat+Transpose(Qmat);
	Qmat *:= LCM([Denominator(q): q in Eltseq(Qmat)]);
	Qmat := ChangeRing(Qmat, Integers());
        Qmat /:= GCD(Eltseq(Qmat));
	Qmat1, S := MinimizeReduceQuadric(Qmat);
        X := Matrix([[R4.i: i in [1..4]]]);
	Gamma1 := Evaluate(Gamma, Eltseq(X*ChangeRing(S, R4)));
	Gamma1 *:= LCM([Denominator(q): q in Coefficients(Gamma1)]);
        Gamma1 /:= GCD([Numerator(q): q in Coefficients(Gamma1)]);
	Q1 := (X * ChangeRing(Qmat1, R4) * Transpose(X))[1,1];
        Q1 /:= GCD([Numerator(q): q in Coefficients(Q1)]);

	mons3 := MonomialsOfDegree(R4, 3);
	LMat := Matrix([[MonomialCoefficient(pos,m ): m in mons3]: pos in [Gamma1] cat [Q1*R4.i: i in [1..4]]]);
	QL := Matrix(LMat[2..5]);
	D, U, V := SmithForm(ChangeRing(QL, Integers()));
	SNF := V^-1;
	vec := ChangeRing(LMat[1], Integers()) * V;
	Gcdvec := GCD([vec[i]: i in [5..20]]);
	vec := Vector([0,0,0,0] cat [vec[i] div Gcdvec: i in [5..20]]) * SNF;
	L := Lattice(QL);
	proj := Projection(L, vec);
	ratio := Norm(vec)/Norm(vec-proj);
	bound := Floor(Sqrt(ratio));
	vecopt := vec;
	i := 1;
	while i lt bound do
		clos := ClosestVector(L, i * proj);
		if Norm(i*vec - clos) lt Norm(vecopt) then
			vecopt := i * vec - clos;
			ratio := Norm(vecopt)/Norm(vec - proj);
		        bound := Floor(Sqrt(ratio));
		end if;
		i+:=1;
	end while;
	Gammanew := &+[vecopt[i]*m: i->m in mons3];
	return Q1, Gammanew, S;
end function;

for i in [1..1] do
	R4 := PolynomialRing(Rationals(), 4);
        mons2 := MonomialsOfDegree(R4, 2);
	mons3 := MonomialsOfDegree(R4, 3);
	Q := &+[Random(1000)*m: m in mons2];
        Gamma := &+[Random(1000)*m: m in mons3];
	Q1, Gamma1, S := MinimizeG4(Q, Gamma);
	print Q, Gamma;
	print Q1, Gamma1, "\n";
end for;

