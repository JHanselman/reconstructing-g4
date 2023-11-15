function RecillasConstruction(X, points)
	//Given a plane quartic X and 4 points on it (not on a line). Constructs a genus 4 curve C with a two torsion point eta such that the Prym is X 
	K:=BaseField(X);
	R<x0,x1,x2>:=CoordinateRing(Ambient(X));
	Quart:=DefiningEquations(X)[1];
	Mat:=Matrix( &cat[[[(Matrix(3,1, points[i])*Matrix(1,3,points[i]))[mu,nu]: i in [1..4]]: mu in [nu..3]]: nu in [1..3]]  );
	Bas:=Basis(Kernel(Mat));
	count:=1;
	quad:=[ZeroMatrix(K,3,3): i in [1..4]];
	for nu in [1..3] do
		for mu in [nu..3] do
			quad[1][mu,nu]:=Bas[1][count];
			quad[2][mu,nu]:=Bas[2][count];
			count+:=1;
		end for;
	end for;
	quad[1]:=quad[1]+Transpose(quad[1]);
        quad[2]:=quad[2]+Transpose(quad[2]);
	xs:=Matrix(R, 3,1, [x0, x1, x2]);
	quads:= [(Transpose(xs)*ChangeRing(qua,R)*xs)[1,1]: qua in quad[1..2]];
	I:=ideal<R|[Quart, quads[1]]>;
	pointsSch:=&*[ideal<R| &cat[[poi[i]*R.j-poi[j]*R.i: i in [j+1..3]]: j in [1..3]]> : poi in points];
	J:= Saturation(I, pointsSch);
	J:= Saturation(J, ideal<R|[R.1, R.2, R.3]>);
	gen2 := [g: g in Generators(J)| (Degree(g) eq 2) and (not IsDivisibleBy(g, quads[1]))];
	Append(~quads, gen2[1]);
	I2:=ideal<R|[Quart, quads[2]*quads[3]]>;
	J2 := Saturation(I2, quads[1]);
        J2:= Saturation(J2, ideal<R|[R.1, R.2, R.3]>);
        gen2 := [g: g in Generators(J2)| (Degree(g) eq 2)];
	Append(~quads, gen2[1]);
	for nu in [1..3] do
                for mu in [nu..3] do
                        quad[3][mu,nu] := MonomialCoefficient(quads[3], R.nu * R.mu);
                        quad[4][mu,nu] := MonomialCoefficient(quads[4], R.nu * R.mu);
                        count +:= 1;
                end for;
        end for;
        quad[3]:=quad[3] + Transpose(quad[3]);
        quad[4]:=quad[4] + Transpose(quad[4]);
	mons4:=MonomialsOfDegree(R, 4);
	Mat := Matrix([[ MonomialCoefficient(poly, m): m in mons4]: poly in [quads[1]*quads[4], quads[2]*quads[3], Quart]]);
	Bas:= Basis(Kernel(Mat));
	assert #Bas eq 1;
	R4<t0,t1,t2,t3>:=PolynomialRing(K,4);
	A:= &+[R4.i*ChangeRing(quad[i], R4): i in [1..4]]; 
	return [Determinant(A), Bas[1][1]^(-1)*t0*t3 + Bas[1][2]^(-1)*t1*t2], quad;
end function;
load "universal_test.m";

function RationalTrits(p)
	K := GF(p);
	while true do
		bits := [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [Random(K), Random(K), 1], [Random(K), Random(K), 1], [Random(K), Random(K), 1]];
		if &or[mi eq 0: mi in Minors(Matrix(bits),3)] then
			continue;
		end if;
		Quart, bitangents := QuarticFromAronhold(bits);
		if not IsSmoothHyperSurface(Quart) then
			continue;
		end if;
		X := Scheme(ProjectiveSpace(Parent(Quart)), Quart);
		print "found quartic!";
		points := Setseq(Points(X));
		for S in Subsets(Seqset(points), 4) do
			poi := [Coordinates(P): P in S];
			if Rank(Matrix(poi)) lt 3 then
				continue;
			end if;
			Eqs, quad := RecillasConstruction(X, poi);
			print "found genus 4 curve";
			R4:= Universe(Eqs);
			C:=Curve(Scheme(ProjectiveSpace(R4), Eqs));
			if Genus(C) ne 4 then
				continue;
			end if;
			R3:= PolynomialRing(K, 3);
			adjo := Matrix(R3, [[&+[quad[nu][i,j]*R3.j: j in [1..3]]: nu in [1..4]]: i in [1..3]]);
			c := [Determinant(VerticalJoin(adjo, Matrix(RSpace(R3, 4).i))): i in [1..4]];
			nodes := PointsOverSplittingField(SingularSubscheme(Scheme(ProjectiveSpace(R4), Eqs[1])));
			F := Parent(nodes[1][1]);
			R3F := ChangeRing(R3, F);
			ni := [Basis(Image(&+[nodes[i][j]*quad[j]: j in [1..4]]))[1]: i->n in nodes];
			if #ni ne 4 then
				continue;
			end if;
			print "generic symmetroid";
			Ni := R3!(&*[&+[ni[i][j]*R3F.j: j in [1..3]]: i in [1..4]]);
			Sex := [Ni * (&+[bitangents[i][j]*R3.j : j in [1..3]])^2: i->b in bitangents];
			mons := MonomialsOfDegree(R3, 6);
			Mat := Matrix(&cat[[[MonomialCoefficient(c[i]*c[j], m): m in mons]: i in [j..4]]: j in [1..4]]);
			Rhs := Matrix([[MonomialCoefficient(sex, m): m in mons]: sex in Sex]);
			sol := Solution(Mat, Rhs);
			envcon := [ZeroMatrix(K, 4,4): nu in [1..28]];
			qc := ZeroMatrix(K,4,4);
			count:=1;
			for i in [1..4] do
				for j in [i..4] do
					for nu in [1..28] do
						envcon[nu][i,j]+:=sol[nu][count];
					end for;
					qc[i,j] +:= MonomialCoefficient(Eqs[2], R4.i*R4.j);
					count +:=1 ;
				end for;
			end for;
			envcon := [Transpose(env)+env: env in envcon];
			qc := qc + Transpose(qc);
			print "found enveloping cones";
			R1<t> := PolynomialRing(K);
			pencil := [ChangeRing(qc, R1)+t*ChangeRing(env, R1): env in envcon];
			detpencil := [Determinant(penc): penc in pencil];
			pols := [Gcd(poly, Derivative(poly)): poly in detpencil];
			roots := [[ro[1]: ro in Roots(po)]: po in pols];
			assert &and[#ro eq 1: ro in roots];
			roots := [ro[1]: ro in roots];
			tritpairs := [Evaluate(pencil[i], roots[i]): i-> ro in roots];
			R4p2 := ChangeRing(R4, GF(p^2));
                        R3p2 := ChangeRing(R3, GF(p^2));
			tritpairspol := [(Matrix(1, 4, [R4p2.i: i in [1..4]]) * ChangeRing(tri, R4p2) * Matrix(4, 1, [R4p2.i: i in [1..4]]))[1,1]: tri in tritpairs];
			tris := [[e[1]: e in Factorization(tri)]: tri in tritpairspol];
			tricoeffs := [[[MonomialCoefficient(tri, R4p2.i) : i in [1..4]]: tri in trispa]: trispa in tris];
			Hi := [[&+[coeff[i]* R3p2!c[i]: i in [1..4]]: coeff in tricoeff]: tricoeff in tricoeffs];
			return Hi, bitangents, Ni;
		end for;
	end while;
end function;

function Enriques(p)
	Hi, bitangents, Ni := RationalTrits(p);
	R3 := Universe(Hi[1]);
	IntPoints := AssociativeArray();
	for S in Subsets(Set(CartesianProduct([1..28], [1..2])), 2) do
	 	Slist:= Setseq(S);
		i1 := Slist[1][1];
		j1 := Slist[1][2];
		i2 := Slist[2][1]; 
		j2 := Slist[2][2];
		bit1 := bitangents[i1];
		bit2 := bitangents[i2];
		Hi1 := Hi[i1][j1];
                Hi2 := Hi[i2][j2];
		I:= Ideal([Hi1, Hi2]);
		I := Saturation(I, Ideal(R3!Ni));
		Pois := Setseq(PointsOverSplittingField(Scheme(ProjectiveSpace(R3), I)));
		assert #Pois eq 3;
		IntPoints[{i1+(j1-1)*28, i2+(j2-1)* 28}] := Pois;
	end for;
	Gra := AssociativeArray();
	for i1 in [1..55] do
		for i2 in [i1+1..56] do
			if i2 eq i1+1 then
				Gra[{[i1,1], [i2,1]}]:=1;
                		Gra[{[i1,2], [i2,2]}]:=1;
			else
				Gra[{[i1,1], [i2,1]}]:=0;
		                Gra[{[i1,2], [i2,2]}]:=0;
			end if;
			Gra[{[i1,1], [i2,2]}]:=0;
                        Gra[{[i1,2], [i2,1]}]:=0;
			for nu in [1..3] do
				if (nu eq 1) and (i2 eq i1+1) then
					continue;
				end if;
				swa :=	&+[bitangents[i2 mod 28][mu]*IntPoints[{i1, i2}][nu][mu]: mu in [1..3]]/&+[bitangents[i1 mod 28][mu]*IntPoints[{i1, i2}][nu][mu]: mu in [1..3]];
				for loo in [i1..i2-1] do
					swa *:= &+[bitangents[loo mod 28][mu]*IntPoints[{loo, loo+1}][1][mu]: mu in [1..3]]/&+[bitangents[(loo+1) mod 28][mu]*IntPoints[{loo, loo+1}][1][mu]: mu in [1..3]];
				end for;
				if swa eq 1 then
					Gra[{[i1,1], [i2,1]}]+:=1;
                			Gra[{[i1,2], [i2,2]}]+:=1;
				elif swa eq -1 then
					Gra[{[i1,1], [i2,2]}]+:=1;
			                Gra[{[i1,2], [i2,1]}]+:=1;
				else 
					error "swa not pm 1";
				end if;
			end for;
		end for;
	end for;

	return IntPoints, Gra;
end function;

for p in PrimesInInterval(37,1000) do
	for i in [1..100] do
		try
			_, Gra := Enriques(101);
		catch e
			continue;
		end try;
	end for;
end for;
A := ZeroMatrix(Integers(),112,112);
for i1 in [1..112] do
	for i2 in [1..112] do
	A[i1,i2] := Gra[{[i1 mod 56, i1 div 56+1], [i2 mod 56, i2 div 56+1]}];
	end for;
end for;
PrintFile("IntMat.txt", A);

