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
	        print bitangents, bits;			
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
                        C := Curve(Scheme(ProjectiveSpace(R4), Eqs));
                        if Genus(C) ne 4 then
                                continue;
                        end if;
			R4:= Universe(Eqs);
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
			print "found envelopping cones";
			R1<t> := PolynomialRing(K);
			pencil := [ChangeRing(qc, R1)+t*ChangeRing(env, R1): env in envcon];
			detpencil := [Determinant(penc): penc in pencil];
			pols := [Gcd(poly, Derivative(poly)): poly in detpencil];
			assert &and[Degree(po) eq 1: po in pols];
			roots := [Roots(po)[1][1]: po in pols];
			tritpairs := [Evaluate(pencil[i], roots[i]): i-> ro in roots];
			return Eqs, tritpairs;
		end for;
	end while;
end function;
