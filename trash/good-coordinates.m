/*
 Tools for arithmetic construction from small period matrices
 */


function MapBinaryForm(f1, f2)
//Given two equivalent binary forms of even degree f1, f2. They should be univariate polynomials. It returns a transformation matrix.     
    J1 := AnalyticJacobian(f1);
    J2 := AnalyticJacobian(f2);
    bool, _, M := IsIsomorphic(J1, J2);
    if bool then
	    return Submatrix(M, [2,1], [2,1]);
    else 
	   error "Two forms are not equivalent.";
    end if;
    
    
end function;

function mapPGl2(Tup1, Tup2)
	X := Matrix(&cat[  [[Tup1[i],1, 0,0] cat [0: j in [1..(i-1)]] cat  [Tup2[i]] cat [0: j in [1..(3-i)]],   [0,0,Tup1[i],1] cat [0: j in [1..(i-1)]] cat  [1] cat [0: j in [1..(3-i)]]]: i in [1..3]]);
	K := NumericalKernel(Transpose(X));
	return Eltseq(K)[1..4];
end function;

function MapBinaryForm2(f1, f2)
	R:=Parent(f1);
        t:=R.1;
	N:=Degree(f1);
	prec:=Precision(BaseRing(R));
        f1 /:= LeadingCoefficient(f1);
	f2 /:= LeadingCoefficient(f2);
	print f1,f2;
	roo1:=Roots(f1);
	roo2:=Roots(f2);
	for i in [1..N] do
		for j in [1..i-1] cat [i+1..N] do 
			for k in [1..j-1] cat [j+1..N] do
			      a,b,c,d:=Explode(mapPGl2([roo1[nu][1]: nu in [1..3]], [roo2[nu][1]: nu in [i,j,k]]));
  			      ftry:=&+[Coefficient(f2, i)*(a*t+b)^i*(c*t+d)^(N-i): i in [0..N]];
			      ftry /:= LeadingCoefficient(ftry);
			      err := &+[Abs(c): c in  Coefficients(ftry-f1)];
			      if err lt 10^(-prec/2) then
			             return Matrix(2,2, [a,b,c,d]);
		              end if;
		         end for;
	         end for;
         end for;
end function;

function ApplyAut(f, sig)
	coeffs, mons := CoefficientsAndMonomials(f);
	return &+[sig(coeffs[i]) * mons[i]: i in [1..#mons]];
end function;

function ConjugateMatrix(M, sig)
	n := Nrows(M);
	m := Ncols(M);
	return Matrix(n, m, [sig(el): el in Eltseq(M)]);
end function;

function ConjugateMap(map, sig)
	X := Domain(map);
	Y := Codomain(map);
	Equ := DefiningPolynomials(map);
	EquInv := DefiningPolynomials(Inverse(map));
	return iso< X-> Y | [ApplyAut(f, sig) : f in Equ], [ApplyAut(f, sig) : f in EquInv] >;
end function;


function IsIsomorphicConic1(C1, C2)
	Leg1, iso1 := LegendreModel(C1);
	Leg2, iso2 := LegendreModel(C2);
	A1 := QuaternionAlgebra(Leg1);
	A2 := QuaternionAlgebra(Leg2); 
	bool, m := IsIsomorphic(A1, A2);
	if not bool then
	       return false;
	else
		Mat := Matrix([Coordinates(m(A1.i))[2..4]: i in [1..3]]);

		return true, ...;
	end if;
end function;

function IsIsomorphicConic(C1, C2)
	K:=BaseRing(C1);
        P1:=Curve(ProjectiveSpace(K, 1));
	bool, p1 := HasRationalPoint(C1);
        if bool then
            bool2, p2 := HasRationalPoint(C2);
	    if bool2 then
	       map1 := Parametrization(C1, p1, P1);
               map2 := Parametrization(C2, p2, P1);
 	       Equ1 := DefiningPolynomials(map1);
               Equ2 := DefiningPolynomials(map2);
               Mat1 := Matrix([[MonomialCoefficient(f,m): f in Equ1]:  m in [P1.1^2, P1.1*P1.2, P1.2^2]] );
               Mat2 := Matrix([[MonomialCoefficient(f,m): f in Equ2]:  m in [P1.1^2, P1.1*P1.2, P1.2^2]] );
               return true, Mat2*Inverse(Mat1);
             else
               return false;
             end if;
         else
             R<t> := PolynomialRing(K);
             f := Evaluate(DefiningEquation(C1), [0, t, 1]);
             L := ext<K| f>;
             Labs := AbsoluteField(L);
             emb := EmbeddingMap(Labs, L);
             sig := hom<L->L| L.1-Trace(L.1)>;
             C1L := ChangeRing(C1, L);
             C2L := ChangeRing(C2, L);
             C2Labs := ChangeRing(C2, Labs);
             P1L := ChangeRing(P1, L);
             p1:= C1L![0, L.1, 1];
             bool, p2 := HasRationalPoint(C2Labs);
             if bool then
                p2 := C2L![emb(a): a in Eltseq(p2)];
                p1conj := C1L![0, sig(L.1), 1];
                p2conj := C2L![sig(a): a in Eltseq(p2)];
                map1 := Parametrization(C1L, p1, P1L);
                map2 := Parametrization(C2L, p2, P1L);
                map1conj := ConjugateMap(map1, sig);
                map2conj := ConjugateMap(map2, sig);
                phi1 := map1conj * Inverse(map1);
                phi2 := map2conj * Inverse(map2);
		Equ1 := DefiningPolynomials(phi1);
                Equ2 := DefiningPolynomials(phi2);
		gcd1 := Gcd(Equ1);
                gcd2 := Gcd(Equ2);
                Equ1 := [f div gcd1: f in Equ1];
                Equ2 := [f div gcd2: f in Equ2];

		X1 := Matrix([[MonomialCoefficient(f, m): f in Equ1]:  m in [P1L.i: i in [1..2]]]);
                X2 := Matrix([[MonomialCoefficient(f, m): f in Equ2]:  m in [P1L.i: i in [1..2]]]);
		bool1, diag1 := IsDiagonal(X1 * ConjugateMatrix(X1, sig));
                bool2, diag2 := IsDiagonal(X2 * ConjugateMatrix(X2, sig));
		assert(bool1 and bool2);
		assert(#Seqset(diag1) eq 1 and #Seqset(diag2) eq 1);
		lambda1:= diag1[1];
		lambda2:= diag2[1];




             else
                return false;
             end if;
          end if; 
end function;




/*
 *
 * For testing:
 R<t> := PolynomialRing(Rationals());
 RC := PolynomialRing(ComplexField(100));
 f1 := t^6 + &+[Random(-10,10) * t^i: i in [0..5]];
 a,b,c,d := Explode([Random(-10,10): i in [1..4]]);
 f2 := (c*t+d)^6 * Evaluate(f1, (a*t+b)/(c*t+d));
 f2 := R!f2;                                   
 M:=MapBinaryForm(RC!f1, RC!f2);
 Minv:=M^(-1);
 print a*Minv/Minv[1,1];
 print Matrix(2,2, [a,b,c,d]);
 *
 *
 */
