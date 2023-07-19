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
