/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */


function is_azygetic(chars)
	id:=IdentityMatrix(GF(2), 4);
        zer:=ZeroMatrix(GF(2), 4,4);
	J:=BlockMatrix(2,2, [zer,id, id,zer]);
	Nchars:=Ncols(chars);
	i:=1;
	for j in [i+1..Nchars] do
		for k in [j+1..Nchars] do
			submat:=Submatrix(chars, [1..8], [i,j,k]);
			mat:= Transpose(submat)*J*submat;
			if mat[1,2]+mat[1,3]+mat[2,3] eq 0 then
				return false;
			end if;
		end for;
	end for;
	return true;
end function;


function Add1at2(vec)
	id:=IdentityMatrix(GF(2), 4);
	zer:=ZeroMatrix(GF(2), 4,4);
	if vec[2] eq 1 then
		return IdentityMatrix(GF(2), 8);
	else
		rhs:=Matrix(GF(2), [[1]]);
		lhs:=Matrix(GF(2),2,1, Eltseq(vec)[7..8]);
		sol:=Solution(lhs, rhs);
                Bx:=Matrix(GF(2), [[0,0,0,0],[0,0, sol[1,1], sol[1,2]],[0,sol[1,1],0,0],[0,sol[1,2], 0,0]]);
                Sx:=BlockMatrix(2,2,[[id, zer  ], [Bx, id]]);
                return Sx;
	end if;
end function;

function map_orthogonal(vec)
	v1:=Vector(Eltseq(vec)[3..4]);
        v2:=Vector(Eltseq(vec)[7..8]);
	id2:=IdentityMatrix(GF(2),2);
	zer2:=ZeroMatrix(GF(2),2,2);
	V2:=VectorSpace(GF(2),2);
	if v1 eq V2!0 and v2 eq V2!0 then
		B:=Matrix(GF(2), [[1,0],[0,0]]);
                ret:=BlockMatrix([[id2, zer2],[B, id2]]);
		retA:=Submatrix(ret, [1..2], [1..2]);
                retB:=Submatrix(ret, [1..2], [3..4]);
                retC:=Submatrix(ret, [3..4], [1..2]);
                retD:=Submatrix(ret, [3..4], [3..4]);
                return BlockMatrix(4,4, [[id2, zer2, zer2, zer2],[zer2, retA, zer2, retB],[zer2, zer2, id2, zer2], [zer2, retC, zer2, retD]]);
	elif v2 eq Vector([GF(2)|0,0]) then
		S1:=IdentityMatrix(GF(2),4);
	elif v1 eq Vector([GF(2)|0,0]) then
		S1:=BlockMatrix([[zer2, id2],[id2, zer2]]);
		v1:=v2;
	else
		B:=Matrix(GF(2), [[0,1],[1,0]]);
	        S1:=BlockMatrix([[id2, B],[zer2, id2]]);
	end if;
	i:=Min([i: i in [1..2]| v1[i] ne 0]);
	A:=Matrix(GF(2), 2,2, [ v1, V2.(3-i)]);
	S2:=BlockMatrix([[A^(-1), zer2], [zer2, Transpose(A)]]);
	ret:=S1*S2;
	retA:=Submatrix(ret, [1..2], [1..2]);
        retB:=Submatrix(ret, [1..2], [3..4]);
        retC:=Submatrix(ret, [3..4], [1..2]);
        retD:=Submatrix(ret, [3..4], [3..4]);
	return BlockMatrix(4,4, [[id2, zer2, zer2, zer2],[zer2, retA, zer2, retB],[zer2, zer2, id2, zer2], [zer2, retC, zer2, retD]]);
end function;

function map_azygetic(azy1, azy2)
	id:=IdentityMatrix(GF(2), 4);
        zer:=ZeroMatrix(GF(2), 4,4);
	J:=BlockMatrix(2,2, [[zer, id], [id, zer]]);
	Rhs:=Vector([GF(2)|0, 0,1]);
        M1:=Matrix([azy1[1]+azy1[2], azy1[1]+azy1[3], azy1[1]+azy1[2]+azy1[3]+azy1[4]]);
	M1:=VerticalJoin(M1, Solution(J*Transpose(M1), Rhs));
	K1:=Matrix(Basis(Kernel(J*Transpose(M1))));
	V1:=VectorSpace(GF(2), 4, K1*J*Transpose(K1));
	HypDec1:=HyperbolicSplitting(V1);
	hypmat1:=Matrix([HypDec1[1][1][1], HypDec1[1][1][2], HypDec1[1][2][1], HypDec1[1][2][2]]);
	M1:=VerticalJoin(M1, hypmat1*K1);
	M1:=Submatrix(M1, [1,3,5,7,2,4,6,8], [1..8]);

	M2:=Matrix([azy2[1]+azy2[2], azy2[1]+azy2[3], azy2[1]+azy2[2]+azy2[3]+azy2[4]]);
	M2:=VerticalJoin(M2, Solution(J*Transpose(M2), Rhs));
	K2:=Matrix(Basis(Kernel(J*Transpose(M2))));
	V2:=VectorSpace(GF(2), 4, K2*J*Transpose(K2));
	HypDec2:=HyperbolicSplitting(V2);
	hypmat2:=Matrix([HypDec2[1][1][1], HypDec2[1][1][2], HypDec2[1][2][1], HypDec2[1][2][2]]);
	M2:=VerticalJoin(M2, hypmat2*K2);
	M2:=Submatrix(M2, [1, 3,5,7,2,4,6,8], [1..8]);
	S:=M1^(-1)*M2;
	error if &or[(Matrix(azy1)*S-Matrix(azy2))[1] ne (Matrix(azy1)*S-Matrix(azy2))[i]: i in [2..4]] , "Error: First symplectic Transformation wrong";
	vec1:=(Matrix(azy1)*M1^(-1))[1];
	vec2:=(Matrix(azy2)*M2^(-1))[1];
	A1:=Submatrix(Transpose(M1^(-1)), [1..4], [1..4]);
	B1:=Submatrix(Transpose(M1^(-1)), [1..4], [5..8]);
	C1:=Submatrix(Transpose(M1^(-1)), [5..8], [1..4]);
	D1:=Submatrix(Transpose(M1^(-1)), [5..8], [5..8]);
    A2:=Submatrix(Transpose(M2^(-1)), [1..4], [1..4]);
	B2:=Submatrix(Transpose(M2^(-1)), [1..4], [5..8]);
	C2:=Submatrix(Transpose(M2^(-1)), [5..8], [1..4]);
	D2:=Submatrix(Transpose(M2^(-1)), [5..8], [5..8]);
	vec1+:= Vector(Diagonal(B1*Transpose(A1)) cat Diagonal(D1*Transpose(C1)));
	vec2+:= Vector(Diagonal(B2*Transpose(A2)) cat Diagonal(D2*Transpose(C2)));
	trans:= vec1+vec2;
	S1add1at3:=Add1at2(vec1);
        S2add1at3:=Add1at2(vec2);
	vec1:=vec1*S1add1at3;
        vec2:=vec2*S2add1at3;
	error if vec1[2] eq 0, "Mistake in function Add1at2()";
        error if vec2[2] eq 0, "Mistake in function Add1at2()";
	Ortho1:=map_orthogonal(vec1);
        Ortho2:=map_orthogonal(vec2);
	/*print "\n vec1, Ortho1, vec1*Ortho1=";
	print vec1, Ortho1, vec1*Ortho1;
        print "\n vec2, Ortho2, vec2*Ortho2=";
	print vec2, Ortho2, vec2*Ortho2;*/
	vec1:=vec1*Ortho1;
	vec2:=vec2*Ortho2;
	A1:=Submatrix(Transpose(Ortho1), [1..4], [1..4]);
        B1:=Submatrix(Transpose(Ortho1), [1..4], [5..8]);
        C1:=Submatrix(Transpose(Ortho1), [5..8], [1..4]);
        D1:=Submatrix(Transpose(Ortho1), [5..8], [5..8]);
	A2:=Submatrix(Transpose(Ortho2), [1..4], [1..4]);
        B2:=Submatrix(Transpose(Ortho2), [1..4], [5..8]);
        C2:=Submatrix(Transpose(Ortho2), [5..8], [1..4]);
        D2:=Submatrix(Transpose(Ortho2), [5..8], [5..8]);
	vec1+:= Vector(Diagonal(B1*Transpose(A1)) cat Diagonal(D1*Transpose(C1)));
        vec2+:= Vector(Diagonal(B2*Transpose(A2)) cat Diagonal(D2*Transpose(C2)));
        error if vec1 ne vec2, "Mistake in map_orthogonal()", Ortho1, Ortho2, vec1, vec2;
	return Transpose(M1^(-1)*S1add1at3* Ortho1*Ortho2^(-1) *S2add1at3^(-1)   *M2);
end function;

intrinsic special_fundamental_system(azy::SeqEnum) -> SeqEnum, SeqEnum, AlgMatElt
{Given an azygetic system of 4 odd characteristics, computes the two sets of 6 even characteristics such that the union is azygetic. It returns a symplectic transformation into the standard system as well  }
	zer:=ZeroMatrix(GF(2), 4,4);
	zer1:=ZeroMatrix(GF(2), 4,1);
	zer2:=ZeroMatrix(GF(2), 4,2);
	id:=IdentityMatrix(GF(2),4);
	J:=BlockMatrix(2,2, [zer,id, id,zer]);
	J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
	triang:=zer;
	for i in [1..4] do
		for j in [i..4] do
			triang[i,j]:=1;
		end for;
	end for;
	M:=VerticalJoin(id, triang);

	N:=VerticalJoin(HorizontalJoin(id, zer2), HorizontalJoin(zer1, HorizontalJoin(triang, zer1)) );
	G:=(Matrix(GF(2),6,6, [1: i in [1..36]])+IdentityMatrix(GF(2), 6));
	N2:=N*G;
	S:=map_azygetic(Transpose(M)[1..4], azy);
	A:=Submatrix(S, [1..4], [1..4]);
	B:=Submatrix(S, [1..4], [5..8]);
	C:=Submatrix(S, [5..8], [1..4]);
	D:=Submatrix(S, [5..8], [5..8]);
	vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
	return [n+vec: n in Transpose(S*N)[1..6]], [n+vec: n in Transpose(S*N2)[1..6]], S;
end intrinsic;

function liftSLN(A, n)
    ZZ:=Integers();
    p:=Characteristic(BaseRing(A));
    assert(BaseRing(A) eq GF(p));
    N:=Nrows(A);
    Ainv:=A^(-1);
    col:=[i : i in [1..N]|Ainv[1,i] ne 0][1];
    entry:=ZZ!(Ainv[1,col]^(-1));
    Alift:=A;
    for i in [2..n] do
	Alift:=ChangeRing(ChangeRing(Alift, ZZ), Integers(p^i));
	eps:=ZZ!(Determinant(Alift)-1);
	Alift[col,1] -:= eps*entry;
    end for;
    return Alift;
end function;




function decompose_symplectic(S)
     ZZ:=Integers();
     ZZ4:=Integers(4);
     N:=Nrows(S) div 2;
     id:=IdentityMatrix(GF(2), N);
     zer:=ZeroMatrix(GF(2), N, N);
     idZ:=IdentityMatrix(ZZ, N);
     zerZ:=ZeroMatrix(ZZ, N, N);
     J:=BlockMatrix(2,2,[[zerZ, idZ], [-idZ, zerZ]]);
     A:=Submatrix(S, [1..N], [1..N]);
     B:=Submatrix(S, [1..N], [N+1..2*N]);
     C:=Submatrix(S, [N+1..2*N], [1..N]);
     D:=Submatrix(S, [N+1..2*N], [N+1..2*N]);
     Ech, T1:=EchelonForm(C);
     Diag, T2:=EchelonForm(Transpose(Ech));
     T1lift:=liftSLN(T1, 3);
     T2lift:=liftSLN(T2,3);
     ret1:=[ChangeRing(DiagonalJoin(Transpose(T1lift), T1lift^(-1)), ZZ)];
     ret2:=[ChangeRing(DiagonalJoin(Transpose(T2lift)^(-1),  T2lift), ZZ)];
     trafo1:=DiagonalJoin(Transpose(T1) , T1^(-1));
     trafo2:=DiagonalJoin(Transpose(T2)^(-1), T2);
     S:=trafo1^(-1)*S*trafo2^(-1);
     A:=Submatrix(S, [1..N], [1..N]);
     B:=Submatrix(S, [1..N], [N+1..2*N]);
     C:=Submatrix(S, [N+1..2*N], [1..N]);
     D:=Submatrix(S, [N+1..2*N], [N+1..2*N]);
     nu:=Max([0] cat [i: i in [1..N]| Diag[i,i] ne 0]);
     X:=IdentityMatrix(GF(2), nu)-Submatrix(A, [1..nu], [1..nu]);
     X:=DiagonalJoin(X, ZeroMatrix(GF(2), N-nu, N-nu));
     trafo:=BlockMatrix(2,2,[[id, X], [zer, id]]);
     Append(~ret1, ChangeRing(trafo, ZZ));
     S:=trafo^(-1)*S;
     A:=Submatrix(S, [1..N], [1..N]);
     Alift:=liftSLN(A,3);
     trafo:=DiagonalJoin(A , Transpose(A)^(-1));
     Append(~ret1, ChangeRing(DiagonalJoin(Alift, Transpose(Alift)^(-1)), ZZ));
     S:=trafo^(-1)*S;
     C:=Submatrix(S, [N+1..2*N], [1..N]);
     trafo:=BlockMatrix(2,2,[[id, zer], [C, id]]);
     Append(~ret1, J);
     Append(~ret1, ChangeRing(BlockMatrix(2,2,[[id, C], [zer, id]]), ZZ));
     Append(~ret1, J);
     kappa := (-1)^N;
     S:=trafo^(-1)*S;
     Append(~ret2, ChangeRing(S, ZZ));
     return ret1 cat Reverse(ret2), kappa;
end function;

intrinsic signs_in_derivative_formula(S:: AlgMatElt)-> SeqEnum
{Computes the signs on the right hand side of the generalized Jacobi derivative formula}
	dec, kappa:=decompose_symplectic(Transpose(S));
	psi1:=0;
	psi2:=0;
	ZZ:=Integers();
	zer:=ZeroMatrix(ZZ, 4,4);
	zer1:=ZeroMatrix(ZZ, 4,1);
	zer2:=ZeroMatrix(ZZ, 4,2);
	id:=IdentityMatrix(ZZ,4);
	J1:=BlockMatrix([[zer, id],[zer,zer]]);
	triang:=zer;
	for i in [1..4] do
		for j in [i..4] do
			triang[i,j]:=1;
		end for;
	end for;
	M:=VerticalJoin(id, triang);
	N:=VerticalJoin(HorizontalJoin(id, zer2), HorizontalJoin(zer1, HorizontalJoin(triang, zer1)) );
        G:=(Matrix(ZZ,6,6, [1: i in [1..36]])+IdentityMatrix(ZZ, 6));
        N2:=N*G;
        for de in dec do
        	Si:=Transpose(de);
        	A:=Submatrix(Si, [1..4], [1..4]);
     		B:=Submatrix(Si, [1..4], [5..8]);
	     	C:=Submatrix(Si, [5..8], [1..4]);
	     	D:=Submatrix(Si, [5..8], [5..8]);
	     	
	     	Mpre:=Si*M;
	     	Npre:=Si*N;
	     	N2pre:=Si*N2;
	     	vecpre:= -2*Vector([0,0,0,0] cat Diagonal(D*Transpose(C)));
	     	Mpre +:= Transpose(Matrix([vecpre: i in [1..4]]));
		Npre +:= Transpose(Matrix([vecpre: i in [1..6]]));
		N2pre +:= Transpose(Matrix([vecpre: i in [1..6]]));
		psi1+:=Trace(Transpose(Mpre)*J1*Mpre)-Trace(Transpose(M)*J1*M)-Trace(Transpose(Npre)*J1*Npre)+Trace(Transpose(N)*J1*N);
		psi2+:=Trace(Transpose(Mpre)*J1*Mpre)-Trace(Transpose(M)*J1*M)-Trace(Transpose(N2pre)*J1*N2pre)+Trace(Transpose(N2)*J1*N2);		
		M:=Si*M;
	     	N:=Si*N;
	     	N2:=Si*N2;
	     	vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
		M +:= Transpose(Matrix([vec: i in [1..4]]));
		N +:= Transpose(Matrix([vec: i in [1..6]]));
		N2 +:= Transpose(Matrix([vec: i in [1..6]]));
	end for;
	carryM:= M div 2;
	carryN:= N div 2;
	carryN2:= N2 div 2;	
	switch1:=(-1)^(Trace(Transpose(M)*J1*carryM)+Trace(Transpose(N)*J1*carryN));
	switch2:=(-1)^(Trace(Transpose(M)*J1*carryM)+Trace(Transpose(N2)*J1*carryN2));
        bool1, quo1:=IsDivisibleBy(psi1,4);
        bool2, quo2:=IsDivisibleBy(psi2,4);
        if bool1 and bool2 then
	       	return [kappa*(-1)^(quo1)*switch1, -kappa*(-1)^(quo2)*switch2]; //TODO: Find initial signs
        else
       	 	error("Not divisible by 4");
	end if;
end intrinsic;





	
/* For testing:

V:=VectorSpace(GF(2), 8);
A:=Matrix(GF(2), [[1: j in [1..2*i-1]] cat [0: j in [1..8-2*i+1]  ]: i in [1..4] ] cat [[1: j in [1..2*i-2]] cat [0, 1] cat [0: j in [1..8-2*i]  ]: i in [1..4] ]  );
zer:=ZeroMatrix(GF(2), 4,4);
id:=IdentityMatrix(GF(2),4);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
U:=Transpose(A)^(-1);
Transpose(U)*A*U;
v:=Matrix(GF(2), 8,1, [1: i in [1..8]]);
U:=HorizontalJoin(U, v);
new:=[];
for v in V do
        X:=Transpose(Matrix([v: i in [1..9 ]]   ));
        diag:=Diagonal(Transpose(X+U)*J1*(X+U)  );
        Append(~new, &+[Integers()!diag[i]: i in [1..9]]);
end for;
_, ind:=Max(new);
v:=[v: v in V][ind];
X:=Transpose(Matrix([v: i in [1..9 ]]   ));
N:=X+U;

N3:=Transpose(N)[1..3];
q4:=Transpose(N)[4];
q6:=Transpose(N)[6];
for i in [1..3] do
        N3[i] +:= q4;
end for;
M:=Transpose(Matrix(N3));
K:=Kernel(J*M);
K1:=Kernel(Transpose(M)*J*M);

cos:=[Matrix(8,1, Eltseq(q6+k)): k in K];
sys:=[x: x in cos| (Transpose(x)*J1*x)[1,1] eq 1];
eta:=M*Matrix(3, 1, Basis(K1)[1]);
sys1:=[s:s in sys | s[1,1] eq 0];
w:=Basis(K1)[1]*Transpose(M);
T:=Matrix(GF(2), [[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
 Blo:=BlockMatrix([[(Transpose(T))^(-1), zer],[zer, T]]);

T2:=Matrix(GF(2), [[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]]);
Blo2:=BlockMatrix([[id, zer],[T2, id]]);
Ortho:=Blo2*Blo;
sysnew:=[Transpose(Ortho)*sys[i]: i in [1..20]];
sys1new:=[Eltseq(s): s in sysnew| s[5,1] eq 0];
proj:=[[s[2..4],s[6..8]]: s in sys1new];




for i in [1..20] do
	for k in [1..4] do
		temp_list:=Transpose(N)[1..k-1] cat [Vector(sys[i])] cat Transpose(N)[k+1..4];
		N1, N2, S:=special_fundamental_system(temp_list);
		Special1:=Transpose(Matrix(temp_list cat N1));
		Special2:=Transpose(Matrix(temp_list cat N2));
		print is_azygetic(Special1);
		print is_azygetic(Special2);
		if not is_azygetic(Special2) then
			print i,k;
		end if;
		print "\n";
	end for;
end for;



S:=map_azygetic(Transpose(N)[1..4], Transpose(N)[1..3] cat [Vector(sys[1])]);
A:=Submatrix(S, [1..4], [1..4]);
B:=Submatrix(S, [1..4], [5..8]);
C:=Submatrix(S, [5..8], [1..4]);
D:=Submatrix(S, [5..8], [5..8]);
vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
print [n+vec: n in Transpose(S*N)[1..4]];


S:=map_azygetic(Transpose(N)[1..4], Transpose(N)[1..3] cat [Vector(sys[3])]);
A:=Submatrix(S, [1..4], [1..4]);
B:=Submatrix(S, [1..4], [5..8]);
C:=Submatrix(S, [5..8], [1..4]);
D:=Submatrix(S, [5..8], [5..8]);
vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
print [n+vec: n in Transpose(S*N)[1..4]];
*/

