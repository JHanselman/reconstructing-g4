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



/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Integers(), [[ 0 : i in [1..2*g] ]]));
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;

f := 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;

function ComputeGaloisAction(f)
  
  L := SplittingField(f);
  M, mM := OptimisedRepresentation(L);
  G,S, phi := AutomorphismGroup(M);
  gens := Setseq(Generators(G));
  roots := [r[1]: r in Roots(f, L)];
  Sort(~roots);
  ActionOnTheta := [];
  Lt := [6,9,12,15,16];
  U:={1,3,5};
  eta := EtaFunction0(2);
  
  for g in gens do
    sigma_g := mM * phi(g) * mM^(-1);
    g_roots := [sigma_g(r): r in roots];
    L := [1..#roots];
    ParallelSort(~g_roots, ~L);
    if #L eq 5 then
      Append(~L, 6);
    end if;
    
    Append(~ActionOnTheta, [0,0,0,0,0]);
    
    for i in [1..5] do
      for j in [i+1..5] do
        T := {i,j, 6};	
        S := T sdiff U;
        cha := EtaValue(eta, Setseq(S));
        cha_i := TCharToIndex(cha);
        if cha_i in Lt then
          pos := Position(Lt, cha_i); 
          Tnew := {L[i], L[j], L[6]};
          Snew := Tnew sdiff U;
          chanew := EtaValue(eta, Setseq(Snew));
          cha_inew := TCharToIndex(chanew);
          ActionOnTheta[#ActionOnTheta][pos]:= cha_inew;
        end if;
      end for;
    end for;
    
  end for;
  return ActionOnTheta;
end function;
