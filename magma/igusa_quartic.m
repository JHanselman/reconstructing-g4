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
    return [[ 0 : i in [1..g] ], [ 0 : i in [1..g] ]];
end if;
res:=Eltseq(&+[ eta[i] : i in S ]);
return [res[1..g], res[g+1..2*g]];
end function;

//f := 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;

function ComputeGaloisAction(f)
  L := SplittingField(f);
  if Characteristic(L) eq 0 then
          M, mM := OptimisedRepresentation(L);
          G,S, phi := AutomorphismGroup(M);
          gens := Setseq(Generators(G));
  else
          p:= Characteristic(L);
          K := BaseRing(Parent(f));
          q:= p^Degree(K);
	  if L eq BaseRing(Parent(f)) then
		  gens := [];
	  else
	          gens := [hom<L -> L | L.1^q>];
	  end if;
  end if;
  roots := [r[1]: r in Roots(f, L)];
  Sort(~roots);
  ActionOnTheta := [];
  signs := [];
  Lt := [6,9,12,15,16];
  U:={1,3,5};
  eta := EtaFunction0(2);
  sigmas:=[];
  for g in gens do
    if p eq 0 then
            sigma_g := mM * phi(g) * mM^(-1);
    else
            sigma_g := g;
    end if;
    Append(~sigmas, sigma_g);    
    g_roots := [sigma_g(r): r in roots];
    indices := [1..#roots];
    ParallelSort(~g_roots, ~indices);
    if #indices eq 5 then
      Append(~indices, 6);
    end if;
    
    Append(~ActionOnTheta, [0,0,0,0,0]);
    Append(~signs, [0,0,0,0,0]);

    for i in [1..5] do
      for j in [i+1..5] do
        T := {i,j, 6};	
        S := T sdiff U;
        cha := EtaValue(eta, Setseq(S));
        cha_i := TCharToIndex(cha);
        if cha_i in Lt then
          pos := Position(Lt, cha_i); 
          Tnew := {indices[i], indices[j], indices[6]};
          Snew := Tnew sdiff U;
          signswitch:=(-1)^( #(S meet U)- #(Snew meet U));
          chanew := EtaValue(eta, Setseq(Snew));
          cha_inew := TCharToIndex(chanew);
          ActionOnTheta[#ActionOnTheta][pos]:= cha_inew;
          signs[#ActionOnTheta][pos]:= signswitch;
        end if;
      end for;
    end for;
    
  end for;
  return ActionOnTheta,signs,  sigmas;
end function;

intrinsic IgusaCoordinates(K::Fld) -> Any
  {}
  M := ZeroMatrix(K, 5, 16);
  M[1, 8]:=1;
  M[1, 12]:=-1;
  M[1, 15]:=1;
  M[1, 9]:=-1;
  M[2, 16]:=1;
  M[2, 1]:=-1;
  M[2, 6]:=-1;
  M[2, 12]:=-1;
  M[3, 6]:=1;
  M[3, 2]:=-1;
  M[3, 15]:=-1;
  M[3, 3]:=1;
  M[4, 4]:=1;
  M[4, 16]:=-1;
  M[4, 9]:=1;
  M[4, 3]:=-1;
  M[4, 3]:=1;
  M[5, 4]:=1;
  M[5, 8]:=-1;
  M[5, 1]:=1;
  M[5, 2]:=-1;
  Ech := EchelonForm(M);
  pivots := [Min([i: i in [1..16]| Ech[j,i] ne 0  ]): j in [1..5]];
  even := [ 1, 2, 3, 4, 6, 8, 9, 12, 15, 16 ];
  R16 := PolynomialRing(K, 16);
  t8 := &+[R16.i^2: i in even];
  t16 := &+[R16.i^4: i in even];
  equ := t8^2-4*t16;

  subs := Setseq(MonomialsOfDegree(R16,1));

  for i in [1..5] do
  piv:= pivots[i];
  subs[piv]:= -&+[Ech[i, j]*R16.j: j in [piv+1..16]];
  end for;
  equnew := Evaluate(equ, subs);


  R4 := PolynomialRing(K,5);
  subsfin := [R4!0: i in [1..16]];
  subsfin[6] := R4.1;
  subsfin[9] := R4.2;
  subsfin[12] := R4.3;
  subsfin[15] := R4.4;
  subsfin[16] := R4.5;
  equfin:=Evaluate(equnew, subsfin);

  evals:=[Evaluate(sub, subsfin): sub in subs];
  vecs := [Vector([MonomialCoefficient(ev, R4.i): i in [1..5]]): ev in evals];
  t6 := R16.1;
  t7 := R16.6;
  t8 := R16.2;
  t9 := R16.3;
  t10 := R16.4;
  x1 := Evaluate(t6 - t7 + t9 + 2*t10, evals);
  x2 := Evaluate(-2* t6 - t7 + t9 - t10, evals);
  x3 := Evaluate(t6 - t7 - 2*t9 - t10, evals);
  x4 := Evaluate(t6 + 2*t7 - 3*t8 + t9 + 2*t10, evals);
  x5 := Evaluate(-2*t6 - t7 + 3*t8 - 2*t9 - t10, evals);
  x6 := Evaluate(t6 + 2*t7 + t9 - t10, evals);
  xi:=[x1,x2,x3,x4,x5,x6];

  return vecs, equfin, xi;
  end intrinsic;

  function ComputeCocycle(f)
      K := BaseRing(Parent(f));
      ActionOnTheta, signs, sigmas := ComputeGaloisAction(f);
      vecs, Igusa := IgusaCoordinates(K);
      return [Transpose(Matrix([signs[j][i] * vecs[ActionOnTheta[j][i]]: i in [1..5]]  )): j in [1..#ActionOnTheta]], sigmas, Igusa;
  end function;

  intrinsic IgusaTwist(f::RngUPolElt) -> Any
    {}
      K := BaseRing(Parent(f));
      coc, sigmas, Igusa := ComputeCocycle(f);
      if #sigmas eq 0 then
        U0 := IdentityMatrix(K,5);
        Igtest := Igusa;
 	L := K;
      else
        L := Domain(sigmas[1]);
        K := BaseRing(Parent(f));
        p := Characteristic(K);
        coc := [ChangeRing(co, L): co in coc];
        if p eq 0 then
                cocV := [BlockMatrix(5,5,[RepresentationMatrix(a, K): a in Eltseq(co)]): co in coc];
        else
                _, mat := MatrixAlgebra(L, K);
                cocV := [BlockMatrix(5,5,[mat(a): a in Eltseq(co)]): co in coc];
        end if;
        V, m := KSpace(L, K);
        n:= Dimension(V);
        sigmasV:= [m^(-1) * sig * m: sig in sigmas];
        Matsigmas := [Transpose(Matrix( [ sig(V.i): i in [1..n]]  )): sig in sigmasV];
        Matsigmas5 := [DiagonalJoin([mat: i in [1..5]]): mat in Matsigmas];
        IntKer := &meet[Kernel(Transpose(Matsigmas5[i])^(-1)-cocV[i]): i in [1..#coc]];
        U0 := Matrix([ [Inverse(m)(V!(Eltseq(b)[n*i+1..n*(i+1)])): i in [0..4]]:      b in Basis(IntKer)]);

        R4 := PolynomialRing(L, 5);
        UR4:=ChangeRing(U0, R4);
        Igusa:=Evaluate(Igusa, [R4.i: i in [1..5]]);
        X := Matrix(R4,5,1, [R4.i: i in [1..5]]);
        Igtest:=Evaluate(Igusa, Eltseq(Transpose(UR4)*X));
      end if;

      roots := Roots(f, L);
      roots := [ro[1]: ro in roots];
      Sort(~roots);
      U:={1,3,5};
      eta := EtaFunction0(2);

    theta4 := [L!0: i in [1..16]];

    for i in [1..5] do
      for j in [i+1..5] do
        T := {i,j, 6};
        Tcomp := {1..6} diff T;
        S := T sdiff U;
        sign := (-1)^#(S meet U);
        cha := EtaValue(eta, Setseq(S));
        theta4[TCharToIndex(cha)] := sign * &*[ &*[(roots[nu]-roots[mu])^(-1): nu in T]: mu in Tcomp];
      end for;
    end for;
      Point :=[ theta4[j]: j in [6,9,12,15,16]];
      PointOnTwist:=[el: el in Eltseq(Vector(Point)*U0^(-1))];
      PointOnTwist:=[K!(el/PointOnTwist[1]): el in PointOnTwist];
      return Igtest, PointOnTwist, U0;
end intrinsic;

function TakaseQuotientSq(theta4, eta, k, l, m)
  g := #Eltseq(eta[1]) div 2;
  U := {1,3,5};
  Bm := { 1..2*g+1 }; L := [ bp : bp in (Bm diff { k, l, m }) ];
  V := { L[i]: i in [1..g-1] }; W := { L[i]: i in [g..2*(g-1)] };
  num1 := theta4[TCharToIndex(EtaValue(eta, U sdiff (V join { k, l })))];
  num2 := theta4[TCharToIndex(EtaValue(eta, U sdiff (W join { k, l })))];
  den1 := theta4[TCharToIndex( EtaValue(eta, U sdiff (V join { k, m })))];
  den2 := theta4[TCharToIndex(EtaValue(eta, U sdiff (W join { k, m })))];
  return (num1*num2)/(den1*den2);
end function;

intrinsic HyperellipticCurveFromTheta4(theta4::SeqEnum) -> Any
  {}
  g:=2;
  L:=Parent(theta4[1]);
  eta:=EtaFunction0(2);
  //ros := [ Sqrt(TakaseQuotientSq(theta4, eta, 1, l, 2)) : l in [3..2*g+1] ];
  ros := [];
  for l in [3..2*g+1] do
    bool, s := IsSquare(TakaseQuotientSq(theta4, eta, 1, l, 2));
    assert bool;
    Append(~ros,s);
  end for;
  U:={1,3,5};
  R<X>:=PolynomialRing(L);
  for loo in [0..7] do
    bin:=Intseq(loo,2,3);
    rostest:= [ros[water]*(-1)^bin[water]: water in [1..3]];
    repeat wei := L!Random(100);
    until &and[((wei-1)*ro-wei) ne 0: ro in rostest] and not wei in [L|0, 1];
    rootstest := [0,1] cat [-ro/((wei-1)*ro-wei): ro in rostest] cat [-1/(wei-1)];
    if #Seqset(rootstest) ne #rootstest then
      continue;
    end if; 
    ftest:=X*(X-1)*&*[(X-L!ro): ro in rostest];
    //print IgusaInvariantsEqual(IgusaInvariants(f), IgusaInvariants(ftest));
    theta4test:= [L!0: i in [1..16]];
    for i in [1..5] do
      for j in [i+1..5] do
        T := {i,j, 6};
        Tcomp := {1..6} diff T;
        S := T sdiff U;
        sign := (-1)^#(S meet U);
        cha := EtaValue(eta, Setseq(S));
        theta4test[TCharToIndex(cha)] := sign * &*[ &*[(rootstest[nu]-rootstest[mu])^(-1): nu in T]: mu in Tcomp];
      end for;
    end for;
    if &and[theta4test[i]*theta4[1] eq theta4[i]*theta4test[1]: i in [1..16]] then
      return ftest, rostest;
    end if;
  end for;
end intrinsic;

