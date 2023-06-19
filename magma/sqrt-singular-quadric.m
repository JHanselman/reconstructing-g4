function SquareRootOnCone(f6)
  //CC<I> := ComplexField();
  //R<z,x,y,w> := PolynomialRing(CC,4);
  R<z,x,y,w> := Parent(f6);
  CC<I> := BaseRing(R);
  RQ<Z,X,Y,W> := PolynomialRing(QQ,4);
  S, mp := quo< RQ | X*Y - Z^2 >;
  mons6 := [mp(el) : el in MonomialsOfDegree(RQ,6) | not IsDivisibleBy(el, Z^2)];
  mons3 := [mp(el) : el in MonomialsOfDegree(RQ,3) | not IsDivisibleBy(el, Z^2)];

  /*
    prods := AssociativeArray();
    for i := 1 to #mons3 do
      m1 := mons3[i];
      for j := i to #mons3 do
        m2 := mons3[j];
        prods[[m1,m2]] := m1*m2;
      end for;
    end for;
  */

  prods := AssociativeArray();
  for i := 1 to #mons3 do
    m1 := mons3[i];
    for j := i to #mons3 do
      m2 := mons3[j];
      if IsDefined(prods, m1*m2) then
        prods[m1*m2] := Append(prods[m1*m2], [m1,m2]);
      else
        prods[m1*m2] := [[m1,m2]];
      end if;
    end for;
  end for;

  // TODO: Fix problems with going from QQ to CC and to and from the quotient ring
  //known := [mp(X^3)];
  known := [mp(W^3)];
  f3 := Sqrt(CC!MonomialCoefficient(f6,R!known[1]^2))*known[1]; // use @@ ?
  while #known ne #mons3 do
    for k->v in prods do
      all_pairs := &cat v;
      unknown := [el : el in all_pairs | not el in known];
      if #unknown eq 1 then
        partner := [el : el in [pair : pair in v | unknown[1] in pair][1] | el ne unknown[1]];
        if partner[1] eq known[1] then
          Append(~known, unknown[1]);
          //printf "now %o known\n", unknown[1];
          //printf "now #known = %o\n", #known;
          Remove(~prods, k);
          c := CC!MonomialCoefficient(f6,R!k);
          for pair in [pair : pair in v | not unknown[1] in pair] do
            if pair[1] eq pair[2] then
              c -:= MonomialCoefficient(f3,pair[1])*MonomialCoefficient(f3,pair[2]);
            else
              c -:= 2*MonomialCoefficient(f3,pair[1])*MonomialCoefficient(f3,pair[2]);
            end if;
          end for;
          c /:= 2*MonomialCoefficient(f3,known[1]);
          f3 +:= c*unknown[1];
        end if;
      end if;
    end for;
  end while;
end function;
