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
  h:= hom< RQ -> R | z, x, y, w>;
  f3 := Sqrt(MonomialCoefficient(f6,h(known[1]^2)))*h(known[1]); // use @@ ?
  while #known ne #mons3 do
    found := false;
    for k->v in prods do
      all_pairs := &cat v;
      print v;
      unknown := {el : el in all_pairs | not el in known};
      if #unknown eq 1 then
        print unknown;
        unknown_el := Random(unknown);
        partner := [el : el in [pair : pair in v | unknown_el in pair][1] | el ne unknown_el];
        
        //Added this to see if it would help.
        //if IsEmpty(partner) then
         // Append(~known, unknown_el);
         // f3 +:= Sqrt(MonomialCoefficient(f6,h(unknown_el^2)))*h(unknown_el);
         // found := true;
        //end if;
        
        print "partner";
        print partner;
        if (not IsEmpty(partner)) and partner[1] in known[1] then
          Append(~known, unknown_el);
          //printf "now %o known\n", unknown_el;
          //printf "now #known = %o\n", #known;
          Remove(~prods, k);
          c := MonomialCoefficient(f6,h(k));
          for pair in [pair : pair in v | not unknown_el in pair] do
            if pair[1] eq pair[2] then
              c -:= MonomialCoefficient(f3,h(pair[1]))*MonomialCoefficient(f3,h(pair[2]));
            else
              c -:= 2*MonomialCoefficient(f3,h(pair[1]))*MonomialCoefficient(f3,h(pair[2]));
            end if;
          end for;
          c /:= 2*MonomialCoefficient(f3,h(known[1]));
          f3 +:= c*h(unknown_el);
          found := true;
        end if;
      end if;
    end for;
    if not found then
      print "Error";
      break;
    end if;
  end while;
  return f3;
end function;
