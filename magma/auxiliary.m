/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 */



//function SplitBigPeriodMatrix(Pi)
intrinsic SplitBigPeriodMatrix(Pi::ModMatFldElt) -> Any
  {}
  g := Nrows(Pi);
  Pi1 := Submatrix(Pi,1,1,g,g);
  Pi2 := Submatrix(Pi,1,g+1,g,g);
  return Pi1, Pi2;
end intrinsic;
//end function;

//For converting characteristics to integers and back
intrinsic TCharToIndex(chara::SeqEnum) -> RngIntElt
{}
  ZZ:=Integers();
  eps := Reverse(chara[1]);
  delta := Reverse(chara[2]);
  g := #eps;
  assert #delta eq g;
  
  sum := &+[(ZZ!eps[i] mod 2)*2^(g + i-1) : i in [1..g]] + &+[(ZZ!delta[i] mod 2)*2^(i-1) : i in [1..g]];
  if sum ne 0 then
    return sum;
  else
    return 2^(2*g);
  end if;
end intrinsic;

intrinsic IndexToTChar(index::RngIntElt, g::RngIntElt) -> SeqEnum
{}
  s := Intseq(index mod 2^(2*g), 2, 2*g);
  s := Reverse(s);
  return [s[1..g], s[g+1..2*g]];
end intrinsic;


intrinsic TritangentPlane(Pi::ModMatFldElt, char::SeqEnum) -> SeqEnum
  {Given a big period matrix Pi for a genus 4 curve and an odd theta characteristic char, return the corresponding tritangent plane}

  CC := BaseRing(Pi);
  prec := Precision(CC);
  Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
  tau := Pi1^-1*Pi2;
  g := Nrows(tau);
  cs := [];
  for i := 1 to g do
    dz := [0,0,0,0];
    dz[i] := 1;
    Append(~cs, Theta([CC | 0,0,0,0], tau : char := char, dz := [dz], prec := prec));
  end for;
  cs := Eltseq(Matrix(1,g,cs)*(Pi1^-1));
  //cs := [cs[i]/cs[g] : i in [1..g]];
  return cs;
end intrinsic;

intrinsic CanonicalImage(S::RieSrf : diffs := []) -> Crv
  {Given a Riemann surface S, return the canonical image and canonical embedding. If differentials are provided, these are used to define the canonical embedding}

  f := DefiningPolynomial(S);
  K := BaseRing(Parent(f));
  Cplane := Curve(Spec(Parent(f)), f);
  // make differentials if not provided
  if #diffs eq 0 then
    basis, M := Explode(HolomorphicDifferentials(S)); // need to deal with superelliptic vs generic case...
    diffs := [];
    //T<X,Y> := FunctionField(QQ,2); // QQ or RationalsAsNumberField()?
    T<X,Y> := FunctionField(K,2); // QQ or RationalsAsNumberField()?
    for j := 1 to Ncols(M) do
      pows := Rows(Transpose(M))[j];
      Append(~diffs, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
    end for;
  end if;
  // make canonical embedding
  phi := map< Cplane -> ProjectiveSpace(K,3) | Eltseq(diffs)>;
  Ccan := Image(phi);
  return Ccan, phi;
end intrinsic;

// copied from /Applications/Magma/package/Geometry/RieSrf/riesrfclass.m, line 661
function LexOrdering(x,y)
  prec := Precision(Parent(x));
  eps := 10^-(prec+1);
  if Abs(Re(x)-Re(y)) gt eps then
    if Re(x) lt Re(y) then 
      return -1;
    elif Re(x) gt Re(y) then 
      return 1;
    end if;
  elif Abs(Im(x)-Im(y)) gt eps then
    if Im(x) lt Im(y) then 
      return -1;
    elif Im(x) gt Im(y) then 
      return 1; 
    end if; 
  else
    return 0;
  end if;
end function;

function MatchRootPairs(roots)
  CC := Parent(roots[1]);
  prec := Precision(CC);
  pairs := [];
  for root in roots do
    new_bool := true;
    for i := 1 to #pairs do
      //printf "i = %o\n", i;
      //printf "difference = %o\n", Abs(root - pairs[i][1]);
      //printf "error tolerance = %o\n", 10^(-prec/3);
      if Abs(root - pairs[i][1]) lt 10^(-prec/3) then
        Append(~pairs[i], root);
        new_bool := false;
        //printf "new_bool = %o\n", new_bool;
      end if;
    end for;
    if new_bool then
      Append(~pairs, [root]);
    end if;
    //printf "now pairs = %o\n", pairs;
  end for;
  return pairs;
end function;

intrinsic TritangentSanityCheck(Ccan::Crv, cs_new::SeqEnum) -> BoolElt
  {Given a canonically embedded curve Ccan and the coefficients of a tritangent plane, check that the plane is indeed tritangent to precision}

  g := Dimension(Ambient(Ccan)) + 1;
  require g eq 4: "Curve must have genus 4";
  eqns0 := DefiningEquations(Ccan);
  eqns := [el : el in eqns0 | Degree(el) in [2,3]];
  J := ideal< Parent(eqns[1]) | eqns >;
  for h in [el : el in eqns0 | not el in eqns] do
    assert h in J; 
  end for;
  CC := Parent(cs_new[1]);
  prec := Precision(CC);
  CC4 := PolynomialRing(CC,g);
  xx := (1/cs_new[1])*&+[-cs_new[i]*CC4.i : i in [2..g]];
  //xx := Evaluate(xx, [CC4.1,CC4.2,CC4.3,1]);
  xx := Evaluate(xx, [CC4.i : i in [1..g-1]] cat [CC4!1]);
  //evals := [Evaluate(el, [xx,CC4.2,CC4.3,1]) : el in eqns];
  evals := [Evaluate(el, [xx] cat [CC4.i : i in [2..g-1]] cat [CC4!1]) : el in eqns];
  r := Resultant(evals[1],evals[2],CC4.3);
  //r;
  _<t> := PolynomialRing(CC);
  r := Evaluate(r,[0,t,0,0]);
  printf "resultant = %o\n", r;
  roots0 :=  Roots(r,CC);
  roots := [];
  for pair in roots0 do
    for i := 1 to pair[2] do
      Append(~roots, pair[1]);
    end for;
  end for;
  assert #roots eq 6;
  Sort(~roots, LexOrdering);
  printf "roots = %o\n", roots;
  // want 3 pairs
  pairs := MatchRootPairs(roots);
  assert #pairs eq 3;
  for pair in pairs do
    assert #pair eq 2;
  end for;
  return true;
end intrinsic;
