/*
  Adapted from Labrande's fast theta code: https://github.com/hlabrand/phdcode
*/

/*
  The duplication formula only gives us theta^2. We can take a square root, but then sign is ambiguous. Functions below use low precision approximations to determine + or -
*/

// given integers n and g, write n as a binary expansion of length 2*g
function IntegerToThetaCharacteristic(n,g)
  QQ := Rationals();
  return (1/2)*Matrix(2*g,1,Reverse([QQ!el : el in Intseq(n,2,2*g)]));
end function;

// ThetaGenus3(12,t) = Theta_12(0,t) & so on
function ThetaMagma(n,t)
  CC := Parent(t[1][1]);
  g := Nrows(t);
  c := IntegerToThetaCharacteristic(n,g);
  //return Theta(c, ZeroMatrix(CC,2*g,1), t);
  return Theta(c, ZeroMatrix(CC,g,1), t);
end function;

// Old version: recompute every time
//     this is wasteful and accounts for a big part of the computation
//     (ex with prec 400 : 35s with, 1.6s with the new version)
// SignTheta(n,t)*Sqrt(thetasquared/theta0squared) = theta/theta0

function SignTheta(n,t)
  Clow := ComplexField(10);
  g := Nrows(t);
  taulow := Matrix(g,g, [ Clow!t[i][j] : i,j in [1..g]]);
  num := ThetaMagma(n, taulow);
  den := ThetaMagma(0, taulow);
  if (Real(num/den) lt 0) then return -1; else return 1; end if;
end function;

function ComputeOneTable(tau)
  Clow := ComplexField(10);
  g := Nrows(tau);
  taulow := Matrix(g, g, [ Clow!tau[i][j] : i,j in [1..g]]);
  mat := ElementToSequence(ZeroMatrix(Integers(),2^(2*g),1));
  den := ThetaMagma(0, taulow);
  for i:= 0 to (2^g-1) do
    num := ThetaMagma(i, taulow);
    if (Real(num/den) lt 0) then mat[i+1] := -1; else mat[i+1] := 1; end if;
  end for;
  return mat;
end function;

// TODO: finish converting this...what is delta?
// used in DiffFiniesOneStep and toInverse
/*
function ComputeTableOfSigns(tau)
  g := Nrows(tau);
  Clow := ComplexField(10);
  taulow := Matrix(g,g, [ Clow!tau[i][j] : i,j in [1..g]]);
  tt := 2*taulow;

  Ig := IdentityMatrix(Clow,g);
  delta11 := Matrix(3,3, [Clow!1,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta22 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta33 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!1]);
  delta12 := Matrix(3,3, [Clow!0,Clow!1,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta13 := Matrix(3,3, [Clow!0,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!0]);
  delta23 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!1,Clow!0]);
  jm1 := (-tt-delta11)*(delta11*tt+delta11-Ig)^(-1);
  jm2 := (-tt-delta22)*(delta22*tt+delta22-Ig)^(-1);
  jm3 := (-tt-delta33)*(delta33*tt+delta33-Ig)^(-1);
  jm12 := (-tt-delta12)*(delta12*tt+delta11+delta22-Ig)^(-1);
  jm13 := (-tt-delta13)*(delta13*tt+delta11+delta33-Ig)^(-1);
  jm23 := (-tt-delta23)*(delta23*tt+delta33+delta22-Ig)^(-1);
  ttinv := -(tau)^(-1);

  return [ ComputeOneTable(taulow), ComputeOneTable(jm1), ComputeOneTable(jm2), ComputeOneTable(jm3), ComputeOneTable(jm12), ComputeOneTable(jm13), ComputeOneTable(jm23), ComputeOneTable(ttinv)];
end function;
*/

// for g=4, 2^g=2^4=16 fundamental theta constants: a,b in ZZ/2^g with a=0
// for g=4, 2^(2*g)=2^8=256 theta constants: a,b in ZZ/2^g

// see Proposition 3.3
function HadamardMatrix(fi, n)
  m := Matrix(2,2, [fi | 1, 1, 1, -1]);
  res := m;
  for i:=2 to n do
    res := TensorProduct(res,m);
  end for;
  return res;
end function;

//TODO: finish this
// Given theta_{0,b}(0,t) compute theta_{a,b}(0,t)
// we still need to find a more genus g way to do this (Makarov = sum ab, not sum (-1)ab ab)
// see Proposition 3.3
function AllDuplication(fund_thetas)
  n := #fund_thetas;
  g := Integers()!Log(2,n);

  ThetaProducts := Matrix(n,n,
    [[fund_thetas[j+1]*fund_thetas[BitwiseXor(j,i)+1] : i in [0..n-1]] : j in [0..n-1]]
		);
  hadam := HadamardMatrix(Parent(fund_thetas[1]), g);

  ThetaProducts := hadam*ThetaProducts;
  ThetaProducts := ElementToSequence(ThetaProducts/n); // why divide?

  return ThetaProducts;
end function;

// duplication formula test
/*
  prec := 60;
  C := ComplexField(prec);
  val := C!0.866;
  tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
  thconstants := [ ThetaMagma(i,tau) : i in [0..7]];
  time allthetas2tau := AllDuplication(thconstants);
  for i := 1 to 64 do
    ComplexField(10)!Abs(allthetas2tau[i]-ThetaMagma(i-1,2*tau)^2);
  end for;
*/

