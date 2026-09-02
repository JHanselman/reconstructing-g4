/***
 *  Reconstructing genus 4 curves
 *
 *  Jeroen Hanselman
 *  Andreas Pieper
 *  Sam Schiavone
 *
 *  See LICENSE.txt for license details.
 *
 *  Theta functions via FLINT's acb_theta module, called through the small
 *  external C program Cext/acb_theta_cli (no Python / python-flint needed).
 *  Build it first:  cd Cext && make    (or make mac / make server; see
 *  Cext/README.md).  The binary is located as follows, first hit wins:
 *    1. the environment variable THETA_FLINT_BIN,
 *    2. <this package>/../Cext/acb_theta_cli,
 *    3. acb_theta_cli on $PATH.
 *
 *  Conventions (identical to CHIMP's Theta.magma and to FLINT):
 *    * a characteristic is a (2g) x 1 matrix with entries in {0, 1/2}
 *      (top half a, bottom half b), numbered ab = 2^g * a + b with the
 *      first entry of a (resp. b) as most significant bit;
 *    * ThetaFlint(z, tau : ord) returns a sequence of length 2^(2g),
 *      entry ab+1 being the sequence [ theta_ab(z,tau), d/dz_1 theta_ab, ...,
 *      d/dz_g theta_ab, <higher partials up to total order ord> ] -- so with
 *      the default ord := 0 each entry is a one-element sequence and
 *      [ th[1] : th in ThetaFlint(z, tau) ] is the flat list of theta values.
 *      Higher partials (ord >= 2) come in FLINT's tuple order (by total order,
 *      then reverse-lexicographic) and are genuine partial derivatives, not
 *      Taylor coefficients.
 *    * The precision of the result is decided by FLINT's certified error
 *      bounds: the returned complex numbers live in ComplexField(D) with D the
 *      number of digits FLINT guarantees for every returned value, given the
 *      input precision.
 */

declare verbose ThetaFlint, 2;
theta_flint_store := NewStore();

/* ---------- characteristics <-> integers (copied from CHIMP so that this
 *            package does not depend on CHIMP for theta computations) ------ */

function IntegerToHalfCharacteristic(n, g)
  return Matrix(Rationals(), g, 1, [ QQ!e/2 : e in Reverse(Intseq(n, 2, g)) ]) where QQ := Rationals();
end function;

function PairToCharacteristic(pair, g)
  return VerticalJoin(IntegerToHalfCharacteristic(pair[1], g), IntegerToHalfCharacteristic(pair[2], g));
end function;

function HalfCharacteristicToInteger(char)
  assert Ncols(char) eq 1;
  char := ChangeRing(2*char, Integers());
  return SequenceToInteger(Reverse([elt mod 2 : elt in Eltseq(char)]), 2);
end function;

function CharacteristicToPair(char)
  assert Ncols(char) eq 1;
  assert Nrows(char) mod 2 eq 0;
  g := Nrows(char) div 2;
  return <HalfCharacteristicToInteger(Submatrix(char, 1 + g*i, 1, g, 1)) : i in [0,1]>;
end function;

function CharacteristicToInteger(char)
  g := Nrows(char) div 2;
  a, b := Explode(CharacteristicToPair(char));
  return 2^g*a + b;
end function;

function IntegerToCharacteristic(n, g)
  return PairToCharacteristic(<n div 2^g, n mod 2^g>, g);
end function;

intrinsic ThetaCharacteristicToInteger(char::Mtrx) -> RngIntElt
  {The FLINT/CHIMP index 2^g*a + b (bits most significant first) of the
   characteristic char, a (2g) x 1 matrix with entries in \{0, 1/2\}.}
  return CharacteristicToInteger(char);
end intrinsic;

intrinsic ThetaCharacteristicFromInteger(n::RngIntElt, g::RngIntElt) -> Mtrx
  {The characteristic (a (2g) x 1 matrix with entries in \{0, 1/2\}) with
   FLINT/CHIMP index n, 0 <= n < 2^(2g).}
  require 0 le n and n lt 2^(2*g): "n must satisfy 0 <= n < 2^(2g)";
  return IntegerToCharacteristic(n, g);
end intrinsic;

/* ---------- locating the binary ------------------------------------------ */

function FileExists(path)
  try
    _ := Pipe(Sprintf("test -x \"%o\"", path), "");
    return true;
  catch e
    return false;
  end try;
end function;

function AcbThetaBinary()
  if StoreIsDefined(theta_flint_store, "bin") then
    return StoreGet(theta_flint_store, "bin");
  end if;
  candidates := [];
  env := GetEnv("THETA_FLINT_BIN");
  if env ne "" then
    Append(~candidates, env);
  end if;
  // <repo>/magma/FlintWrapper.m  ->  <repo>/Cext/acb_theta_cli
  filenames := GetFilenames(ThetaCharacteristicToInteger);
  if #filenames ge 1 then
    parts := Split(filenames[1,1], "/");
    if #parts ge 3 then
      repo := "/" cat Join(parts[1..#parts-2], "/");
      Append(~candidates, repo cat "/Cext/acb_theta_cli");
    end if;
  end if;
  try
    onpath := Split(Pipe("command -v acb_theta_cli 2>/dev/null", ""), "\n");
    if #onpath ge 1 and onpath[1] ne "" then
      Append(~candidates, onpath[1]);
    end if;
  catch e
    _ := true;
  end try;
  for c in candidates do
    if FileExists(c) then
      vprintf ThetaFlint, 2: "ThetaFlint: using binary %o\n", c;
      StoreSet(theta_flint_store, "bin", c);
      return c;
    end if;
  end for;
  error Sprintf("ThetaFlint: cannot find the acb_theta_cli binary (tried %o).\n" cat
    "Build it with `cd Cext && make` (or `make mac` / `make server`, see Cext/README.md), " cat
    "or set the environment variable THETA_FLINT_BIN to its path.", candidates);
end function;

/* ---------- Magma -> text ------------------------------------------------ */

// A Magma real is a bare floating-point number; we hand it to FLINT as a ball
// whose radius reflects the working precision, so that FLINT's certified
// output precision honestly accounts for the input precision (same choice as
// CHIMP's python-flint wrapper).
function RealToArbString(x, digits)
  RR := RealField(digits);
  x := RR!x;
  rad := Max(Abs(x), RR!1) * RR!10^(1-digits);
  return ReplaceCharacter(Sprintf("%o +/- %.*o", x, 3, rad), "E", "e");
end function;

function ComplexToArbLines(x, digits)
  return RealToArbString(Real(x), digits) cat "\n" cat RealToArbString(Imaginary(x), digits) cat "\n";
end function;

function WriteInput(z, tau, digits, ord)
  g := Nrows(tau);
  bits := Ceiling(digits*Log(2, 10)) + 64;
  s := Sprintf("%o %o %o\n", g, bits, ord);
  for i in [1..g] do
    s cat:= ComplexToArbLines(z[i,1], digits);
  end for;
  for i in [1..g], j in [1..g] do
    s cat:= ComplexToArbLines(tau[i,j], digits);
  end for;
  return s;
end function;

function TmpName()
  return Split(Pipe("mktemp", ""), "\n")[1];
end function;

/* ---------- the call ----------------------------------------------------- */

function CallAcbTheta(z, tau, ord)
  CC := BaseRing(Parent(tau));
  digits := Precision(CC);
  try
    digits := Max(digits, Precision(BaseRing(Parent(z))));
  catch e
    _ := true;
  end try;
  z := ChangeRing(z, CC);
  bin := AcbThetaBinary();
  infile := TmpName();
  outfile := TmpName();
  Write(infile, WriteInput(z, tau, digits, ord) : Overwrite := true);
  cmd := Sprintf("\"%o\" \"%o\" \"%o\"", bin, infile, outfile);
  vprintf ThetaFlint, 2: "ThetaFlint: %o\n", cmd;
  vprintf ThetaFlint: "ThetaFlint: calling acb_theta_cli (g = %o, %o digits, ord = %o)...", Nrows(tau), digits, ord;
  vtime ThetaFlint:
  ret := System(cmd);
  if ret ne 0 then
    System(Sprintf("rm -f \"%o\" \"%o\"", infile, outfile));
    error Sprintf("ThetaFlint: acb_theta_cli failed with status %o (see stderr above)", ret);
  end if;
  out := Read(outfile);
  System(Sprintf("rm -f \"%o\" \"%o\"", infile, outfile));
  return eval out;
end function;

/* ---------- cache -------------------------------------------------------- */

intrinsic CacheClearThetaFlint()
  {Clear the cache of theta values computed by ThetaFlint.}
  StoreSet(theta_flint_store, "cache", AssociativeArray());
end intrinsic;

function CacheKey(z, tau, ord)
  try
    precz := Precision(BaseRing(Parent(z)));
  catch e
    precz := 0;
  end try;
  return <Nrows(tau), Precision(BaseRing(Parent(tau))), precz, ord, Eltseq(tau), Eltseq(z)>;
end function;

function GetCached(z, tau, ord)
  if not StoreIsDefined(theta_flint_store, "cache") then
    return false, _;
  end if;
  cache := StoreGet(theta_flint_store, "cache");
  return IsDefined(cache, CacheKey(z, tau, ord));
end function;

procedure SetCached(z, tau, ord, v)
  if StoreIsDefined(theta_flint_store, "cache") then
    cache := StoreGet(theta_flint_store, "cache");
  else
    cache := AssociativeArray();
  end if;
  cache[CacheKey(z, tau, ord)] := v;
  StoreSet(theta_flint_store, "cache", cache);
end procedure;

/* ---------- public intrinsics ------------------------------------------- */

intrinsic ThetaFlint(z::Mtrx, tau::Mtrx[FldCom] : ord := 0) -> SeqEnum
  {Theta functions theta_ab(z, tau) for all 2^(2g) characteristics, computed
   rigorously with FLINT's acb_theta.  Here z is a g x 1 matrix and tau a
   symmetric g x g complex matrix with positive definite imaginary part.
   Returns a sequence indexed by ab + 1 = 2^g*a + b + 1 (see
   ThetaCharacteristicToInteger); entry ab + 1 is the sequence
   [theta_ab(z,tau), d/dz_1 theta_ab, ..., d/dz_g theta_ab, ...] of all
   partial derivatives in z of total order at most ord (default 0).}
  g := Nrows(tau);
  require Ncols(tau) eq g: "tau must be square.";
  require Nrows(z) eq g and Ncols(z) eq 1: "z must be a g x 1 matrix.";
  require ord ge 0: "ord must be a nonnegative integer.";
  // symmetrise (Magma's tau is only symmetric to working precision)
  tau := (tau + Transpose(tau))/2;
  found, v := GetCached(z, tau, ord);
  if found then
    return v;
  end if;
  // check the imaginary part is positive definite (at reduced precision)
  Y := Matrix(g, g, [ Im(c) : c in Eltseq(tau) ]);
  RRsmall := RealField(Min(15, Precision(BaseRing(tau))));
  require IsPositiveDefinite(ChangeRing(Y, RRsmall)):
    "tau must have positive definite imaginary part.";
  v := CallAcbTheta(z, tau, ord);
  SetCached(z, tau, ord, v);
  return v;
end intrinsic;

intrinsic ThetaFlint(char::Mtrx, z::Mtrx, tau::Mtrx[FldCom] : ord := 0) -> SeqEnum
  {Theta function with characteristic char (a (2g) x 1 matrix with entries in
   \{0, 1/2\}) at z (g x 1) and tau, computed with FLINT's acb_theta.  Returns
   the sequence [theta(z,tau), d/dz_1 theta, ..., d/dz_g theta, ...] of all
   partial derivatives in z of total order at most ord (default 0), so
   ThetaFlint(char, z, tau)[1] is the theta value itself.}
  g := Nrows(tau);
  require Nrows(char) eq 2*g and Ncols(char) eq 1: "char must be a (2g) x 1 matrix.";
  require forall{ c : c in Eltseq(char) | c in {0, 1/2} }: "char must have entries in {0, 1/2}.";
  n := CharacteristicToInteger(char);
  return ThetaFlint(z, tau : ord := ord)[n + 1];
end intrinsic;

intrinsic ThetaFlint(char::SeqEnum, z::Mtrx, tau::Mtrx[FldCom] : ord := 0) -> SeqEnum
  {Same, with the characteristic given as [a, b] with a, b sequences of
   length g with entries in \{0, 1\} (the convention of IndexToTChar), or as a
   flat sequence of length 2g with entries in \{0, 1\} or \{0, 1/2\}.}
  g := Nrows(tau);
  if #char eq 2 and Type(char[1]) eq SeqEnum then
    char := char[1] cat char[2];
  end if;
  require #char eq 2*g: "char must have 2g entries.";
  ch := Matrix(Rationals(), 2*g, 1, [ (Rationals()!c in {0, 1}) select Rationals()!c/2 else Rationals()!c : c in char ]);
  return ThetaFlint(ch, z, tau : ord := ord);
end intrinsic;

intrinsic ThetaConstantsFlint(tau::Mtrx[FldCom]) -> SeqEnum
  {The 2^(2g) theta constants theta_ab(0, tau), as a flat sequence indexed by
   ab + 1, computed with FLINT's acb_theta.}
  g := Nrows(tau);
  CC := BaseRing(Parent(tau));
  return [ th[1] : th in ThetaFlint(ZeroMatrix(CC, g, 1), tau) ];
end intrinsic;
