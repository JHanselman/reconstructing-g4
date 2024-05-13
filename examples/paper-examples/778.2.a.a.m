
// Make newform 778.2.a.a in Magma, downloaded from the LMFDB on 06 May 2024.
// To make the character of type GrpDrchElt, type "MakeCharacter_778_a();"
// To make the coeffs of the qexp of the newform in the Hecke field type "qexpCoeffs();"
// To make the newform (type ModFrm), type "MakeNewformModFrm_778_2_a_a();".
// This may take a long time!  To see verbose output, uncomment the SetVerbose lines below.
// The precision argument determines an initial guess on how many Fourier coefficients to use.
// This guess is increased enough to uniquely determine the newform.
// To make the Hecke irreducible modular symbols subspace (type ModSym)
// containing the newform, type "MakeNewformModSym_778_2_a_a();".
// This may take a long time!  To see verbose output, uncomment the SetVerbose line below.
// The default sign is -1.  You can change this with the optional parameter "sign".
function ConvertToHeckeField(input: pass_field := false, Kf := [])
    if not pass_field then
        poly := [1, 4, -4, -1, 1];
        Kf := NumberField(Polynomial([elt : elt in poly]));
        AssignNames(~Kf, ["nu"]);
    end if;
    Rf_num := [[1, 0, 0, 0], [0, 1, 0, 0], [-2, 0, 1, 0], [0, -3, 0, 1]];
    Rf_basisdens := [1, 1, 1, 1];
    Rf_basisnums := ChangeUniverse([[z : z in elt] : elt in Rf_num], Kf);
    Rfbasis := [Rf_basisnums[i]/Rf_basisdens[i] : i in [1..Degree(Kf)]];
    inp_vec := Vector(Rfbasis)*ChangeRing(Transpose(Matrix([[elt : elt in row] : row in input])),Kf);
    return Eltseq(inp_vec);
end function;


// To make the character of type GrpDrchElt, type "MakeCharacter_778_a();"
function MakeCharacter_778_a()
    N := 778;
    order := 1;
    char_gens := [391];
    v := [1];
    // chi(gens[i]) = zeta^v[i]
    assert UnitGenerators(DirichletGroup(N)) eq char_gens;
    F := CyclotomicField(order);
    chi := DirichletCharacterFromValuesOnUnitGenerators(DirichletGroup(N,F),[F|F.1^e:e in v]);
    return MinimalBaseRingCharacter(chi);
end function;

function MakeCharacter_778_a_Hecke(Kf)
    return MakeCharacter_778_a();
end function;


function ExtendMultiplicatively(weight, aps, character)
    prec := NextPrime(NthPrime(#aps)) - 1; // we will able to figure out a_0 ... a_prec
    primes := PrimesUpTo(prec);
    prime_powers := primes;
    assert #primes eq #aps;
    log_prec := Floor(Log(prec)/Log(2)); // prec < 2^(log_prec+1)
    F := Universe(aps);
    FXY<X, Y> := PolynomialRing(F, 2);
    // 1/(1 - a_p T + p^(weight - 1) * char(p) T^2) = 1 + a_p T + a_{p^2} T^2 + ...
    R<T> := PowerSeriesRing(FXY : Precision := log_prec + 1);
    recursion := Coefficients(1/(1 - X*T + Y*T^2));
    coeffs := [F!0: i in [1..(prec+1)]];
    coeffs[1] := 1; //a_1
    for i := 1 to #primes do
        p := primes[i];
        coeffs[p] := aps[i];
        b := p^(weight - 1) * F!character(p);
        r := 2;
        p_power := p * p;
        //deals with powers of p
        while p_power le prec do
            Append(~prime_powers, p_power);
            coeffs[p_power] := Evaluate(recursion[r + 1], [aps[i], b]);
            p_power *:= p;
            r +:= 1;
        end while;    
    end for;
    Sort(~prime_powers);
    for pp in prime_powers do
        for k := 1 to Floor(prec/pp) do
            if GCD(k, pp) eq 1 then
                coeffs[pp*k] := coeffs[pp]*coeffs[k];
            end if;
        end for;
    end for;
    return coeffs;
end function;


function qexpCoeffs()
    // To make the coeffs of the qexp of the newform in the Hecke field type "qexpCoeffs();"
    weight := 2;
    raw_aps := [[1, 0, 0, 0], [-1, 0, -1, -1], [-2, 0, 1, 0], [0, -1, 0, 2], [-1, 2, 1, 0], [-2, -2, 0, 1], [-4, 2, -1, -3], [-1, 1, 0, 1], [-3, -1, 1, -3], [0, 1, -2, -1], [0, -3, 0, 3], [-1, 0, -2, -2], [-4, 2, -2, -2], [0, 2, 0, -3], [-6, 0, 2, -2], [1, -2, -2, 2], [1, -1, -2, 2], [1, -3, 2, 5], [1, -1, 1, -3], [0, 4, 0, -3], [0, 1, -1, -1], [2, -1, 4, -1], [4, 1, -1, 7], [5, -3, 6, 10], [-4, 3, 0, -1], [-8, 0, -1, -16], [0, -3, -3, -6], [10, -2, 1, 17], [11, -7, 0, 12], [-3, -2, 0, 0], [0, -2, 3, 6], [-1, -7, -1, 0], [-4, -4, -6, 1], [0, 9, 3, -6], [-4, -8, 3, -3], [-2, 4, -1, -10], [-5, 1, 11, 8], [10, 5, 1, 3], [-14, 1, -1, 1], [-4, 1, 4, 2], [6, -2, 8, 3], [5, 1, -11, -5], [-9, 2, -2, -10], [5, -2, 7, 6], [-1, 4, -1, 1], [3, -3, 2, 2], [8, -7, 2, 2], [-2, 2, 7, 3], [0, -1, -12, -3], [-9, 12, 8, 3], [7, -5, -4, 9], [-5, -5, -3, -4], [-13, 2, 4, -8], [-18, 14, -5, -18], [-5, 6, 3, 11], [-1, 11, -8, -5], [12, 9, -4, 8], [-6, -5, -5, -1], [-12, 0, 3, -4], [-5, 10, -13, -10], [14, -7, -10, -5], [-13, -7, 0, 3], [-9, 3, 17, 15], [-3, -2, -10, -3], [-1, -3, 10, 11], [-5, 7, 9, -8], [13, 6, -2, 5], [-17, -2, -3, -13], [-4, -2, 11, 1], [-13, -10, 9, 1], [8, 1, 5, 16], [-1, -9, -4, -6], [-3, 10, -18, -14], [15, -9, 5, 7], [9, -7, -10, 9], [-10, 12, 13, 9], [1, 0, 0, 0], [-10, 7, 15, 13], [-5, -5, -1, 6], [-5, 1, 6, 3], [-6, -6, -2, -1], [12, 0, -10, 1], [7, 8, 14, -2], [10, 2, -2, -3], [-5, 11, 3, -10], [-13, -11, -3, -3], [-11, 9, -10, -12], [-5, 3, 8, 16], [6, -7, 10, 9], [3, 3, 7, 0], [-18, 6, -14, -21], [2, -3, -2, -3], [-16, 12, 11, 7], [2, 2, -11, 1], [29, -4, 4, 15], [-17, 11, 10, -3], [-13, -8, -12, 9], [14, 7, -2, 7], [2, 8, -1, -6], [-4, 13, 18, -5], [9, -7, -1, 3], [-11, 10, 2, -8], [-24, 9, -6, -24], [26, -7, -7, 0], [13, -20, -11, -6], [8, -10, -7, 15], [-20, -8, 3, -2], [14, -19, 3, 17], [-4, -2, -17, 5], [12, -17, -1, 5], [1, 13, -12, -14], [12, 3, -14, -10], [-4, -15, 3, 13], [1, -7, 20, 3], [7, -5, -6, -4], [10, 14, -8, 9], [29, 9, -2, 10], [14, -4, 3, 35], [-14, -1, -8, -30], [18, 10, -7, -6], [-5, 7, -14, -22], [13, -8, -11, -15], [-2, 14, 14, -14], [-10, 2, -7, -18], [29, -3, 1, 11], [-11, -3, -4, -16], [-5, -12, 1, -10], [20, -13, 15, 29], [18, -16, -13, 6], [6, 6, 0, -4], [19, -10, 4, 31], [15, -6, -12, 11], [35, -7, -3, -4], [27, -12, -11, -2], [-10, 21, -2, -6], [-16, 27, -5, -16], [21, -21, -3, 9], [5, -2, 4, -18], [-30, 9, 10, 8], [5, 5, 15, -5], [6, 11, -14, -14], [2, 0, 9, -2], [23, -19, 7, 26], [-23, 9, -10, -30], [-11, 2, -2, -23], [-13, 2, 17, 18], [8, -10, -9, -3], [-2, -9, -12, -5], [-5, 13, 1, -10], [-18, -4, -1, -22], [-13, 16, 1, -24], [21, -19, -3, 19], [-25, 10, -8, -14], [-27, 7, 9, -21], [6, -12, 18, 12], [10, -15, 0, 8], [3, 3, -17, -5], [12, -9, 10, 7], [-30, 14, -7, -27], [0, 12, -5, -7], [-8, 4, -8, 6], [-10, -5, -19, -24], [2, -7, -12, -4], [17, -16, 10, 16], [-19, 5, -5, -5], [-12, -11, -10, -12], [-22, -1, 5, 9], [-11, -7, -10, -9]];
    aps := ConvertToHeckeField(raw_aps);
    chi := MakeCharacter_778_a_Hecke(Universe(aps));
    return ExtendMultiplicatively(weight, aps, chi);
end function;


// To make the newform (type ModFrm), type "MakeNewformModFrm_778_2_a_a();".
// This may take a long time!  To see verbose output, uncomment the SetVerbose lines below.
// The precision argument determines an initial guess on how many Fourier coefficients to use.
// This guess is increased enough to uniquely determine the newform.
function MakeNewformModFrm_778_2_a_a(:prec:=4)
    chi := MakeCharacter_778_a();
    f_vec := qexpCoeffs();
    Kf := Universe(f_vec);
    // SetVerbose("ModularForms", true);
    // SetVerbose("ModularSymbols", true);
    S := CuspidalSubspace(ModularForms(chi, 2));
    S := BaseChange(S, Kf);
    maxprec := NextPrime(997) - 1;
    while true do
        trunc_vec := Vector(Kf, [0] cat [f_vec[i]: i in [1..prec]]);
        B := Basis(S, prec + 1);
        S_basismat := Matrix([AbsEltseq(g): g in B]);
        if Rank(S_basismat) eq Min(NumberOfRows(S_basismat), NumberOfColumns(S_basismat)) then
            S_basismat := ChangeRing(S_basismat,Kf);
            f_lincom := Solution(S_basismat,trunc_vec);
            f := &+[f_lincom[i]*Basis(S)[i] : i in [1..#Basis(S)]];
            return f;
        end if;
        error if prec eq maxprec, "Unable to distinguish newform within newspace";
        prec := Min(Ceiling(1.25 * prec), maxprec);
    end while;
end function;


// To make the Hecke irreducible modular symbols subspace (type ModSym)
// containing the newform, type "MakeNewformModSym_778_2_a_a();".
// This may take a long time!  To see verbose output, uncomment the SetVerbose line below.
// The default sign is -1.  You can change this with the optional parameter "sign".
function MakeNewformModSym_778_2_a_a( : sign := -1)
    R<x> := PolynomialRing(Rationals());
    chi := MakeCharacter_778_a();
    // SetVerbose("ModularSymbols", true);
    Snew := NewSubspace(CuspidalSubspace(ModularSymbols(chi,2,sign)));
    Vf := Kernel([<3,R![1, -3, -1, 3, 1]>],Snew);
    return Vf;
end function;