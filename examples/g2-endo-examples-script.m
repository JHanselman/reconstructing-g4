AttachSpec("~/github/CHIMP/CHIMP.spec");

SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 300;

F := RationalsExtra(prec);
//P3<x,y,z,w> := ProjectiveSpace(F, 3);
R<x> := PolynomialRing(F);

// TODO: Compare dim of geometric endomorphism algebra with dimension given by Sato-Tate group
file := Open("filtered2.txt","r");
eof := false;
while not eof do
  line := Gets(file);
  if IsEof(line) then
    eof := true;
    break;
  end if;

  line := "cs_list := "*line*"; return cs_list;";
  cs_list := eval line;
  f, h := Explode([R!el : el in cs_list]);
  X := HyperellipticCurve(f,h);

  /*
    f1 := x^2 + 28/5*x*y + 18/5*y^2 - 224/5*z^2 + 102/5*z*w - 33/10*w^2;
    f2 := -27/226*x*y^2 - 205/1808*x*z^2 + 15/113*x*z*w - 185/3616*x*w^2 - 39/452*y^3 + y*z^2 - 271/452*y*z*w + 17/452*y*w^2;
    X := Curve(P3, [f1, f2]);
  */

  print "";
  print "Curve:";
  print X;

  time P := PeriodMatrix(X);

  try 
    print "";
    print "Endomorphism lattice:";
    time lat := HeuristicEndomorphismLattice(X);
    print lat;

    print "";
    print "Endomorphisms over base:";
    print HeuristicEndomorphismDescription(X);

    print "";
    print "Endomorphisms over closure:";
    descript := HeuristicEndomorphismDescription(X : Geometric := true);
    print descript;
    if not "M_2 (QQ)" in descript then
      print "Non-split example found!";
      Write("interesting-g2.txt", Sprintf("%m",X));
    end if;

    print "";
    print "Decomposition:";
    dec := HeuristicDecomposition(X);
    print dec;
  catch e
    Write("g2-errors.txt", Sprintf("%m",X));
    Write("g2-errors.txt", Sprintf("%m",e`Object));
  end try;
end while;
