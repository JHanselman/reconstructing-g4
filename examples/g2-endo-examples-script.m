
SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 300;

F := RationalsExtra(prec);
//P3<x,y,z,w> := ProjectiveSpace(F, 3);
R<x> := PolynomialRing(F);

file := Open("filtered2.txt","r");
eof := false;
while not eof do
  line := Gets(file);
  if IsEof(line) then
    eof := true;
    break;
  end if;

  line *:=";";
  cs_list := eval line
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

  print "";
  print "Endomorphism lattice:";
  time lat := HeuristicEndomorphismLattice(X);
  print lat;

  print "";
  print "Endomorphisms over base:";
  print HeuristicEndomorphismDescription(X);

  print "";
  print "Endomorphisms over closure:";
  print HeuristicEndomorphismDescription(X : Geometric := true);

  print "";
  print "Decomposition:";
  dec := HeuristicDecomposition(X);
  print dec;
end while;
