AttachSpec("~/github/quartic_reconstruction/magma/spec");
AttachSpec("~/github/curve_reconstruction/magma/spec");
AttachSpec("~/github/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/github/quartic_isomorphisms/magma/spec");

/* Arithmetic reconstruction */

SetVerbose("QuarticIso", 1);
SetVerbose("QuarticRec", 1);
SetVerbose("Gluing", 1);
//SetVerbose("CurveRec", 1);
SetVerbose("EndoFind", 1);

prec := 500;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f1 := x^5 - x;
f2 := x^5 + 20*x^3 + 36*x;

print "";
print "Can we glue over QQ?";
print ExistsGaloisStableSubgroupFor22(f1, f2);

X1 := HyperellipticCurve(f1);
X2 := HyperellipticCurve(f2);

print "";
print "All arithmetic 2-gluings:";
Ys := AllArithmetic2GluingsCCFor22(X1, X2, F);
Ys;
