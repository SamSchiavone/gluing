AttachSpec("spec");
Attach("~/github/gluing/magma/analytic/reconstruction.m");
SetDebugOnError(true);
prec := 40;
SetDefaultRealFieldPrecision(prec);
//RR := RealField(prec);
//CC<I> := ComplexField(prec);
R<x,y,z,w> := PolynomialRing(QQ,4);
Q := x^2 + 2*y^2 + z^2 - w^2;
F := 3*x^3 + y^3 - 7*z^3 + w^3;
C := Curve(Proj(R),[F,Q]);
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
assert IsIrreducible(p);
S := RiemannSurface(p : Precision := prec);
Cplane := Curve(Spec(R2), p);
g := Genus(S);
Pi := BigPeriodMatrix(S);
//Pi1, Pi2 := SplitBigPeriodMatrix(Pi)
//tau := SmallPeriodMatrix(S);
//assert Max([Abs(Eltseq(Pi1^-1*Pi2)[i] - Eltseq(tau)[i]) : i in [1..#Eltseq(tau)]]) lt 10^(-prec/3);

chars := OddThetaCharacteristics(g);

tritangents := [];
for char in chars[1..4] do
  Append(~tritangents, TritangentPlane(Pi, char));
end for;

// PP^3 attempt
Ccan, phi := CanonicalImage(S);
for cs_new in tritangents do
  print cs_new;
  TritangentSanityCheck(Ccan, cs_new);
end for;
