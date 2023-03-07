AttachSpec("spec");
SetDebugOnError(true);
prec := 300;
R<x,y> := PolynomialRing(QQ,2);
f := 12^2*(x^4+y^4) - 15*(x^2+y^2) + 350*x^2*y^2 + 81;
C := Curve(Spec(R), f);
S := RiemannSurface(f);
g := Genus(S);
Pi := BigPeriodMatrix(S);
Pi1 := Submatrix(Pi,1,1,g,g);
Pi2 := Submatrix(Pi,1,g+1,g,g);
tau := SmallPeriodMatrix(S);
assert Max([Abs(Eltseq(Pi1^-1*Pi2)[i] - Eltseq(tau)[i]) : i in [1..#Eltseq(tau)]]) lt 10^-15;

cs := [];
for i := 1 to g do
  dz := [0 : i in [1..g]];
  dz[i] := 1;
  Append(~cs,Theta([CC!0 : i in [1..g]], tau : char := [[1,0,0],[1,0,0]], dz := [dz]));
end for;

