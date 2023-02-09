AttachSpec("spec");
SetDebugOnError(true);
prec := 300;
R<x,y,z,w> := PolynomialRing(QQ,4);
Q := x^2 + 2*y^2 + z^2 - w^2;
F := 3*x^3 + y^3 - 7*z^3 + w^3;
C := Curve(Proj(R),[F,Q]);
R2<u,v> := PolynomialRing(QQ,2);
p := Resultant(Q,F,w);
p := Evaluate(p, [R2.1, R2.2, 1, 0]);
S := RiemannSurface(p);
Cplane := Curve(Spec(R2), p);
Genus(S);
Pi := BigPeriodMatrix(S);
tau := SmallPeriodMatrix(S);
//omegas := S`DFF;
// or omegas := HolomorphicDifferentials(S); ?
cs := [];
for i := 1 to 4 do
  dz := [0,0,0,0];
  dz[i] := 1;
  Append(~cs,Theta([CC | 0,0,0,0], tau : char := [[1,0,0,0],[1,0,0,0]], dz := [dz]));
end for;

basis := HolomorphicDifferentials(S);
basis, M := Explode(basis);
KCplane<uu,vv> := FunctionField(Cplane);
sections := [Evaluate(el, [uu,vv]) : el in basis];
diffs := Matrix(1,#sections,sections)*ChangeRing(M,Parent(sections[1]));
phi := map< Cplane -> ProjectiveSpace(QQ,3) | Eltseq(diffs)>;
Ccan := Image(phi);
// Ccan; // LOL look at this thing

Qs := [Evaluate(el, Eltseq(diffs)) : el in MonomialsOfDegree(R,2)];



// next loop over all odd theta chars
