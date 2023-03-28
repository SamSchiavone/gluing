AttachSpec("spec");
AttachSpec("~/github/CHIMP/CHIMP.spec");
Attach("~/github/gluing/magma/analytic/reconstruction.m");
SetDebugOnError(true);
//prec := 40;
prec := 15;
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
g := Genus(S);
Pi := BigPeriodMatrix(S);
Pi1, Pi2 := SplitBigPeriodMatrix(Pi);
tau := SmallPeriodMatrix(S);
//assert Max([Abs(Eltseq(Pi1^-1*Pi2)[i] - Eltseq(tau)[i]) : i in [1..#Eltseq(tau)]]) lt 10^(-prec/3);

chars := OddThetaCharacteristics(g);

// make Steiner complex
v := Vector([GF(2) | 0, 0, 0, 0, 0, 0, 1, 0 ]);
steiner := [];
for c1 in chars do
  for c2 in chars do
    cc1 := [GF(2)!el : el in c1[1] cat c1[2]];
    cc2 := [GF(2)!el : el in c2[1] cat c2[2]];
      if Vector(cc1) - Vector(cc2) eq v then
        Append(~steiner, [c1,c2]);
      end if;
  end for;
end for;

tritangents := [];
//for char in chars[1..4] do
//for char in chars do
for s in steiner[1..7] do
  new := [];
  for char in s do
    Append(~new, TritangentPlane(Pi, char));
  end for;
  Append(~tritangents, new);
end for;

chars_even := EvenThetaCharacteristics(3);
eps := chars_even[1];
eta_sqs := AssociativeArray();
for delta in chars_even do
  print delta;
  eps_new1 := [[0] cat el : el in eps];
  eps_new2 := [[0] cat eps[1], [1] cat eps[2]];
  delta_new1 := [[0] cat el : el in delta];
  delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
  eta_sqs[[eps, delta]] := (Theta([CC!0 : i in [1..g]], tau : char := eps_new1)*Theta([CC!0 : i in [1..g]], tau : char := eps_new2))/(Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2)); 
end for;

thetas := [];
for i := 1 to 63 do
  s := Intseq(i,2,6);
  s := Reverse(s);
  delta := [s[1..3], s[4..6]];
  if IsDefined(eta_sqs, [eps,delta]) then
    Append(~thetas, Sqrt(eta_sqs[[eps,delta]])); // sign???
  else
    Append(~thetas, 0);
  end if;
end for;

// copy-pasta from curve_reconstruction riemann.m
function ModuliFromTheta(thetas);
  I:=Parent(thetas[1]).1;
  a1:=I*thetas[33]*thetas[5]/(thetas[40]*thetas[12]);
  a2:=I*thetas[21]*thetas[49]/(thetas[28]*thetas[56]);
  a3:=I*thetas[7]*thetas[35]/(thetas[14]*thetas[42]);
  ap1:=I*thetas[5]*thetas[54]/(thetas[27]*thetas[40]);
  ap2:=I*thetas[49]*thetas[2]/(thetas[47]*thetas[28]);
  ap3:=I*thetas[35]*thetas[16]/(thetas[61]*thetas[14]);
  as1:=-thetas[54]*thetas[33]/(thetas[12]*thetas[27]);
  as2:=thetas[2]*thetas[21]/(thetas[56]*thetas[47]);
  as3:=thetas[16]*thetas[7]/(thetas[42]*thetas[61]);
  return [a1,a2,a3,ap1,ap2,ap3,as1,as2,as3];
end function;

ModuliFromTheta(thetas);

CC4<X,Y,Z,W> := PolynomialRing(CC,4);
ells := [&+[c[1][i]*CC4.i : i in [1..g] ]: c in tritangents];
ell_ps := [&+[c[2][i]*CC4.i : i in [1..g] ]: c in tritangents];

// PP^3 attempt
Ccan, phi := CanonicalImage(S);
for cs_new in tritangents do
  print cs_new;
  TritangentSanityCheck(Ccan, cs_new);
end for;
