AttachSpec("spec");
SetDebugOnError(true);
prec := 300;
R<x,y> := PolynomialRing(QQ,2);
f := 12^2*(x^4+y^4) - 15^2*(x^2+y^2) + 350*x^2*y^2 + 81;
C := Curve(Spec(R), f);
S := RiemannSurface(f);
g := Genus(S);
Pi := BigPeriodMatrix(S);
Pi1 := Submatrix(Pi,1,1,g,g);
Pi2 := Submatrix(Pi,1,g+1,g,g);
tau := SmallPeriodMatrix(S);
assert Max([Abs(Eltseq(Pi1^-1*Pi2)[i] - Eltseq(tau)[i]) : i in [1..#Eltseq(tau)]]) lt 10^-15;

function DotProductSeq(v1,v2)
  assert #v1 eq #v2;
  return &+[v1[i]*v2[i] : i in [1..#v2]];
end function;

Pow := CartesianPower(GF(2),4);
chars := [];
for v1 in Pow do
  for v2 in Pow do
    w1 := [ZZ!i : i in v1];
    w2 := [ZZ!i : i in v2];
    if DotProductSeq(w1, w2) mod 2 eq 1 then
      Append(~chars, [w1, w2]);
    end if;
  end for;
end for;

//Append(~cs,Theta([CC!0 : i in [1..g]], tau : char := [[1,0,0],[1,0,0]], dz := [dz]));

bitangents := [];
for char in chars do
  cs := [];
  for i := 1 to g do
    dz := [0 : i in [1..g]];
    dz[i] := 1;
    Append(~cs,Theta([CC!0 : i in [1..g]], tau : char := char, dz := [dz]));
  end for;
  //cs := Eltseq(Pi1*Matrix(g,1,cs));
  //cs := Eltseq((Pi1^-1)*Matrix(g,1,cs));
  cs := Eltseq(Matrix(1,g,cs)*(Pi1^-1));
  cs := [cs[i]/cs[3] : i in [1..3]];
  Append(~bitangents, cs);
end for;

