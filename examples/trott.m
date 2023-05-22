AttachSpec("spec");
SetDebugOnError(true);
prec := 300;
R<x,y> := PolynomialRing(QQ,2);
f := 12^2*(x^4+y^4) - 15^2*(x^2+y^2) + 350*x^2*y^2 + 81;
C := Curve(Spec(R), f);
S := RiemannSurface(f);
g := Genus(S);
Pi_big := BigPeriodMatrix(S);
Pi1 := Submatrix(Pi_big,1,1,g,g);
Pi2 := Submatrix(Pi_big,1,g+1,g,g);
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

/*
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
*/

aronhold := [
  [[1,1,1],[1,1,1]],[[0,0,1],[0,1,1]],[[0,1,1],[0,0,1]],[[1,0,1],[1,0,0]],[[1,0,0],[1,0,1]],[[1,1,0],[0,1,0]],[[0,1,0],[1,1,0]], [[0,1,0],[0,1,0]], [[1,0,0],[1,1,0]], [[1,1,0],[1,0,0]]
  ];
q := aronhold;
q := [Vector(GF(2), el[1] cat el[2]) : el in aronhold];
odd_thetas := [Vector(GF(2), el[1] cat el[2]) : el in aronhold[1..4]];

nums := [];
for q_s in [q[1], q[2], q[3]] do
  num := 1;
  for c in [q[5]+q[6]+q_s, q[5]+q[7]+q_s, q[6]+q[7]+q_s] do
    s := Eltseq(c);
    s := [ZZ!el : el in s];
    c_seq := [s[1..g], s[g+1..(2*g)]];
    num *:= Theta([CC!0 : i in [1..g]], tau : char := c_seq);
  end for;
  Append(~nums, num);
end for;

den := 1;
for c in [q[5]+q[6]+q[4], q[5]+q[7]+q[4], q[6]+q[7]+q[4]] do
  s := Eltseq(c);
  s := [ZZ!el : el in s];
  c_seq := [s[1..g], s[g+1..(2*g)]];
  den *:= Theta([CC!0 : i in [1..g]], tau : char := c_seq);
end for;


