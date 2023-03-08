AttachSpec("spec");
SetDebugOnError(true);
prec := 40;
SetDefaultRealFieldPrecision(prec);
RR := RealField(prec);
CC<I> := ComplexField(prec);
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
Pi1 := Submatrix(Pi,1,1,g,g);
Pi2 := Submatrix(Pi,1,g+1,g,g);
tau := SmallPeriodMatrix(S);
assert Max([Abs(Eltseq(Pi1^-1*Pi2)[i] - Eltseq(tau)[i]) : i in [1..#Eltseq(tau)]]) lt 10^-15;
//omegas := S`DFF;
// or omegas := HolomorphicDifferentials(S); ?

function DotProductSeq(v1,v2)
  assert #v1 eq #v2;
  return &+[v1[i]*v2[i] : i in [1..#v2]];
end function;

Pow := CartesianPower(GF(2),g);
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


tritangents := [];
for char in chars[1..4] do
//for char in chars[57..61] do
  cs := [];
  for i := 1 to g do
    dz := [0,0,0,0];
    dz[i] := 1;
    Append(~cs, Theta([CC | 0,0,0,0], tau : char := char, dz := [dz], prec := prec));
  end for;
  cs := Eltseq(Matrix(1,g,cs)*(Pi1^-1));
  cs := [cs[i]/cs[g] : i in [1..g]];
  Append(~tritangents, cs);
end for;

//Append(~cs,Theta([CC | 0,0,0,0], tau : char := [[1,0,0,0],[1,0,0,0]], dz := [dz], prec := prec));

basis := HolomorphicDifferentials(S);
basis, M := Explode(basis);
KCplane<uu,vv> := FunctionField(Cplane);
sections := [Evaluate(el, [uu,vv]) : el in basis];
T<X,Y> := FunctionField(QQ,2);
diffs := [];
for j := 1 to Ncols(M) do
  pows := Rows(Transpose(M))[j];
  Append(~diffs, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
end for;
//diffs := Matrix(1,#sections,sections)*ChangeRing(M,Parent(sections[1]));
phi := map< Cplane -> ProjectiveSpace(QQ,3) | Eltseq(diffs)>;
Ccan := Image(phi);
// Ccan; // LOL look at this thing

CCUV<U,V> := PolynomialRing(CC,2);
//cs_new := Eltseq(Pi1*Matrix(g,1,cs));
//cs_new := Eltseq(Pi1^-1*Matrix(g,1,cs));
//tritangents0 := tritangents;
//tritangents := [Eltseq(Matrix(1,g,cs)*(Pi1^-1)) : cs in tritangents0];

// PP^3 attempt

eqns := DefiningEquations(Ccan);
eqns := eqns[2..3];

CC4 := PolynomialRing(CC,4);
for cs_new in tritangents do
  xx := (1/cs_new[1])*&+[-cs_new[i]*CC4!R.i : i in [2..4]];
  xx := Evaluate(xx, [CC4.1,CC4.2,CC4.3,1]);
  evals := [Evaluate(el, [xx,CC4.2,CC4.3,1]) : el in eqns];
  #evals;
  r := Resultant(evals[1],evals[2],CC4.3);
  //r;
  CCt<t> := PolynomialRing(CC);
  r := Evaluate(r,[0,t,0,0]);
  roots :=  Roots(r,CC);
  roots := [el[1] : el in roots];
  Sort(~roots,S`Ordering);
  print roots;
  print "------";
end for;

// pullback attempt
phi_eqns := DefiningEquations(phi);
f := &+[ cs[i]*CCUV!phi_eqns[i] : i in [1..#cs]];
fplane := CCUV!DefiningEquation(Cplane);
r := Resultant(f,fplane,V);
CCX<X> := PolynomialRing(CC);
r := Evaluate(r,[X,0]);
roots := Roots(r,CC);
[el : el in roots | el[2] eq 2];
roots := [el[1] : el in roots];
Sort(~roots,S`Ordering);

pairs := [];
for i := 1 to #roots do
  for j := i+1 to #roots do
    if Abs(roots[i]-roots[j]) lt 10^-5 then
      Append(~pairs, [roots[i], roots[j]]);
    end if;
  end for;
end for;
#pairs;

//Qs := [Evaluate(el, Eltseq(diffs)) : el in MonomialsOfDegree(R,2)];



// next loop over all odd theta chars

/*
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
phi;
#DefiningEquations(Ccan);
Qs;
CoefficientsAndMonomials(Qs[1]);
phi;
Coefficients(Qs[1]);
[Coefficients(el) : el in Coefficients(Qs[1])];
Coefficients(Qs[1])[1];
Eltseq($1);
uu;
Parent(uu);
Parent(u);
Coefficients(Qs[1])[1];
Evaluate($1,[u,v]);
Evaluate($1,v);
Coefficients($1);
test := Coefficients(Qs[1])[1];
Evaluate(test,v);
Eltseq($1);
v;
Parent(v);
Evaluate(test,v);
Numerator($1);
Coefficients($1);
Qs;
phi;
cs;
[el/cs[1] : el in cs];
[cs[i]*DefiningEquations(phi)[i] : i in [1..#cs]];
[cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs]];
[*cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs]*];
&+$1;
[Type(el) : el in DefiningEquations(phi)];
cs;
[ cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs] ];
[*cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs]*];
s := 0;
for el in [*cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs]*] do
s +:= el;
end for;
s;
s := 0;
for el in [*cs[i]*ChangeRing(DefiningEquations(phi)[i],Parent(cs[1])) : i in [1..#cs]*] do
s := s + el;
end for;
phi_eqns := DefiningEquations(phi);
Parent(phi_eqns[1]) eq Parent(phi_eqns[2]);
CC;
Parent(cs[i]);
Parent(cs[1]);
cs[1];
cs;
phi_eqns[1];
R;
CCUV<U,V> := PolynomialRing(CC,2);
[ cs[i]*CCUV!phi_eqns[i] : i in [1..#cs]];
&+$1;
f := $1;
C;
S;
DefiningEquation(S);
Cplane;
DefiningEquation($1);
fplane := CCUV!$1;
fplane;
f;
Resultant(f,fplane,V);
r := $1;
Factorization(r);
Roots(r,CC);
r;
CCX<X> := PolynomialRing(CC);
r := Evaluate(r, [X,0]);
r;
Roots(r,CC);
roots := $1;
Order(roots);
for i := 1 to #roots do
pairs := [];
for i := 1 to #roots do
for j := i+1 to #roots do
if Abs(roots[i][1]-roots[j][1]) lt 10^-10 then
Append(~pairs, [roots[i][1], roots[j][1]]);
end if;
end for;
end for;
#pairs;
roots := [el[1] : el in roots];
roots;
#roots;
for i := 1 to #roots do
for j := i+1 to #roots do
if Abs(roots[i]-roots[j]) lt 10^-5 then
end if;
end for;
end for;
pairs;
Sort(roots);
Sort;
S;
ListAttributes(Type(S) : Isa := false);
ListAttributes(Type(S));
S`Ordering;
S`Ordering(roots[1],roots[2]);
Sort(roots,S`Ordering);
HolomorphicDifferentials;
HolomorphicDifferentials(S);
S`DFF;
basis;
basis*M;
S`DFF;
diffs;
basis;
M;
S`DFF[1];
omega := $1;
Type(omega);
Eltseq(omega);
new_diffs := [];
for j in Ncols(M) do
pows := Rows(Transpose(M))[j]);
for j in Ncols(M) do
pows := Rows(Transpose(M))[j];
Append(~new_diffs, &*[basis[i]^pows[i] : i in [1..#basis]]);
end for;
for j := 1 to Ncols(M) do
pows := Rows(Transpose(M))[j];
Append(~new_diffs, &*[basis[i]^pows[i] : i in [1..#basis]]);
end for;
basis[1];
basis[2];
basis[2]^-1;
Type(basis[2]);
F;
S;
T<X,Y> := FunctionField(QQ,2);
T;
diffs_new := [];
for j := 1 to Ncols(M) do
Append(~diffs_new, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
end for;
diffs_new;
for j := 1 to Ncols(M) do
pows := Rows(Transpose(M))[j];
Append(~diffs_new, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
end for;
diffs_new := [];
for j := 1 to Ncols(M) do
pows := Rows(Transpose(M))[j];
Append(~diffs_new, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
end for;
diffs_new := [];
for j := 1 to Ncols(M) do
pows := Rows(Transpose(M))[j];
Append(~diffs_new, &*[(T!basis[i])^pows[i] : i in [1..#basis]]);
end for;
diffs_new;
omega;
Parent(omega);
Id($1);
1*d(x);
1*Differential(x);
Differential(1);
Parent(omega).1;
Cplane;
S;
S.1;
Generators(S);
diffs_new;
uu;
vv;
phi := < Cplane -> ProjectiveSpace(QQ,3) | [Evaluate(el,[uu,vv]) : el in diffs_new]>;
phi := < Cplane -> ProjectiveSpace(QQ,3) | [Evaluate(el,[uu,vv]) : el in diffs_new]>;
Evaluate(diffs_new[1],[uu,vv]);
phi := map< Cplane -> ProjectiveSpace(QQ,3) | [Evaluate(el,[uu,vv]) : el in diffs_new]>;
phi;
Ccan := Image(phi);
#DefiningEquations(Ccan);
Ccan;
cs;
CcanCC := ChangeRing(Ccan,CC);
[el/cs[1] : el in cs];
MinimalPolynomial(2.75302387156539483569431256918 - 1.61655716799668399772865791068*I,20);
DefiningEquations(Ccan);
eqns := $1;
eqns_CC := [ChangeRing(el,CC) : el in eqns];
CCUV;
CC4<X,Y,Z,W> := PolynomialRing(CC,4);
eqns_CC := [CC4!el : el in eqns];
eqns_CC;
CcanCC := Curve(Proj(CC4), eqns_CC);
&+[cs[i]*CC4.i : i in [1..#cs];
&+[cs[i]*CC4.i : i in [1..#cs]];
plane := $1;
gens := eqns_CC cat [plane];
I := ideal< CC4 | gens >;
SolveZeroDimIdeal;
DefiningEquations(phi);
phi_eqns := $1;
[ cs[i]*CCUV!phi_eqns[i] : i in [1..#cs]];
&+[ cs[i]*CCUV!phi_eqns[i] : i in [1..#cs]];
f := $1;
fplane := DefiningEquation(Cplane);
fplane;
Resultant(f,CCUV!fplane,V);
*/
