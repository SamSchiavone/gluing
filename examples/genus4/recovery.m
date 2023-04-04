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
Pi_big := BigPeriodMatrix(S);
Pi1, Pi2 := SplitBigPeriodMatrix(Pi_big);
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
    Append(~new, TritangentPlane(Pi_big, char));
  end for;
  Append(~tritangents, new);
end for;

chars_even := EvenThetaCharacteristics(3);
eps := chars_even[1];
eta_sqs := AssociativeArray(); // now just etas
for delta in chars_even do
  print delta;
  eps_new1 := [[0] cat el : el in eps];
  eps_new2 := [[0] cat eps[1], [1] cat eps[2]];
  delta_new1 := [[0] cat el : el in delta];
  delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
  // formula from Lemma 1 (p. 148) of Farkas
  eta_sqs[[eps, delta]] := Sqrt((Theta([CC!0 : i in [1..g]], tau : char := eps_new1)*Theta([CC!0 : i in [1..g]], tau : char := eps_new2))/(Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2)));  // TODO: fix signs here
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

mods := ModuliFromTheta(thetas);

function RiemannModelFromModuli(mods);
  a1:=mods[1];a2:=mods[2];a3:=mods[3];
  ap1:=mods[4];ap2:=mods[5];ap3:=mods[6];
  as1:=mods[7];as2:=mods[8];as3:=mods[9];
  F:=Parent(a1);
  P<x1,x2,x3>:=PolynomialRing(F,3);
  k:=1;kp:=1;ks:=1;
  M:=Matrix([[1,1,1],[k*a1,k*a2,k*a3],[kp*ap1,kp*ap2,kp*ap3]]);
  Mb:=Matrix([[1,1,1],[1/a1,1/a2,1/a3],[1/ap1,1/ap2,1/ap3]]);
  U:=-Mb^(-1)*M;
  u1:=U[1];
  u2:=U[2];
  u3:=U[3];
  u1:=u1[1]*x1+u1[2]*x2+u1[3]*x3;
  u2:=u2[1]*x1+u2[2]*x2+u2[3]*x3;
  u3:=u3[1]*x1+u3[2]*x2+u3[3]*x3;
  return (x1*u1+x2*u2-x3*u3)^2-4*x1*u1*x2*u2, u1, u2, u3;
end function;


// compute bitangents of genus 3 curve
// formulas from Theorem 6.1.9 (p. 230) of Dolgachev
bitangents := [ [CC | 1, 0, 0], [CC | 0,1,0], [CC | 0,0,1], [CC | 1,1,1]];
bitangents cat:= mods_mat;
F, u0, u1, u2 := RiemannModelFromModuli(mods);
bitangents cat:= [Coefficients(el) : el in [u0, u1, u2]];
CC3<t0,t1,t2> := Parent(u0);
bitangents cat:= [Coefficients(el) : el in [t0+t1+u2, t0+u1+t2, u0+t1+t2]];

// (3)
for i := 1 to 3 do
  new := u0/mods_mat[1,1] + ks[i,1]*(mods_mat[2,i]*t1 + mods_mat[3,i]*t2);
  Append(~bitangents, Coefficients(new));
end for;
// (4)
for i := 1 to 3 do
  new := u1/mods_mat[2,1] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[3,i]*t2);
  Append(~bitangents, Coefficients(new));
end for;
// (5)
for i := 1 to 3 do
  new := u2/mods_mat[3,1] + ks[i,1]*(mods_mat[1,i]*t0 + mods_mat[2,i]*t1);
  Append(~bitangents, Coefficients(new));
end for;
// (6)
for i := 1 to 3 do
  new := u0/(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i]) + u1/(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i]) + u2/(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]);
  Append(~bitangents, Coefficients(new));
end for;
// (7)
for i := 1 to 3 do
  new := u0/(mods_mat[1,i]*(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + u1/(mods_mat[2,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + u2/(mods_mat[3,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
  Append(~bitangents, Coefficients(new));
end for;

fsq_mat := [];
for f in fs do  
  cs := [];
  for m in mons do
    Append(~cs, MonomialCoefficient(f^2,m)); 
  end for;
  Append(~fsq_mat, cs);
end for;
fsq_mat := Matrix(fsq_mat);
K := NumericalKernel(fsq_mat);

// TODO: clean this up
// copy-paste from Yuwei's example Magma-Schottky-Igusa-form file
/*using magma to compute the Schottky form*/
//C := ComplexField(prec);
prec := 30;
C := CC;
char := Matrix(C, 8, 1, [0,0,0,0,0,0,0,0]);
z := Matrix(C, 4, 1, [0,0,0,0]);

m1 := 1/2*Matrix(C, 8, 1, [1,0,1,0,1,0,1,0]);
m2 := 1/2*Matrix(C, 8, 1, [0,0,0,1,1,0,0,0]);
m3 := 1/2*Matrix(C, 8, 1, [0,0,1,1,1,0,1,1]);
n0 := Matrix(C, 8, 1, [0,0,0,0,0,0,0,0]);
n1 := 1/2*Matrix(C, 8, 1, [0,0,0,1,1,1,1,0]);
n2 := 1/2*Matrix(C, 8, 1, [0,0,1,1,0,0,0,1]);
n3 := 1/2*Matrix(C, 8, 1, [0,0,1,0,1,0,1,1]);
n4 := n1+n2;
n5 := n1+n3;
n6 := n2+n3;
n7 := n1+n2+n3;
SchottkyN := [n0,n1,n2,n3,n4,n5,n6,n7];
M1 := [m1 + n: n in SchottkyN];
M2 := [m2 + n: n in SchottkyN];
M3 := [m3 + n: n in SchottkyN];
pi1 := 1;
pi2 := 1;
pi3 := 1;

function CharacteristicMatrixToPair(c)
  ZZ := Integers();
  QQ := Rationals();
  c *:= 2;
  c := [QQ!(ZZ!(GF(2)!(ZZ!el))) : el in Eltseq(c)];
  return [c[1..4], c[5..8]];
end function;

M1 := [CharacteristicMatrixToPair(el) : el in M1];
M2 := [CharacteristicMatrixToPair(el) : el in M2];
M3 := [CharacteristicMatrixToPair(el) : el in M3];

z := Eltseq(z);
for m in M1 do
    //pi1 := pi1* Theta(m, z, tau);
  pi1 := pi1*Theta(z, tau : char := m, prec := prec);
end for;

for m in M2 do
    //pi2 := pi2* Theta(m, z, tau);
  pi2 := pi2*Theta(z, tau : char := m, prec := prec);
end for;

for m in M3 do
    //pi3 := pi3* Theta(m, z, tau);
  pi3 := pi3*Theta(z, tau : char := m, prec := prec);
end for;

Schottky := pi1^2 + pi2^2 + pi3^2 - 2*(pi1*pi2 + pi2*pi3 + pi1*pi3);
Schottky;

// next, take derivative of Schottky modular form wrt to tau
// this will give the quadric


/* ------------------------------------------ */

CC4<X,Y,Z,W> := PolynomialRing(CC,4);
ells := [&+[c[1][i]*CC4.i : i in [1..g] ]: c in tritangents];
ell_ps := [&+[c[2][i]*CC4.i : i in [1..g] ]: c in tritangents];

// PP^3 attempt
Ccan, phi := CanonicalImage(S);
for cs_new in tritangents do
  print cs_new;
  TritangentSanityCheck(Ccan, cs_new);
end for;

/*
A := Matrix(3,3,mods);
A := Matrix(3,3,[1/el : el in mods]);
A^-1;
$1*Matrix(3,1,[-1,-1,-1]);
A^-1;
Ainv := $1;
Ainv*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);
lambdas := $1;
mods;
mods_mat := [[mods[i], mods[i+1], mods[i+2]] : i in [1,4,7]];
mods_mat;
mods_mat_scaled := [];
for i := 1 to 3 do
DiagonalMatrix(lambdas);
DiagonalMatrix(3,3,lambdas);
DiagonalMatrix;
lambdas;
DiagonalMatrix(Eltseq(lambdas));
L := $1;
B := Matrix(3,3,mods)*L;
Binv := Inverse(B);
Binv*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
ks := $1;
bitangents := [];
*/

/*
bitangents[1..7];
M;
M := Matrix(bitangents[1..7]);
NumericalKernel(Transpose(M));
NumericalKernel(M);
K := $1;
K;
Rows(K);
Krows := $1;
Krows[1];
[el/Krows[1][3] : el in Krows[1]];
[el/Krows[1][3] : el in Eltseq(Krows[1])];
[el/Krows[1][#Krows[1]] : el in Eltseq(Krows[1])];
[el/Krows[1][#Eltseq(Krows[1])] : el in Eltseq(Krows[1])];
[el/Krows[1][1] : el in Eltseq(Krows[1])];
bitangents[1..7];
CC3;
CC3.1;
fs := [[el[i]*CC3.i : i in [1..3] : el in bitangents]];
fs := [[el[i]*CC3.i : i in [1..3]] : el in bitangents];
fs[1];
fs := [&+[el[i]*CC3.i : i in [1..3]] : el in bitangents];
fs[1];
fs[7]^2;
fs_mat := [];
for f in fs do
fsq_math := [];
for f in fs do
fs[1];
Coefficients($1);
CoefficientsAndMonomials(fs[1]);
Coefficients;
AttachSpec("~/github/CHIMP/CHIMP.spec");
Abseltseq;
AbsEltseq;
AbsEltseq(fs[1]);
AbsEltseq:Maximal;
for f in fs do
cs := Coefficients(f);
while #cs lt 4 do
for f in fs do
fsq_math;
fsq_mat := [];
delete fsq_math;
for f in fs do
cs := Coefficients(f^2);
fs[4];
$^1;
fs[4]^2;
Coefficients($1);
#$1;
for f in fs do
cs := Coefficients(f^2);
while #cs lt 6 do
Append(~cs,0);
Coefficients(fs[2]^2);
fs[2]^2;
AbsEltseq;
MonomialsOfDegree(CC3,2);
mons := Eltseq($1);
MonomialsOfDegree(CC3,2);
mons := $1;
mons[3];
for m in mons do
print Coefficient(fs[2]^2,m);
end for;
Coefficient(fs[2]^2,CC3.2);
MonomialCoefficient;
for m in mons do
MonomialCoefficient(fs[2]^2,CC3.2);
end for;
for m in mons do
MonomialCoefficient(fs[2]^2,m);
end for;
for f in fs do
cs := [];
for m in mons do
Append(~cs, MonomialCoefficients(f^2,m));
end for;
end for;
for f in fs do
cs := [];
for m in mons do
Append(~cs, MonomialCoefficient(f^2,m));
end for;
end for;
#fsq_mat[1];
fsq_mat;
fsq_mat := [];
for f in fs do  
  cs := [];
  for m in mons do
    Append(~cs, MonomialCoefficient(f^2,m)); 
  end for;
  Append(~fsq_mat, cs);
end for;
fsq_mat[1];
fsq_mat := Matrix(fsq_mat);
NumericalKernel(fsq_mat);
ss := $1;
#ss;
ss[1];
ss := Eltseq(ss);
#ss;
fsq_mat;
Nrows(fsq_mat);
#fs;
fs := f[1..7];
fs := fs[1..7];
fsq_mat := [];
for f in fs do  
  cs := [];
  for m in mons do
    Append(~cs, MonomialCoefficient(f^2,m)); 
  end for;
  Append(~fsq_mat, cs);
end for;
fsq_mat := [];
for f in fs do  
  cs := [];
  for m in mons do
    Append(~cs, MonomialCoefficient(f^2,m)); 
  end for;
  Append(~fsq_mat, cs);
end for;
fsq_mat := Matrix(fsq_mat);
NumericalKernel(fsq_mat);
K := $1;
Dimension(K);
Eltseq(K);
ss := $1;
tritangents[1];
*/
