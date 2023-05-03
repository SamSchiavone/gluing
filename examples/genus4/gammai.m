AttachSpec("/home/apieper/gluing-g4/magma/spec");
prec := 30;
SetDefaultRealFieldPrecision(prec);
QQ:=Rationals();
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
CC:=BaseRing(tau);
Quadric:= Matrix(CC, 4,4, [[1,0,0,4/7],[0,4/7,0,0], [0,0, 4/7, 0], [4/7, 0,0,2/7]]);

v := Vector([GF(2) | 0, 0, 0, 0, 1, 0, 0, 0 ]);
steiner := [];
steinerrie:=[];
//1st 7 form an aronhold set.
aronhold:=[[[1,1,1],[1,1,1]],[[0,0,1],[0,1,1]],[[0,1,1],[0,0,1]],[[1,0,1],[1,0,0]],[[1,0,0],[1,0,1]],[[1,1,0],[0,1,0]],[[0,1,0],[1,1,0]], [[0,1,0],[0,1,0]], [[1,0,0],[1,1,0]], [[1,1,0],[1,0,0]]];
rie:=[[[0,0,0],[0,0,1]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,1]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]], [[0,0,1],[1,1,0]],[[0,0,0],[1,1,0]],[[0,0,0],[0,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[1,0,0]]];
r:=10;



for aro in aronhold[1..r] do

    cc1 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!0] cat [GF(2)!el: el in  aro[2]];
    cc2 := [GF(2)!0] cat  [GF(2)!el : el in aro[1]] cat [GF(2)!1] cat [GF(2)!el: el in  aro[2]];
    c1:=[[0] cat aro[1] , [0] cat aro[2]];
    c2:=[[0] cat aro[1] , [1] cat aro[2]];

    Append(~steiner, [c1,c2]);


end for;



tritangents := [];
//for char in chars[1..4] do
//for char in chars do
for s in steiner do
  new := [];
  for char in s do
      print(char);
    Append(~new, TritangentPlane(Pi_big, char));
    //Append(~new, TritangentPlane(HorizontalJoin(IdentityMatrix(CC, 4), tau), char));
  end for;
  Append(~tritangents, new);
end for;



chars_even := EvenThetaCharacteristics(3);
/*eps := chars_even[1];
eta_sqs := AssociativeArray(); // now just etas
for delta in chars_even do
  print delta;
  delta_new1 := [[0] cat el : el in delta];
  delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
  // formula from Lemma 1 (p. 14) of Farkas
  eta_sqs[delta] := Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2);  // TODO: fix signs here
end for;
*/

thetas := [];
for i := 1 to 64 do
    s := Intseq(i mod 64,2,6);
    s := Reverse(s);
    delta := [s[1..3], s[4..6]];
    if delta in chars_even then
        print delta;
        delta_new1 := [[0] cat el : el in delta];
        delta_new2 := [[0] cat delta[1], [1] cat delta[2]];
      // formula from Lemma 1 (p. 148) of Farkas
        thetas[i] := Sqrt(Theta([CC!0 : i in [1..g]], tau : char := delta_new1)*Theta([CC!0 : i in [1..g]], tau : char := delta_new2));
    else
        Append(~thetas, 0);
    end if;
end for;
thetas:=correct_signs(thetas);



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

A := Matrix(3,3,[1/el : el in mods]);
Ainv:=A^-1;
Ainv*Matrix(3,1,[BaseRing(Parent(Ainv)) | -1,-1,-1]);
lambdas := $1;

mods_mat := [[mods[i], mods[i+1], mods[i+2]] : i in [1,4,7]];
L:=DiagonalMatrix(Eltseq(lambdas));
B := Matrix(3,3,mods)*L;
Binv := Inverse(B);
ks:=Binv*Matrix(3,1,[BaseRing(Parent(Binv)) | -1,-1,-1]);
bitangents := [];
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
for i := 1 to 3 do
  new := u0/(mods_mat[1,i]*(1-ks[i,1]*mods_mat[2,i]*mods_mat[3,i])) + u1/(mods_mat[2,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[3,i])) + u2/(mods_mat[3,i]*(1-ks[i,1]*mods_mat[1,i]*mods_mat[2,i]));
  Append(~bitangents, Coefficients(new));
end for;




fs := [&+[el[i]*CC3.i : i in [1..3]] : el in bitangents[1..r]];
mons:=MonomialsOfDegree(CC3,2);


fsq_mat := [];
for f in fs do
  cs := [];
  for m in mons do
    Append(~cs, MonomialCoefficient(f^2,m));
  end for;
  Append(~fsq_mat, cs);
end for;
fsq_mat := Matrix(fsq_mat);
si := NumericalKernel(fsq_mat);
sirows:=Nrows(si);
/*sum_i gamma_i s_1i l_i l_i'= Q
sum_i gamma_i s_ji l_i l_i'= a_j Q

Compute gammai with Q:*/
mats1:=[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2]) * si[1][i]): i in [1..r]] cat [ ZeroMatrix(CC, 4,4): j in [1..sirows-1]];
mats1:=[(m +Transpose(m))/2 : m in mats1];
W:= Matrix( [&cat[[CC!Quadric[i,j]: j in [i..4]]: i in [1..4]] cat [CC!0: j in [1..10*(sirows-1)]  ]] );
mats2:=[[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2]) * si[j][i]): i in [1..r]]  cat [ Quadric: nu in [1..sirows-1]] : j in [2..sirows]];
mats2:=[[(m +Transpose(m))/2 : m in mats2[j]] : j in [1..sirows-1]];
mats:=[mats1] cat mats2;
X:=HorizontalJoin([Matrix(  [&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats[j]]  ): j in [1..sirows]]);
gammai:=NumericalSolution(X, W: Epsilon:=RR!10^(-15));
gammai:=Eltseq(gammai1/gammai1[1,1]);


//


//Compute gammai without Q:
mats1new:=[(Matrix(4,1, tritangents[i][1])*Matrix(1,4, tritangents[i][2])): i in [1..r]];
mats1new:=[(m +Transpose(m))/2 : m in mats1new];
Xnew:=Matrix([&cat[[m[i,j]: j in [i..4]]: i in [1..4]] : m in mats1new]  );
vi:=NumericalKernel(Xnew: Epsilon:=RR!10^(-15));
N:=HorizontalJoin([  DiagonalMatrix(Eltseq(vi[i]))*fsq_mat : i in [1..Nrows(vi) ]]);
gammaiinv:=NumericalKernel(N: Epsilon:=RR!(10^(-15)));

