AttachSpec("spec");
AttachSpec("~/github/CHIMP/CHIMP.spec");
Attach("~/github/gluing/magma/analytic/reconstruction.m");
P<x,y,z> := PolynomialRing(QQ,3);
g := 3;
X := KleinQuartic(Proj(P));
f := DefiningEquations(AffinePatch(X,1))[1];
T := RiemannSurface(f);
pi := SmallPeriodMatrix(T);
CC<I> := BaseRing(Parent(pi));
prec := Precision(CC);
z := [CC!0 : i in [1..g]];
SetProfile(true);
t0 := Cputime();
Theta(z, pi : prec := prec);
t1 := Cputime();
printf "Theta took %o\n", t1-t0;
G := ProfileGraph();
ProfilePrintByTotalTime(G : Percentage);

thetas2 := [];
for i := 1 to 64 do
  s := Intseq(i mod 64,2,6);
  s := Reverse(s);
  delta := [s[1..3], s[4..6]];
  //delta := (1/2)*Matrix(QQ,2*g,1,s);
  if &+[delta[1][i]*delta[2][i] : i in [1..3]] mod 2 eq 0 then
    //Append(~thetas2, Theta(delta, Matrix(g,1,z), pi));
    Append(~thetas2, Theta(z, pi : char := delta, prec := prec));
  else
    Append(~thetas2, 0);
  end if;
end for;

thetas3 := [];
for i := 1 to 64 do
  s := Intseq(i mod 64,2,6);
  s := Reverse(s);
  //delta := [s[1..3], s[4..6]];
  delta := (1/2)*Matrix(QQ,2*g,1,s);
  //if &+[delta[1][i]*delta[2][i] : i in [1..3]] mod 2 eq 0 then
    Append(~thetas3, Theta(delta, Matrix(g,1,z), pi));
    //Append(~thetas2, Theta(z, pi : char := delta, prec := prec));
  //else
    //Append(~thetas3, 0);
  //end if;
end for;



for c in EvenThetaCharacteristics(3) do
  Append(~thetas2, Theta(z, pi : char := c, prec := prec));
end for;
