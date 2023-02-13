AttachSpec("~/github/gluing/magma/spec");
SetDebugOnError(true);
prec := 500;
 Q := RationalField();
 P<x> := PolynomialRing(Rationals());
R<s> := NumberField(P![1, 3, -6, -1, 1]); //change the totally real quartic field here

 PP<x> := PolynomialRing(R);
 RF<im> := NumberField(PP![1,0,1]); //change the integer in "CM by sqrt(-integer)" here

 CMF<a> := AbsoluteField(RF);

/*generating a principal polarization xi when class number is 1*/
 O := MaximalOrder(CMF);
D:= Different(O);
_, xi:= IsPrincipal(D^-1); //finding generator of inverse difference
xi;

/*making xi principal polarization of weil type*/
auts := Automorphisms(CMF);
cc := auts[2];
cmcounter := 0;
cmtype := [* *];
unitplaces := [1,1,1,1];
if cc(xi) eq xi then
    xi := xi*im;
end if;
if cc(xi*im) eq xi*im then
    for i in [2,4,6,8] do
        if Im(Conjugate(xi, i)) gt 0 then
        cmcounter := cmcounter +1;
        cmtype := Append(cmtype, Floor(i/2));
        end if;
    end for;

    if cmcounter eq 2 then
        printf "Xi gives a Weil type principal polarization";
        elif cmcounter eq 0 or cmcounter eq 4 then
            b := UnitsWithSigns(R, RealPlaces(R), [1,1,-1,-1])[1];
            b := CMF ! b;
            xi := xi*b;
            xi;
        elif cmcounter eq 1 then
            unitplaces[cmtype[1]] := -1;
            b_0 := UnitsWithSigns(R, RealPlaces(R), [1,1,-1,-1])[1];
            b_1 := UnitsWithSigns(R, RealPlaces(R), unitplaces)[1];
            b := CMF ! (b_0 * b_1);
            xi := xi*b;
            xi;
        else
            unitplaces[cmtype[1]] := -1;
            b := UnitsWithSigns(R, RealPlaces(R), unitplaces)[1];
            b := CMF ! b;
            xi := xi*b;
            xi;
        end if;
   else printf "xi is not totally imaginary";
end if;

/*Generating period matrix*/

Z := IntegerRing();
E := Matrix(Z, [[Trace(xi*cc(a)*b) : b in Basis(O)] : a in Basis(O)]);
D, C := FrobeniusFormAlternating(E);

g := Degree(CMF) div 2;

newb := ElementToSequence(Matrix(O,C)*Matrix(O,2*g,1,Basis(O)));

conj := [i  : i in [1..2*g] | Abs(Re(Conjugate(xi,i))) lt 10^-10 and Im(Conjugate(xi,i)) gt 0];
BigPM := Matrix([[Conjugate(b,i : Precision:=prec) : b in newb] : i in conj]);
tau := Submatrix(BigPM, 1, g + 1, g, g)^-1*Submatrix(BigPM, 1, 1, g, g);

function IntegerReduceMatrixG4(tau);
tauZZ := Matrix([ [ Round(Re(c)) : c in Eltseq(row) ] : row in Rows(tau) ]);
tauZZ[2,1] := tauZZ[1,2];
tauZZ[3,1] := tauZZ[1,3];
tauZZ[3,2] := tauZZ[2,3];
tauZZ[4,1] := tauZZ[1,4];
tauZZ[4,2] := tauZZ[2,4];
tauZZ[4,3] := tauZZ[3,4];
return tau - ChangeRing(tauZZ, BaseRing(tau));
end function;

function LLLReduceMatrixG4(tau);
//assert IsSmallPeriodMatrix(tau);
tau[2,1] := tau[1,2];
tau[3,1] := tau[1,3];
tau[3,2] := tau[2,3];
tau[4,1] := tau[1,4];
tau[4,2] := tau[2,4];
tau[4,3] := tau[3,4];
Imtau := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
_, T := LLLGram(Imtau);
return T*tau*Transpose(T);
end function;

function LeftActionHg(M, tau)
//assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
A := Submatrix(M, 1,  1, g,g); B := Submatrix(M, 1,  g+1, g,g);
C := Submatrix(M, g+1,1, g,g); D := Submatrix(M, g+1,g+1, g,g);
return (A*tau + B)*((C*tau + D)^(-1));
end function;

function ReduceSmallPeriodMatrixG4(tau)
N0 := Matrix(Integers(), [
[  0,  0,  0,  0,  1,  0,  0,  0],
[  0,  1,  0,  0,  0,  0,  0,  0],
[  0,  0,  1,  0,  0,  0,  0,  0],
[  0,  0,  0,  1,  0,  0,  0,  0],
[ -1,  0,  0,  0,  0,  0,  0,  0],
[  0,  0,  0,  0,  0,  1,  0,  0],
[  0,  0,  0,  0,  0,  0,  1,  0],
[  0,  0,  0,  0,  0,  0,  0,  1]
]);
counter := 0;
while true do
    counter +:= 1;
    /*
    if counter mod 10^3 eq 0 then
        vprint CurveRec : "";
        vprint CurveRec : "Counter for period matrix reduction is high:", counter;
        break;
    end if;
    */

    tau := LLLReduceMatrixG4(tau);



    tau := IntegerReduceMatrixG4(tau);


    if Abs(tau[1,1]) gt 1 then
        return tau;
    end if;
    tau := LeftActionHg(N0, tau);
end while;
end function;

tau := ReduceSmallPeriodMatrixG4(tau);
tau;



/*using magma to compute the Schottky form*/
C := ComplexField(prec);
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
  pi2 := pi1*Theta(z, tau : char := m, prec := prec);
end for;

for m in M3 do
    //pi3 := pi3* Theta(m, z, tau);
  pi3 := pi1*Theta(z, tau : char := m, prec := prec);
end for;

Schottky := pi1^2 + pi2^2 + pi3^2 - 2* (pi1*pi2 + pi2*pi3 + pi1*pi3);
Schottky;
