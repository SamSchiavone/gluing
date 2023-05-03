IO:=Open("output.txt", "r");
thetabath := ReadObject(IO);
IO:=Open("tritangents-save.txt", "r");
tritangents := ReadObject(IO);
Q:=Matrix(thetabath[1]) + Transpose(Matrix(thetabath[1]));
//for i in [1..4] do
//  Q[i,i] /:=2;
//end for;
A:= Matrix([tr[1] : tr in tritangents[1..5]]); 
h:=[];
for i in [1..4] do
  h[i]:= Minor(A, [1..i-1] cat [5] cat [i+1..4], [1..4]);
end for;
A:= Submatrix(A, 1, 1, 4, 4);
D:= DiagonalMatrix(h);
A:=Transpose(A)*D;
A:= Inverse(A);
Coeff := Transpose(A)*Q*A;

IO2:=Open("output1.txt", "r");
thetabath2 := ReadObject(IO2);
IO2:=Open("tritangents2", "r");
tritangents2 := ReadObject(IO2);
Q2:=Matrix(thetabath2[1]) + Transpose(Matrix(thetabath2[1]));
//for i in [1..4] do
//  Q2[i,i] /:=2;
//end for;
A2:= Matrix([tr[1] : tr in tritangents2[1..5]]);
h2:=[];
for i in [1..4] do
  h2[i]:= Minor(A2, [1..i-1] cat [5] cat [i+1..4], [1..4]);
end for;
A2:= Submatrix(A2, 1, 1, 4, 4);
D2:= DiagonalMatrix(h2); 
A2:=Transpose(A2)*D2;
A2:= Inverse(A2);
Coeff2 := Transpose(A2)*Q2*A2;

