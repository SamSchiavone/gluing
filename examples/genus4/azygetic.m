V:=VectorSpace(GF(2), 8);
A:=Matrix(GF(2), [[1: j in [1..2*i-1]] cat [0: j in [1..8-2*i+1]  ]: i in [1..4] ] cat [[1: j in [1..2*i-2]] cat [0, 1] cat [0: j in [1..8-2*i]  ]: i in [1..4] ]  );
zer:=ZeroMatrix(GF(2), 4,4);
id:=IdentityMatrix(GF(2),4);
J:=BlockMatrix(2,2, [zer,id, id,zer]);
J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
U:=Transpose(A)^(-1);
Transpose(U)*A*U;
v:=Matrix(GF(2), 8,1, [1: i in [1..8]]);
U:=HorizontalJoin(U, v);
new:=[];
for v in V do
	X:=Transpose(Matrix([v: i in [1..9 ]]   ));
	diag:=Diagonal(Transpose(X+U)*J1*(X+U)  );
	Append(~new, &+[Integers()!diag[i]: i in [1..9]]);
end for;
_, ind:=Max(new);
v:=[v: v in V][ind];
X:=Transpose(Matrix([v: i in [1..9 ]]   ));
N:=X+U;

N3:=Transpose(N)[1..3];
q4:=Transpose(N)[4];
q6:=Transpose(N)[6];
for i in [1..3] do
	N3[i] +:= q4;
end for;
M:=Transpose(Matrix(N3));
K:=Kernel(J*M);
K1:=Kernel(Transpose(M)*J*M);

cos:=[Matrix(8,1, Eltseq(q6+k)): k in K];
sys:=[x: x in cos| (Transpose(x)*J1*x)[1,1] eq 1];
eta:=M*Matrix(3, 1, Basis(K1)[1]);
sys1:=[s:s in sys s[1,1] eq 0];
w:=Basis(K1)[1]*Transpose(M);
T:=Matrix(GF(2), [[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
 Blo:=BlockMatrix([[(Transpose(T))^(-1), zer],[zer, T]]);

T2:=Matrix(GF(2), [[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]]);
Blo2:=BlockMatrix([[id, zer],[T2, id]]);
Ortho:=Blo2*Blo;
sysnew:=[Transpose(Ortho)*sys[i]: i in [1..20]];
sys1new:=[Eltseq(s): s in sysnew| s[5,1] eq 0];
proj:=[[s[2..4],s[6..8]]: s in sys1new];
