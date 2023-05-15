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
sys1:=[s:s in sys | s[1,1] eq 0];
w:=Basis(K1)[1]*Transpose(M);
T:=Matrix(GF(2), [[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
 Blo:=BlockMatrix([[(Transpose(T))^(-1), zer],[zer, T]]);

T2:=Matrix(GF(2), [[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]]);
Blo2:=BlockMatrix([[id, zer],[T2, id]]);
Ortho:=Blo2*Blo;
sysnew:=[Transpose(Ortho)*sys[i]: i in [1..20]];
sys1new:=[Eltseq(s): s in sysnew| s[5,1] eq 0];
proj:=[[s[2..4],s[6..8]]: s in sys1new];



function is_azygetic(chars)
	Nchars:=Ncols(chars);
	i:=1;
	for j in [i+1..Nchars] do
		for k in [j+1..Nchars] do
			submat:=Submatrix(chars, [1..8], [i,j,k]);
			mat:= Transpose(submat)*J*submat;
			if mat[1,2]+mat[1,3]+mat[2,3] eq 0 then
				return false;
			end if;
		end for;
	end for;
	return true;
end function;

function map_azygetic(azy1, azy2)
	Rhs:=Vector([GF(2)|0, 0,1]);
	M1:=Matrix([azy1[1]+azy1[2], azy1[1]+azy1[3], azy1[1]+azy1[2]+azy1[3]+azy1[4]]);
	M1:=VerticalJoin(M1, Solution(J*Transpose(M1), Rhs));
	K1:=Matrix(Basis(Kernel(J*Transpose(M1))));
	V1:=VectorSpace(GF(2), 4, K1*J*Transpose(K1));
	HypDec1:=HyperbolicSplitting(V1);
	hypmat1:=Matrix([HypDec1[1][1][1], HypDec1[1][1][2], HypDec1[1][2][1], HypDec1[1][2][2]]);
	M1:=VerticalJoin(M1, hypmat1*K1);
	M1:=Submatrix(M1, [1, 3,5,7,2,4,6,8], [1..8]);

	M2:=Matrix([azy2[1]+azy2[2], azy2[1]+azy2[3], azy2[1]+azy2[2]+azy2[3]+azy2[4]]);
	M2:=VerticalJoin(M2, Solution(J*Transpose(M2), Rhs));
	K2:=Matrix(Basis(Kernel(J*Transpose(M2))));
	V2:=VectorSpace(GF(2), 4, K2*J*Transpose(K2));
	HypDec2:=HyperbolicSplitting(V2);
	hypmat2:=Matrix([HypDec2[1][1][1], HypDec2[1][1][2], HypDec2[1][2][1], HypDec2[1][2][2]]);
	M2:=VerticalJoin(M2, hypmat2*K2);
	M2:=Submatrix(M2, [1, 3,5,7,2,4,6,8], [1..8]);
	S:=M1^(-1)*M2;
	if &or[(Matrix(azy1)*S-Matrix(azy2))[1] ne (Matrix(azy1)*S-Matrix(azy2))[i]: i in [2..4]] then
		error("Error: First symplectic Transformation wrong");
	end if;
	vec1:=(Matrix(azy1)*M1^(-1))[1];
	vec2:=(Matrix(azy2)*M2^(-1))[1];
	A1:=Submatrix(Transpose(M1^(-1)), [1..4], [1..4]);
	B1:=Submatrix(Transpose(M1^(-1)), [1..4], [5..8]);
	C1:=Submatrix(Transpose(M1^(-1)), [5..8], [1..4]);
	D1:=Submatrix(Transpose(M1^(-1)), [5..8], [5..8]);
        A2:=Submatrix(Transpose(M2^(-1)), [1..4], [1..4]);
	B2:=Submatrix(Transpose(M2^(-1)), [1..4], [5..8]);
	C2:=Submatrix(Transpose(M2^(-1)), [5..8], [1..4]);
	D2:=Submatrix(Transpose(M2^(-1)), [5..8], [5..8]);
	vec1+:= Vector(Diagonal(B1*Transpose(A1)) cat Diagonal(D1*Transpose(C1)));
	vec2+:= Vector(Diagonal(B2*Transpose(A2)) cat Diagonal(D2*Transpose(C2)));
	trans:= vec1+vec2;
	//TODO still does not work for all cases of the pairing
	if &+[Eltseq(vec1)[i]*Eltseq(trans)[i+4]: i in [1..4]] eq 1 then
		v1:=Matrix(1,2, Eltseq(vec1)[3..4]  );
                v2:=Matrix(1,2, Eltseq(vec1)[7..8]  );
                w2:=Matrix(1,2, Eltseq(trans)[7..8]  );
		if v2+w2 ne 0 then
			min1:=Min([i: i in [1..2]| Eltseq(v1)[i] ne 0]);
	                min2:=Min([i: i in [1..2]| Eltseq(v2+w2)[i] ne 0]);
			Vd2:=VectorSpace(GF(2),2);
			Mat1:=Matrix(2,2, [Vd2!v1, Vd2.min1]);
	                Mat2:=Matrix(2,2, [Vd2!(v2+w2), Vd2.min1]);
		        if (v1*Transpose(v2))[1,1] eq 0 then
				Arot:=Mat2^(-1)*Matrix(GF(2), [[0,1],[1,0]])*Transpose(Mat1^(-1));
			else 
				Arot:=Mat2^(-1)*Transpose(Mat1^(-1));
			end if;
		        id2:=IdentityMatrix(GF(2),2);
			Srot:=DiagonalJoin([id2, Arot, id2, Transpose(Arot)^(-1)]);
		else 
			v1:=Matrix(1,3, Eltseq(vec1)[2..4]  );
	                v2:=Matrix(1,3, Eltseq(vec1)[6..8]  );
			Rhs:=Vector([(v1*Transpose(v2))[1,1]+v1[1,1]]);
			col:=Solution( Matrix(2,1, Eltseq(vec1)[3..4])  , Rhs);
			Srot:=IdentityMatrix(GF(2),8);
			Srot[3,2]:=col[1];
			Srot[4,2]:=col[2];
			Srot[6,7]:=col[1];
			Srot[6,8]:=col[2];
		end if;
		vec1 := vec1 * Srot;
		trans:= vec1+vec2;
	else 
		Srot:=IdentityMatrix(GF(2),8);
	end if;
	matx:=Matrix([[1+vec1[6], vec1[7], vec1[8],0,0,0],[0, vec1[6], 0, 1+vec1[7], vec1[8],0], [0, 0, vec1[6], 0, vec1[7], 1+vec1[8]]] );
	print "\n";
	print matx, vec1, vec2, trans;
	sol:=Solution(Transpose(matx), Vector(Eltseq(trans)[2..4]));
	Bx:=Matrix([[0,0,0,0],[0,sol[1], sol[2], sol[3]],[0,sol[2], sol[4], sol[5]],[0,sol[3], sol[5], sol[6]]]);
	Sx:=BlockMatrix(2,2,[[id, Matrix([[0,0,0,0],[0,sol[1], sol[2], sol[3]],[0,sol[2], sol[4], sol[5]],[0,sol[3], sol[5], sol[6]]])  ], [zer, id]]);

	vec1 +:= Vector(Eltseq(trans)[1..4] cat [GF(2)!0: i in [1..4]]  );
	maty:=Matrix([[1+vec1[3], vec1[4], 0],[0, vec1[3], 1+vec1[4]]] );
	sol:=Solution(Transpose(maty), Vector(Eltseq(trans)[7..8]));
	Sy:=BlockMatrix(2,2,[[id, zer], [Matrix([[0,0,0,0],[0,0,0,0],[0,0,sol[1], sol[2]], [0,0,sol[2], sol[3]]]), id]]);
	return Transpose(M1^(-1)*Srot*Transpose(Sx)*Transpose(Sy)*M2);

end function;

function special_fundamental_system(azy)
	zer:=ZeroMatrix(GF(2), 4,4);
	zer1:=ZeroMatrix(GF(2), 4,1);
	zer2:=ZeroMatrix(GF(2), 4,2);
	id:=IdentityMatrix(GF(2),4);
	J:=BlockMatrix(2,2, [zer,id, id,zer]);
	J1:=BlockMatrix(2,2, [zer,id, zer,zer]);
	triang:=zer;
	for i in [1..4] do
		for j in [i..4] do
			triang[i,j]:=1;
		end for;
	end for;
	M:=VerticalJoin(id, triang);

	N:=VerticalJoin(HorizontalJoin(id, zer2), HorizontalJoin(zer1, HorizontalJoin(triang, zer1)) );
	S:=map_azygetic(Transpose(M)[1..4], azy);
	A:=Submatrix(S, [1..4], [1..4]);
	B:=Submatrix(S, [1..4], [5..8]);
	C:=Submatrix(S, [5..8], [1..4]);
	D:=Submatrix(S, [5..8], [5..8]);
	vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
	return [n+vec: n in Transpose(S*N)[1..4]];
end function;

for i in [9] do
N1:=special_fundamental_system(Transpose(N)[1..3] cat [Vector(sys[i])]);
Special:=Transpose(Matrix(Transpose(N)[1..3] cat [Vector(sys[i])] cat N1));
is_azygetic(Special);
end for;


S:=map_azygetic(Transpose(N)[1..4], Transpose(N)[1..3] cat [Vector(sys[1])]);
A:=Submatrix(S, [1..4], [1..4]);
B:=Submatrix(S, [1..4], [5..8]);
C:=Submatrix(S, [5..8], [1..4]);
D:=Submatrix(S, [5..8], [5..8]);
vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
print [n+vec: n in Transpose(S*N)[1..4]];


S:=map_azygetic(Transpose(N)[1..4], Transpose(N)[1..3] cat [Vector(sys[3])]);
A:=Submatrix(S, [1..4], [1..4]);
B:=Submatrix(S, [1..4], [5..8]);
C:=Submatrix(S, [5..8], [1..4]);
D:=Submatrix(S, [5..8], [5..8]);
vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
print [n+vec: n in Transpose(S*N)[1..4]];


