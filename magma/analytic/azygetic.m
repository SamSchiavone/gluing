
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

function Add1at2(vec)
	id:=IdentityMatrix(GF(2), 4);
	zer:=ZeroMatrix(GF(2), 4,4);
	if vec[2] eq 1 then
		return IdentityMatrix(GF(2), 8);
	else
		rhs:=Matrix(GF(2), [[1]]);
		lhs:=Matrix(GF(2),2,1, Eltseq(vec)[7..8]);
		sol:=Solution(rhs, lhs);
                Bx:=Matrix(GF(2), [[0,0,0,0],[0,0, sol[1,1], sol[2,1]],[0,sol[1,1],0,0],[0,sol[2,1], 0,0]]);
                Sx:=BlockMatrix(2,2,[[id, zer  ], [Bx, id]]);
                return Sx;
	end if;
end function;

function map_orthogonal(vec)
	v1:=Vector(Eltseq(vec)[3..4]);
        v2:=Vector(Eltseq(vec)[7..8]);
	id2:=IdentityMatrix(GF(2),2);
	zer2:=ZeroMatrix(GF(2),2,2);
	V2:=VectorSpace(GF(2),2);
	if v1 eq V2!0 and v2 eq V2!0 then
		B:=Matrix(GF(2), [[1,0],[0,0]]);
                ret:=BlockMatrix([[id2, B],[zer2, id2]]);
		retA:=Submatrix(ret, [1..2], [1..2]);
                retB:=Submatrix(ret, [1..2], [3..4]);
                retC:=Submatrix(ret, [3..4], [1..2]);
                retD:=Submatrix(ret, [3..4], [3..4]);
                return BlockMatrix(4,4, [[id2, zer2, zer2, zer2],[zer2, retA, zer2, retB],[zer2, zer2, id2, zer2], [zer2, retC, zer2, retD]]);
	elif v2 eq Vector([GF(2)|0,0]) then
		S1:=IdentityMatrix(GF(2),4);
	elif v1 eq Vector([GF(2)|0,0]) then
		S1:=BlockMatrix([[zer2, id2],[id2, zer2]]);
		v1:=v2;
	else
		B:=Matrix(GF(2), [[0,1],[1,0]]);
	        S1:=BlockMatrix([[id2, B],[zer2, id2]]);
	end if;
	i:=Min([i: i in [1..2]| v1[i] ne 0]);
	A:=Matrix(GF(2), 2,2, [ v1, V2.(3-i)]);
	S2:=BlockMatrix([[A^(-1), zer2], [zer2, Transpose(A)]]);
	ret:=S1*S2;
	retA:=Submatrix(ret, [1..2], [1..2]);
        retB:=Submatrix(ret, [1..2], [3..4]);
        retC:=Submatrix(ret, [3..4], [1..2]);
        retD:=Submatrix(ret, [3..4], [3..4]);
	return BlockMatrix(4,4, [[id2, zer2, zer2, zer2],[zer2, retA, zer2, retB],[zer2, zer2, id2, zer2], [zer2, retC, zer2, retD]]);
end function;

function map_azygetic(azy1, azy2)
	id:=IdentityMatrix(GF(2), 4);
        zer:=ZeroMatrix(GF(2), 4,4);
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
	S1add1at3:=Add1at2(vec1);
        S2add1at3:=Add1at2(vec2);
	vec1:=vec1*S1add1at3;
        vec2:=vec2*S2add1at3;
	Ortho1:=map_orthogonal(vec1);
        Ortho2:=map_orthogonal(vec2);
	/*print "\n vec1, Ortho1, vec1*Ortho1=";
	print vec1, Ortho1, vec1*Ortho1;
        print "\n vec2, Ortho2, vec2*Ortho2=";
	print vec2, Ortho2, vec2*Ortho1;*/
	return Transpose(M1^(-1)*S1add1at3* Ortho1*Ortho2^(-1) *S2add1at3^(-1)   *M2);
end function;

intrinsic special_fundamental_system(azy::SeqEnum) -> SeqEnum, SeqEnum, AlgMatElt
{Given an azygetic system of 4 odd characteristics, computes the two sets of 6 even characteristics such that the union is azygetic. It returns a symplectic transformation into the standard system as well  }
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
	G:=(Matrix(GF(2),6,6, [1: i in [1..36]])+IdentityMatrix(GF(2), 6));
	N2:=N*G;
	S:=map_azygetic(Transpose(M)[1..4], azy);
	A:=Submatrix(S, [1..4], [1..4]);
	B:=Submatrix(S, [1..4], [5..8]);
	C:=Submatrix(S, [5..8], [1..4]);
	D:=Submatrix(S, [5..8], [5..8]);
	vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
	return [n+vec: n in Transpose(S*N)[1..6]], [n+vec: n in Transpose(S*N2)[1..6]], S;
end intrinsic;

	
/* For testing:

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



for i in [1..20] do
N1, N2:=special_fundamental_system(Transpose(N)[1..3] cat [Vector(sys[i])]);
Special1:=Transpose(Matrix(Transpose(N)[1..3] cat [Vector(sys[i])] cat N1));
Special2:=Transpose(Matrix(Transpose(N)[1..3] cat [Vector(sys[i])] cat N2));
is_azygetic(Special1);
is_azygetic(Special2);
print N1, N2;
S:=map_azygetic(Transpose(N)[1..4], Transpose(N)[1..3] cat [Vector(sys[i])]);
A:=Submatrix(S, [1..4], [1..4]);
B:=Submatrix(S, [1..4], [5..8]);
C:=Submatrix(S, [5..8], [1..4]);
D:=Submatrix(S, [5..8], [5..8]);
vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
print [n+vec: n in Transpose(S*N)[1..4]] eq Transpose(N)[1..3] cat [Vector(sys[i])];
print Transpose(S)*J*S eq J;
print "\n";
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
*/

