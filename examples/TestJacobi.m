

procedure JacobiNullwert(tau, odd_thetas, ~theta_list, ~ret)
  ZZ:=Integers();
  g := Nrows(tau);
  CC<I> := BaseRing(Parent(tau));
  prec := Precision(CC);
  pi := Pi(CC);
  derivs := [];
  for c in odd_thetas do
      chara := [ZZ!v : v in Eltseq(c)];
      chara := [chara[1..g], chara[g+1..(2*g)]];
      for i := 1 to g do
	if not IsDefined(theta_list, <c, i>) then
		dz := [0 : j in [1..g]];
        	dz[i] := 1;
                theta_list[<c, i>] :=  Theta([CC!0 : j in [1..g]], tau : char := chara, dz := [dz], prec := prec);
	 end if;
         Append(~derivs, theta_list[<c, i>]);
      end for;
  end for;
  M := Matrix(g,g,derivs);
  ret:=pi^(-g)*Determinant(M);
end procedure;


procedure TestJacobiDeriv(tau, S, ~theta_list)
        ZZ:=Integers();
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
	A:=Submatrix(S, [1..4], [1..4]);
	B:=Submatrix(S, [1..4], [5..8]);
	C:=Submatrix(S, [5..8], [1..4]);
	D:=Submatrix(S, [5..8], [5..8]);
	vec:= Vector(Diagonal(B*Transpose(A)) cat Diagonal(D*Transpose(C)));
	S1:=[n+vec: n in Transpose(S*N)[1..6]];
	S2:=[n+vec: n in Transpose(S*N2)[1..6]];
	odd_thetas:=[n+vec: n in Transpose(S*M)[1..4]];
	T1 := CC!1;
    	T2 := CC!1;
        for c in S1 do
	      if not IsDefined(theta_list, <c,0>) then
                chara := [ZZ!v : v in Eltseq(c)];
	        chara := [chara[1..4], chara[5..8]];
	        theta_list[<c, 0>] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
	      end if;
	      T1 *:= theta_list[<c, 0>];
        end for;

 	for c in S2 do
           if not IsDefined(theta_list, <c, 0>) then
              chara := [ZZ!v : v in Eltseq(c)];
              chara := [chara[1..4], chara[5..8]];
              theta_list[<c, 0>] :=  Theta([CC | 0,0,0,0], tau : char := chara, prec := prec);
           end if;
           T2 *:= theta_list[<c, 0>];
        end for;
        LHS :=CC!0;
        JacobiNullwert(tau, odd_thetas, ~theta_list, ~LHS);
        signsposs := [[1,1],[1,-1],[-1,1],[-1,-1]];
        RHSposs := [signs[1]*T1+signs[2]*T2: signs in signsposs];
        val, i := Min([Abs(LHS-RHSposs[i]): i in [1..4] ]);
        print val;
        print signsposs[i];        
end procedure;

