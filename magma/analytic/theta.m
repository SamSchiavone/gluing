// this should probably go in curve_reconstruction/magma/theta_derivs.m eventually

// auxiliary functions

intrinsic ComplexConjugate(A::AlgMatElt) -> AlgMatElt
  {}
  return Matrix(Nrows(A), Ncols(A), [[ComplexConjugate(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end intrinsic;

intrinsic Real(S::SeqEnum) -> RngElt
{}
  return [Real(el) : el in S];
end intrinsic;

intrinsic Imaginary(S::SeqEnum) -> RngElt
{}
  return [Imaginary(el) : el in S];
end intrinsic;

intrinsic Real(A::AlgMatElt) -> AlgMatElt
  {}
  return Matrix(Nrows(A), Ncols(A), [[Real(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end intrinsic;

intrinsic Imaginary(A::AlgMatElt) -> AlgMatElt
  {}
  return Matrix(Nrows(A), Ncols(A), [[Imaginary(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end intrinsic;

intrinsic L2Norm(S::SeqEnum[RngElt]) -> RngElt
  {}
  return Sqrt(&+[el^2 : el in S]);
end intrinsic;

intrinsic L2Norm(A::AlgMatElt) -> RngElt
  {}
  return L2Norm(&cat[[A[i,j] : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end intrinsic;

/*
  intrinsic ShortestVectors(M::AlgMatElt)
    {}
    L := Lattice(M);
    return ShortestVectors(L);
  end intrinsic;
*/

// theta series

intrinsic Theta(z::SeqEnum[FldComElt], tau::AlgMatElt : char := [], dz := [], dtau := [], prec := 0) -> SeqEnum
  {}

  g := Nrows(tau);

  if prec gt 0 then
    prec := Min([prec, Precision(Parent(z[1])), Precision(Parent(tau[1,1]))]);
  else
    prec := Min([Precision(Parent(z[1])), Precision(Parent(tau[1,1]))]);
  end if;
  CC<I> := ComplexField(prec);
  RR := RealField(prec);
  pi := Pi(RR);

  if #z ne g then
    error "z must have length g";
  end if;

  if char eq [] then
    char := [Matrix(GF(2),1,4,[0 : i in [1..g]]) : j in [1,2]];
  end if;
  if (#char[1] ne g) or (#char[2] ne g) then
    error "characteristic must have length g";
  end if;

  // tau := X + I*Y
  X := Real(tau);
  Y := Imaginary(tau);
  // Find T upper-triangular with transpose(T)*T = Y
  T := Transpose(Cholesky(Y));

 n := Floor(prec*Log(2)/Log(10));
 eps := RR(10^-n);

  // In Agostini and Chua's code rho = is taken to be the square of the norm of the shortest vector times sqrt(pi)) for some reason. This could affect the error bounds
  rho := L2Norm(ShortestVectors(Lattice(Transpose(T)))*Sqrt(pi));

  N := #dz;
  R0 := (1/2)*(Sqrt(g + 2*N + Sqrt(g^2 + 8*N)) + rho);

  T_inv_norm := L2Norm(Inverse(T));

  // We compute the radius of the ellipsoid over which we take the sum needed to bound the error in the sum by eps (See Theorem 3.1 in Agostini, Chua) 
  function R_function(x, eps)
    Rc := Parent(x);
    return (2*pi)^N * (g/2) * (2/rho)^g * &+[Binomial(N, j) * (pi^(j/2))^-1 * T_inv_norm^j * Sqrt(g)^(N - j) * Gamma(Rc!(g + j)/2, (x - rho/2)^2) : j in [0..N]];
      //init := Rc!0 - eps;
  end function;

  /*
    We want to find the max(R0, x) where x is the solution to R_function(x, eps) = 0
    Taking x bigger will only improve the error bound, but we want x to be as small as possible
    to speed up later computations. We therefore look for an x for which R_function(x, eps) is small and negative.
    As R_function is monotonic descending we first determine an R1 for which R_function(R1, eps) < 0
    and then subdivide intervals until we find a solution that satisfies our requirements
  */

  R1 := R0;
  err := RR!(10^(-n));

  // Find an R1 such that R_function becomes negative
  while R_function(R1, eps) gt 0 do
    R1 := R1 + R0;
  end while;

  if R1 ne R0 then
    while (0 lt R_function(R0, eps)) or (R_function(R0, eps) lt -err) do
      Rmid := (R0 + R1)/2;
      middle := R_function(Rmid, eps);
      if middle lt 0 then
        R1 := Rmid;
      else
        R0 := Rmid;
      end if;
    end while;
  end if;

  radius_ellipsoid := R1;
  error_epsilon := R_function(R1, RR!0);

  ellipsoid_points := [el[1] : el in ShortVectors(Lattice(Y), R1^2/pi)];
  for i := 1 to g do
    Lat1 := [];
    pad := Vector([RR!0 : i in [1..g]]);
    pad[i] := 1;
    for pt in ellipsoid_points do
      Append(~Lat1, pt + pad);
      Append(~Lat1, pt - pad);
    end for;
    ellipsoid_points cat:= Lat1;
  end for;

  factor := CC!1;
  for ij in dtau do
    if ij[1] eq ij[2] then
      factor /:= 4*pi*I;
    else
      factor /:= 2*pi*I;
    end if;
    deriv := [[0 : i in [1..g]] : j in [1,2]];
    deriv[1][ij[1]] := 1;
    deriv[2][ij[2]] := 1;
    dz := VerticalJoin(dz,deriv);
  end for;
  
  // We seem to find more points than Agostini as we also consider lattices centered at points of the form [0,1,-1], etc. This could also affect error bounds

  // We compute the Theta function
  x := Real(z);
  y := Imaginary(z);

  invYy := Inverse(Y)*y; // is y a row or column vector?
  exponential_part := Exp(pi*(Transpose(y)*invYy));

  //eta := ; ? TODO round all entries of invYy

  //pointset := ; ? TODO
  pointset := [];

  //oscillatory_part = (2*piR*i)^N*sum([ prod(transpose(d)*v for d in dz; init = one(Rc)) * exp(piR*i*((transpose(v) * (X * v)) + 2*transpose(v) * (x + char[2]//2))) * exp(-piR* (transpose(v + invYy) * (Y * (v + invYy)))) for v in pointset]; init = zero(Cc))
  oscillatory_part := (2*pi*I)^N*&+[&*[Transpose(d)*v : d in dz] * Exp(pi*I*((Transpose(v) * (X * v)) + 2*Transpose(v) * (x + char[2]/2))) * Exp(-pi * (Transpose(v + invYy) * (Y * (v + invYy)))) : v in pointset];

  result := factor*exponential_part*oscillatory_part;
  error_term := exponential_part*error_epsilon;

  return result, error_term;
end intrinsic;

intrinsic SiegelReduction(tau::AlgMatElt) -> Any
  {}

  g := Nrows(tau);
  CC<I> := BaseRing(tau);
  RR := BaseRing(Real(tau));

  Aq := VerticalJoin(HorizontalJoin(ZeroMatrix(RR,1,1), ZeroMatrix(RR,g-1,1)), HorizontalJoin(ZeroMatrix(RR,g-1,1), IdentityMatrix(RR,g-1)));
  Bq := VerticalJoin(HorizontalJoin(-IdentityMatrix(1), ZeroMatrix(RR,1,g-1)), HorizontalJoin(ZeroMatrix(RR,g-1,1), ZeroMatrix(RR,g-1,g-1)));
  Cq := -Bq;
  Dq := Aq;

  quasi_inversion := VerticalJoin(HorizontalJoin(Aq, Bq), HorizontalJoin(Cq,Dq));

  Gamma := IdentityMatrix(RR, 2*g);
  e := RR!0;

  while e le 1 do
    Y := Imaginary(tau);
    Y := (Y + Transpose(Y))/2; // make sure matrix is symmetric
    T := Cholesky(Imaginary(tau));
    T, U := LLL(T);
    Tt := Transpose(T);

    short := 1;
    //i := 1;
    //while i le g do
    for i := 1 to g do
      if L2Norm(Rows(T)[short]) gt L2Norm(Rows(T)[i]) then
        short := i;
      end if;
      i +:= 1;
    //end while;
    end for;

    if short ne 1 then
      S := SwapColumns(IdentityMatrix(RR,g),1,short);
      T := S*T;
      U := S*U;
    end if;

    Tt := Transpose(T);
    Y := T*Tt;

    Gamma := VerticalJoin(HorizontalJoin(U,ZeroMatrix(RR,g)), HorizontalJoin(ZeroMatrix(RR,g), Inverse(Transpose(U))))*Gamma;
    tau := U*Real(tau)*Transpose(U) + I*Y;
    X := Real(tau);

    B := Parent(X)!0;
    for i := 1 to Nrows(B) do
      for j := 1 to Ncols(B) do
        B[i,j] := Round(X[i,j]);
      end for;
    end for;
    tau -:= B;
    Gamma := VerticalJoin(HorizontalJoin(IdentityMatrix(RR,g), -B), HorizontalJoin(ZeroMatrix(RR,g,g), IdentityMatrix(RR,g)))*Gamma;
    e := Abs(tau[1,1]);
    if e gt 1 then
      return tau, Gamma;
    else
      Gamma := quasi_inversion*Gamma;
      tau := (Aq*tau + Bq)*Inverse(Cq*tau + Bq);
    end if;
  end while;
end intrinsic;
