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

/* TODO finish
intrinsic ShortestVectors(M::AlgMatElt)
  {}
  ZZ := Integers();
  R := BaseRing(Parent(M));
  p := -(Ceiling(Log(Max(radius, M))/Log(2))+4);
  n := Nrows(M);
  d := ZeroMatrix(ZZ, n, n);
  round_scale!(d, M, p); // TODO wut?
  L := Zlattice(d);
  U := shortest_vectors(L);
  return [change_base_ring(R, u)*M for u in U];
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
  //rho := L2Norm(ShortestVectors(Transpose(T))[1]*Sqrt(pi)); // TODO fix
  rho := 1;

  N := #dz;
  R0 := (Sqrt(g + 2*N + Sqrt(g^2 + 8*N)) + rho) div 2;

  T_inv_norm := L2Norm(Inverse(T));

  // We compute the radius of the ellipsoid over which we take the sum needed to bound the error in the sum by eps (See Theorem 3.1 in Agostini, Chua) 
  function R_function(x, eps)
    Rc := Parent(x);
    return (2*pi)^N * (g/2) * (2/rho)^g * &+[Binomial(N, j) * (pi^(j/2))^-1 * T_inv_norm^j * Sqrt(g)^(N - j) * Gamma(Rc!(g + j)/2, (x - rho/2)^2) : j in [0..N]];
      init := Rc!0 - eps; // TODO: what is this?
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

end intrinsic;
