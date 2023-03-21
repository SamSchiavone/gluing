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

declare verbose Theta, 1;

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
  return Sqrt(&+[Abs(el)^2 : el in S]);
end intrinsic;

intrinsic L2Norm(v::ModTupFldElt) -> RngElt
  {}
  return L2Norm(Eltseq(v));
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

  ZZ := Integers();
  QQ := Rationals();
  g := Nrows(tau);

  if prec gt 0 then
    prec := Min([prec, Precision(Parent(z[1])), Precision(Parent(tau[1,1]))]);
  else
    prec := Min([Precision(Parent(z[1])), Precision(Parent(tau[1,1]))]);
  end if;
  SetDefaultRealFieldPrecision(prec);
  CC<I> := ComplexField(prec);
  RR := RealField(prec);
  pi := Pi(RR);

  if #z ne g then
    error "z must have length g";
  end if;

  // coercion stuff for theta chars
  if char eq [] then
    //char := [Matrix(GF(2),1,g,[0 : i in [1..g]]) : j in [1,2]];
    //char := [Matrix(QQ,1,g,[0 : i in [1..g]]) : j in [1,2]];
    char := [Vector([QQ!0 : i in [1..g]]) : j in [1,2]];
  end if;
  //if (#char[1] ne g) or (#char[2] ne g) then
  if Type(char[1]) eq SeqEnum then
    char := [Vector(QQ,el) : el in char];
  end if;
  if BaseRing(Parent(char[1])) eq Integers() then
    char := [ChangeRing(c,QQ) : c in char];
  end if;
  if (Ncols(char[1]) ne g) or (Ncols(char[2]) ne g) then
    error "characteristic must have length g";
  end if;

  // coercion stuff for derivs
  dz := [Vector(QQ,el) : el in dz];
  dtau := [Vector(QQ,el) : el in dtau];

  // tau := X + I*Y
  X := Real(tau);
  Y := Imaginary(tau);
  if not IsSymmetric(Y) then
    Y := (1/2)*(Y + Transpose(Y));
  end if;
  // Find T upper-triangular with transpose(T)*T = Y
  T := Transpose(Cholesky(Y));
  vprintf Theta: "Cholesky decomposition T = %o\n", T;
  printf "Y = %o\n", Y;
  printf "T^t * T = %o\n", Transpose(T)*T;
  printf "T * T^t = %o\n", T*Transpose(T);

  //L := LatticeWithBasis(Transpose(T));
  n := Floor(prec*Log(2)/Log(10));
  eps := RR!(10^-n);

  // In Agostini and Chua's code rho = is taken to be the square of the norm of the shortest vector times sqrt(pi)) for some reason. This could affect the error bounds
  vprint Theta: "Setting radius...";
  rho := L2Norm(ShortestVector(Lattice(Transpose(T)))*Sqrt(pi));
  //rho := L2Norm(ShortestVector(L)*Sqrt(pi));
  vprintf Theta: "rho = %o\n", rho;

  /*
  if #dz eq 0 then
    N := 0;
  else
    N := &+dz;
  end if;
  */
  N := #dz;
  R0 := (1/2)*(Sqrt(CC!(g + 2*N + Sqrt(CC!(g^2 + 8*N)))) + rho);
  printf "initial R0 = %o with precision %o\n", R0, Precision(R0);

  T_inv_norm := L2Norm(Inverse(T));

  // We compute the radius of the ellipsoid over which we take the sum needed to bound the error in the sum by eps (See Theorem 3.1 in Agostini, Chua) 
  function R_function(x, eps)
    Rc := Parent(x);
    //return -eps + (2*pi)^N * (g/2) * (2/rho)^g * &+[Binomial(N, j) * (pi^(j/2))^-1 * T_inv_norm^j * Sqrt(g)^(N - j) * (1-Gamma(Rc!(g + j)/2, (x - rho/2)^2) : Complementary := true) : j in [0..N]]; // Gamma or 1 - Gamma???
    return -eps + (2*pi)^N * (g/2) * (2/rho)^g * &+[Binomial(N, j) * (pi^(j/2))^-1 * T_inv_norm^j * Sqrt(g)^(N - j) * Gamma(Rc!(g + j)/2, (x - rho/2)^2 : Complementary := true) : j in [0..N]]; // Gamma or 1 - Gamma???
    //init := Rc!0 - eps;
  end function;

  /*
    We want to find the max(R0, x) where x is the solution to R_function(x, eps) = 0
    Taking x bigger will only improve the error bound, but we want x to be as small as possible
    to speed up later computations. We therefore look for an x for which R_function(x, eps) is small and negative.
    As R_function is monotonic descending we first determine an R1 for which R_function(R1, eps) < 0
    and then subdivide intervals until we find a solution that satisfies our requirements
  */

  vprint Theta: "Computing radius of ellipsoid";
  R1 := R0;
  err := RR!(10^(-n));
  vprintf Theta: "initializing R1 = %o\n", R1;

  // Find an R1 such that R_function becomes negative
  while R_function(R1, eps) gt 0 do
    R1 := R1 + R0;
    //R1 := 2*R1;
    vprintf Theta: "R1 = %o\n", R1;
    vprintf Theta: "R_function = %o\n", R_function(R1, eps);
  end while;
  vprintf Theta: "Making R_function negative, now R1 = %o\n", R1;

  if R1 ne R0 then
  //if not Abs(R1 - R0) lt eps then
    while (0 lt R_function(R0, eps)) or (R_function(R0, eps) lt -err) do
      Rmid := (R0 + R1)/2;
      middle := R_function(Rmid, eps);
      if middle lt -eps then
        R1 := Rmid;
      else
        R0 := Rmid;
      end if;
      //printf "Rmid = %o, R0 = %o, R1 = %o, middle = %o\n", Rmid, R0, R1, middle;
      //printf "R_function(R0, eps) = %o\n", R_function(R0, eps);
      //printf "eps = %o\n", eps;
    end while;
  end if;
  vprintf Theta: "After while loop, R0 = %o, R1 = %o\n", R0, R1;

  radius_ellipsoid := R1;
  printf "radius ellipsoid = %o\n", R1;
  error_epsilon := R_function(R1, RR!0);

  vprint Theta: "Computing lattice points in translates of ellipsoid";
  //L := Lattice(Y);
  //L := LatticeWithBasis(T);
  L := LatticeWithBasis(Transpose(T));
  //printf "basis = %o\n", Basis(L);
  //printf "R1^2/pi = %o\n", R1^2/pi;
  //printf "short vectors = %o\n", [ComplexField(10)!el[2] : el in ShortVectors(Lattice(T), R1^2/pi)];
  //ellipsoid_points := [Coordinates(L, L!Vector(el[1])) : el in ShortVectors(L, R1^2/pi)];
  ellipsoid_points := [Coordinates(L, L!Vector(el[1])) : el in ShortVectors(L, R1^2/pi)];
  ellipsoid_points cat:= [Eltseq(-Vector(v)) : v in ellipsoid_points];
  printf "initial #ellipsoid points = %o\n", #ellipsoid_points;
  for i := 1 to g do
    Lat1 := [];
    pad := Vector([0 : i in [1..g]]);
    pad[i] := 1;
    for pt in ellipsoid_points do
      Append(~Lat1, Eltseq(Vector(pt) + pad));
      Append(~Lat1, Eltseq(Vector(pt) - pad));
    end for;
    //printf "i = %o\n", i;
    //printf "new points = %o\n", [Vector(el) : el in Lat1];
    ellipsoid_points cat:= [Eltseq(Vector(el)) : el in Lat1];
  end for;
  ellipsoid_points := Setseq(Seqset(ellipsoid_points));
  printf "final #ellipsoid points = %o\n", #ellipsoid_points;
  //printf "ellipsoid points = %o\n", ellipsoid_points;
  if #Eltseq(dz) ne 0 then
    dz := Matrix(dz);
  end if;

  factor := CC!1;
  for ij in dtau do
    if ij[1] eq ij[2] then
      factor /:= 4*pi*I;
    else
      factor /:= 2*pi*I;
    end if;
    deriv := [[0 : i in [1..g]] : j in [1,2]];
    deriv[1][(ZZ!ij[1])] := 1;
    deriv[2][(ZZ!ij[2])] := 1;
    if #Eltseq(dz) ne 0 then
      dz := VerticalJoin([dz] cat [Matrix([el]) : el in deriv]);
    else
      dz := VerticalJoin([Matrix([el]) : el in deriv]);
    end if;
  end for;

  print "derivatives matrix after combining";
  print dz;

  if #Eltseq(dz) ne 0 then
    N := Nrows(dz);
  end if;
  
  // We seem to find more points than Agostini as we also consider lattices centered at points of the form [0,1,-1], etc. This could also affect error bounds

  // We compute the Theta function
  vprint Theta: "Computing theta function";
  x := Real(z);
  x := Vector(x);
  y := Imaginary(z);
  y := Matrix(#y,1,y);

  vprint Theta: "\tComputing exponential part";
  invYy := Inverse(Y)*y; // is y a row or column vector?
  exponential_part := Exp(pi*(Transpose(y)*invYy)[1,1]);
  vprintf Theta: "\t\t= %o\n", exponential_part;

  eta := Vector([Round(el) : el in Eltseq(invYy)]) - char[1]/2;
  pointset := [Vector(el) - eta : el in ellipsoid_points];
  pointset := [Matrix(g,1,Eltseq(el)) : el in pointset];

  //oscillatory_part = (2*piR*i)^N*sum([ prod(transpose(d)*v for d in dz; init = one(Rc)) * exp(piR*i*((transpose(v) * (X * v)) + 2*transpose(v) * (x + char[2]//2))) * exp(-piR* (transpose(v + invYy) * (Y * (v + invYy)))) for v in pointset]; init = zero(Cc))
  vprint Theta: "\tComputing oscillatory part";
  //oscillatory_part := (2*pi*I)^N*&+[CC | &*[CC | Transpose(d)*v : d in dz] * Exp(pi*I*((Transpose(v) * (X * v)) + 2*Transpose(v) * Matrix(g,1,(x + char[2]/2)))[1,1]) * Exp(-pi * (Transpose(v + invYy) * (Y * (v + invYy)))[1,1]) : v in pointset];
  oscillatory_part := RR!0;
  for v in pointset do
    vRR := ChangeRing(v, RR);
    //oscillatory_part +:= &*[CC | Transpose(d)*vRR : d in dz] * Exp(pi*I*((Transpose(vRR) * (X * vRR)) + 2*Transpose(vRR) * Matrix(g,1,(x + char[2]/2)))[1,1]) * Exp(-pi * (Transpose(vRR + invYy) * (Y * (vRR + invYy)))[1,1]);
    if #Eltseq(dz) ne 0 then
      oscillatory_part +:= &*[CC | (d*vRR)[1] : d in Rows(dz)] * Exp(pi*I*((Transpose(vRR) * (X * vRR)) + 2*Transpose(vRR) * Matrix(g,1,(x + char[2]/2)))[1,1]) * Exp(-pi * (Transpose(vRR + invYy) * (Y * (vRR + invYy)))[1,1]);
    else
      oscillatory_part +:= Exp(pi*I*((Transpose(vRR) * (X * vRR)) + 2*Transpose(vRR) * Matrix(g,1,(x + char[2]/2)))[1,1]) * Exp(-pi * (Transpose(vRR + invYy) * (Y * (vRR + invYy)))[1,1]);
    end if;
  end for;
  oscillatory_part *:= (2*pi*I)^N;
  vprintf Theta: "\t\t= %o\n", oscillatory_part;
  print "\n";

  result := factor*exponential_part*oscillatory_part;
  error_term := exponential_part*error_epsilon;

  return result, error_term;
end intrinsic;

intrinsic SiegelReduction(tau::AlgMatElt) -> Any
  {}

  g := Nrows(tau);
  CC<I> := BaseRing(tau);
  RR := BaseRing(Real(tau));
  QQ := Rationals();
  ZZ := Integers();

  vprint Theta: "Setting up block matrices";
  Aq := VerticalJoin(HorizontalJoin(ZeroMatrix(RR,1,1), ZeroMatrix(RR,1,g-1)), HorizontalJoin(ZeroMatrix(RR,g-1,1), IdentityMatrix(RR,g-1)));
  Bq := VerticalJoin(HorizontalJoin(-IdentityMatrix(RR,1), ZeroMatrix(RR,1,g-1)), HorizontalJoin(ZeroMatrix(RR,g-1,1), ZeroMatrix(RR,g-1,g-1)));
  Cq := -Bq;
  Dq := Aq;

  quasi_inversion := VerticalJoin(HorizontalJoin(Aq, Bq), HorizontalJoin(Cq,Dq));

  Gamma := IdentityMatrix(RR, 2*g);
  e := RR!0;

  vprint Theta: "Entering while loop";
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
      //i +:= 1;
    //end while;
    end for;
    vprintf Theta: "short = %o\n", short;

    if short ne 1 then
      S := SwapColumns(IdentityMatrix(RR,g),1,short);
      T := S*T;
      U := S*U;
    end if;

    Tt := Transpose(T);
    Y := T*Tt;

    Gamma := VerticalJoin(HorizontalJoin(U,ZeroMatrix(ZZ,g)), HorizontalJoin(ZeroMatrix(ZZ,g), (Transpose(U)^-1)))*Gamma;
    tau := U*Real(tau)*Transpose(U) + I*Y;
    X := Real(tau);

    B := Parent(X)!0;
    for i := 1 to Nrows(B) do
      for j := 1 to Ncols(B) do
        B[i,j] := Round(X[i,j]);
      end for;
    end for;
    tau -:= ChangeRing(B,CC);
    Gamma := VerticalJoin(HorizontalJoin(IdentityMatrix(RR,g), -B), HorizontalJoin(ZeroMatrix(RR,g,g), IdentityMatrix(RR,g)))*Gamma;
    e := Abs(tau[1,1]);
    vprintf Theta: "Now e = %o\n", e;
    if e gt 1 then
      return tau, Gamma;
    else
      Gamma := quasi_inversion*Gamma;
      tau := (Aq*tau + Bq)*((Cq*tau + Bq)^-1);
    end if;
  end while;
end intrinsic;
