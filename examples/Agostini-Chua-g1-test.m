/*
https://github.com/chualynn/Theta.jl/blob/master/test/theta_test.jl
@testset "Genus 1" begin
    z = [0.85746943+0.30772689im];
    τ = [0.63443932+0.49919024im];
    Y = imag(τ);
    x = real(z);
    y = convert(Array{Float64}, imag(z));
    Yinv = inv.(Y);
    y0 = Yinv.*y;
    R = RiemannMatrix(τ, siegel=false);
    @test Theta.oscillatory_part(R, x, y0, [0]) ≈ -0.257+0.197im atol=ϵ
    @test theta(z, R) ≈ -0.466+0.358im atol=ϵ
    @test theta(z, R, derivs=[[1]]) ≈ 1.724+9.922im atol=ϵ
    @test theta(z, R, derivs=[[1],[1]]) ≈ 67.973-10.642im atol=ϵ
    @test theta(z, R, char=[[1],[0]]) ≈ -1.770-1.225im atol=ϵ
    @test theta(z, R, char=[[1],[1]]) ≈ -1.669+0.298im atol=ϵ
    @test theta(z, R, char=[[1],[0]], derivs=[[1]]) ≈ -2.368+6.931im atol=ϵ
    @test theta(z, R, char=[[1],[1]], derivs=[[1]]) ≈ -0.488+6.871im atol=ϵ
    @test symplectic_transform(siegel_transform(τ)[1], τ) ≈ siegel_transform(τ)[2]
end
*/

//QQ := Rationals();
//R<x,y> := PolynomialRing(QQ,2);
AttachSpec("spec");
SetVerbose("Theta",true);
CC<I> := ComplexField();
z := [0.85746943+0.30772689*I];
tau := [0.63443932+0.49919024*I];
Y := Imaginary(tau);
tau := Matrix(1,1,tau);
x := Real(z);
y := [Imaginary(el) : el in z];
Yinv := [1/el : el in Y];
//y0 := [Yinv*el : el in y];
Theta(z,tau);
