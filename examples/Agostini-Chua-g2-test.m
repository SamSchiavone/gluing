/*
https://github.com/chualynn/Theta.jl/blob/master/test/theta_test.jl
@testset "Genus 2" begin
    z = [0.81149300+0.27027128im; 0.77132834+0.26619567im];
    τ = [0.71106237+1.20021283im 0.57281731+0.89762698im; 0.57281731+0.89762698im 0.22079146+0.68617488im];
    Y = imag(τ);
    x = real(z);
    y = convert(Array{Float64}, imag(z));
    Yinv = inv(Y);
    y0 = Yinv*y
    R = RiemannMatrix(τ, siegel=false);
    @test Theta.oscillatory_part(R, x, y0, [0,0]) ≈ 0.5894+0.3593im atol=ϵ
    @test theta(z, R) ≈ 1.700489+1.03657im atol=ϵ
    @test theta(z, R, derivs=[[1,0]]) ≈ -13.73322+34.5327im atol=ϵ
    @test theta(z, R, derivs=[[0,1]]) ≈ 14.71768-51.3656im atol=ϵ
    @test theta(z, R, derivs=[[1,0], [1,0]]) ≈ -1859.654+898.206im atol=ϵ
    @test theta(z, R, derivs=[[1,0], [0,1]]) ≈ 2533.2922-1397.9940im atol=ϵ
    @test theta(z, R, char=[[1,0],[1,1]]) ≈ -5.91176+4.90798im atol=ϵ
    @test theta(z, R, char=[[1,0],[0,1]]) ≈ -0.68197+0.13723im atol=ϵ
    @test theta(z, R, char=[[1,0],[1,1]], derivs=[[1,0]]) ≈ -33.55188-104.14870im atol=ϵ
    @test theta(z, R, char=[[1,0],[0,1]], derivs=[[0,1]]) ≈ 20.42986+95.749276im atol=ϵ
    @test symplectic_transform(siegel_transform(τ)[1], τ) ≈ siegel_transform(τ)[2]
end
*/

//R<x,y> := PolynomialRing(QQ,2);
AttachSpec("spec");
SetVerbose("Theta",true);
CC<I> := ComplexField();
z := [0.81149300+0.27027128*I, 0.77132834+0.26619567*I];
tau := Matrix(2,2,[0.71106237+1.20021283*I, 0.57281731+0.89762698*I, 0.57281731+0.89762698*I, 0.22079146+0.68617488*I]);
Y := Imaginary(tau);
x := Real(z);
y := [Imaginary(el) : el in z];
Yinv := Inverse(Y);
//y0 := Inverse(Y)*y;
Theta(z,tau);
