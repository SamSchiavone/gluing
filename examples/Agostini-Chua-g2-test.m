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
SetDebugOnError(true);
prec := 300;
SetDefaultRealFieldPrecision(prec);
CC<I> := ComplexField(prec);
z := [CC | 0.81149300+0.27027128*I, 0.77132834+0.26619567*I];
tau := Matrix(2,2,[CC | 0.71106237+1.20021283*I, 0.57281731+0.89762698*I, 0.57281731+0.89762698*I, 0.22079146+0.68617488*I]);
Y := Imaginary(tau);
x := Real(z);
y := [Imaginary(el) : el in z];
Yinv := Inverse(Y);
//y0 := Inverse(Y)*y;
/*
print "Testing theta";
print "Should be 1.700489+1.03657im";
Theta(z,tau);
*/
print "Testing theta with dz = [1,0]";
print "Should be -13.73322+34.5327im";
Theta(z,tau : dz := [[1,0]]);
/*
print "Testing theta with char=[[1,0],[1,1]])";
print "Should be -5.91176+4.90798im";
Theta(z, tau : char := [[1,0], [1,1]]);
*/
Theta(z,tau : dtau := [[1,0]]);

// see https://github.com/chualynn/Theta.jl/blob/master/test/theta_test.jl
z := [0.04134584+0.40910551*I, 0.20972589+0.90269823*I, 0.39996195+0.42432923*I, 0.73063375+0.49945621*I];
tau := Matrix(CC,4,4,[0.95870734+0.73587725*I, 0.22092477+0.76863646*I, 0.53877459+0.87577267*I, 0.68177023+0.867436*I, 0.22092477+0.76863646*I, 0.98812562+1.79674905*I, 0.54859032+1.10626215*I, 0.63310305+1.30158981*I, 0.53877459+0.87577267*I, 0.54859032+1.10626215*I, 0.50173043+1.27729044*I, 0.49163557+1.33147334*I, 0.68177023+0.867436*I, 0.63310305+1.30158981*I, 0.49163557+1.33147334*I, 0.35312207+1.60745975*I]);
Theta(z,tau);

Theta(z, tau : char:=[[0,1,0,1],[0,1,0,0]], dz:=[[0,0,1,0]]);
Theta(z, tau : char:=[[0,1,0,1],[0,1,0,0]], dz:=[[0,0,1,0],[1,0,0,0]]);
Theta(z, tau : char:=[[0,1,0,1],[0,1,0,0]], dz := [[1,0,0,0]], dtau:=[[1,2]]);
Theta(z, tau : dtau:=[[1,2]]);
