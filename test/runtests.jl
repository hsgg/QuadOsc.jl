#!/usr/bin/env julia

using Test
using QuadOsc
using QuadGK
using SpecialFunctions


@testset "SeriesAccelerations Cohen-Villegas-Zagier" begin
	# Check first few coefficients as given in Cohen et al 2000
	ak = [1]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A == 2/3

	ak = [1, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ 16/17
	ak = [0, -1]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ -8/17

	ak = [1, 0, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ 98/99
	ak = [0, -1, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ -80/99
	ak = [0, 0, 1]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ 32/99

	ak = [1, 0, 0, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ 576/577
	ak = [0, -1, 0, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ -544/577
	ak = [0, 0, 1, 0]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ 384/577
	ak = [0, 0, 0, -1]
	A = QuadOsc.accel_cohen_villegas_zagier(ak)
	@test A ≈ -128/577

	# Ex 1 in Cohen et al 2000
	# Note: we actually disagree with the answer given in the paper. What's
	# going on? Since all other tests pass, we probably solving a different
	# problem than them.
	fn(n) = (-1)^(n+1) * loggamma(1 + 1 / (2*n - 1))
	println("n an\tsum(ak)\tA")
	for n=1:20
		ak = @. fn(1:n)
		lnA = QuadOsc.accel_cohen_villegas_zagier(ak)
		A = exp(lnA)
		Adirect = exp(sum(ak))
		println("$n $(ak[end])\t$Adirect\t$A")
	end

	# Ex 3 in Cohen et al 2000
	riemann_zeta(s; n=20) = begin
		k = 0:n-1
		ak = @. (-1)^k / (k + 1)^s / (1 - 2^(1 - s))
		return QuadOsc.accel_cohen_villegas_zagier(ak)
	end
	@test riemann_zeta(1/2) ≈ -1.4603545088
	@test riemann_zeta(-1+im) ≈ 0.0168761517 - 0.1141564804*im
end


function test_quadosc()
	integrand(x) = x == 0.0 ? 0.0 : x * besselj0(x) / (1 + x^2)
	nterm = 15
	a = 0.0
	b = 0.0
	sum = 0.0
	ak = Float64[]
	ans = 0.0
	correct = 0.4210244382
	println("N\tSum(direct)\tSum(Levin)")
	for n in 0:nterm
		b += pi
		s, E = quadgk(integrand, a, b)
		push!(ak, s)
		a = b
		sum += s
		ans = QuadOsc.accel_cohen_villegas_zagier(ak)
		println("$n\t$sum\t$ans")
	end
	directans, E = quadosc(integrand, 0, Inf, n -> n * pi)
	println("ans          = $ans")
	println("directans    = $directans +- $E")
	println("correct + dI = $(correct) + $(directans - correct)")
	@test abs(ans - correct) < 1e-10
	@test abs(ans - directans) < 1e-10
end


function test2()
	thisans, E = quadosc(x -> sphericalbesselj(0, x), 0, Inf, n -> n * pi)
	println("int sphbes0: $thisans +- $E =? $(pi/2)")
	println("I          : $(pi/2) + $(thisans - pi/2)")
	@test abs(thisans - pi/2) < E
end


function test_int_Jν()
	for ℓ=0:100
		Jl_zeros(n,ℓ) = (n + ℓ/2) * π  # just an estimate
		I, E = quadosc(t -> besselj(ℓ+0.5, t), 0, Inf, n -> Jl_zeros(n,ℓ))
		@show ℓ,I,E
		@test I ≈ 1
		@test abs(I - 1) <= E
	end
end


function typetest()
	integrand(x::Float64) = (trunc(Int, x) % 2 == 1) ? -1/x : 1/x
	@time quadosc(integrand, 1, Inf, n -> n)
	@time quadosc(integrand, 1, Inf, n -> n)
end



@testset "QuadOsc" begin
	test_quadosc()
	test2()
	test_int_Jν()
	typetest()
end
