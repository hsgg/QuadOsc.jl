module QuadOsc

export quadosc

include("SeriesAccelerations.jl") using .SeriesAccelerations using QuadGK


# Integrate from 'a' to 'b'.
function quadosc(f::Function, a::Number, b::Number, zeros::Function; pren=2,
		 atol=zero(Float64), rtol=sqrt(eps(Float64)), order=7,
		 nconvergences=ceil(Int,-1.31*log10(rtol))-pren) @assert b ==
	Inf T = Float64

	i1 = findfirst(n -> zeros(n) - a >= 0, 1:typemax(Int)) presum, Err =
	quadgk(f, a, zeros(i1); atol=atol, rtol=rtol, order=order)

	i1max = i1 + pren while i1 < i1max I, E = quadgk(f, zeros(i1),
						 zeros(i1+1); atol=atol,
					 rtol=rtol, order=order) presum += I
	Err += E i1 += 1 end

	z0 = zeros(i1) i = 1 ansI = T(0) ansE = T(0) oldansI = T(0) ak = T[] ek
	= T[] while nconvergences > 0 z1 = zeros(i+i1) I, E = quadgk(f, z0, z1;
								     atol=atol,
		rtol=rtol) push!(ak, I) push!(ek, E) ansI =
		accel_cohen_villegas_zagier(ak) ansE =
		accel_cohen_villegas_zagier(ek)

		adiff = abs(ansI - oldansI) rdiff = adiff * 2 / abs(ansI +
								    oldansI) if
		adiff <= atol || rdiff <= rtol nconvergences -= 1
			#@show i,adiff,rdiff,ansI+presum
		end oldansI = ansI

		z0 = z1 i += 1 end I = ansI + presum E = Err + ansE return I, E
	end


end # module
