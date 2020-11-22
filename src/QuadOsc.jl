module QuadOsc

export quadosc

include("SeriesAccelerations.jl")
using .SeriesAccelerations
using QuadGK


@doc raw"""
    quadosc(fn, a, Inf, fnzeros; ...)

Integrate the function `fn(x)` from `a` to `Inf`. The function `fnzeros(n)`
takes an integer `n` and is such that `fn(fnzeros(n)) == 0`. The algorithm
works by integrating between successive zeros, and accelerating the alternating
series.

The argument `pren` is the number of intervals to integrate before applying the
series acceleration.

`atol` and `rtol` specify the absolute and relative tolerances for determining
convergence.

`order` is passed on to `quadgk()` of the
[QuadGK](https://github.com/JuliaMath/QuadGK.jl) package.

`nconvergences` is the number of iterations before convergence is declared.

See `?QuadOsc.accel_cohen_villegas_zagier` for details on the series
acceleration.
"""
function quadosc(f::Function, a::Number, b::Number, zerosf::Function; pren=2,
		 atol=zero(Float64), rtol=sqrt(eps(Float64)), order=7,
		 nconvergences=ceil(Int,-1.31*log10(rtol)))
    @assert b == Inf
    T = Float64

    i1 = findfirst(n -> zerosf(n) - a >= 0, 1:typemax(Int))
    z1 = zerosf(i1)
    Ipre, Epre = quadgk(f, a, z1; atol=atol, rtol=rtol, order=order)

    z0 = z1
    for i=1:pren
        i1 += 1
        z1 = zerosf(i1)
        I, E = quadgk(f, z0, z1; atol=atol, rtol=rtol, order=order)
        Ipre += I
        Epre += E
        z0 = z1
    end

    I = T(0)
    oldI = T(0)
    ak = T[]
    ek = T[]
    while nconvergences > 0
        i1 += 1
        z1 = zerosf(i1)
        I, E = quadgk(f, z0, z1; atol=atol, rtol=rtol)
        push!(ak, I)
        push!(ek, E)
        z0 = z1
        I = accel_cohen_villegas_zagier(ak)

        adiff = abs(I - oldI)
        rdiff = adiff * 2 / abs(I + oldI)
        if adiff <= atol || rdiff <= rtol
            nconvergences -= 1
        end
        oldI = I
    end
    I = Ipre + I
    E = Epre + accel_cohen_villegas_zagier(ek)
    return I, E
end


end # module


# vim: set sw=4 et sts=4 :
