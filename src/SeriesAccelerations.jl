@doc raw"""
    SeriesAccelerations

This modules provides an algorithm for calculating infinite series efficiently.
Currently, two algorithms are provided. `accel_cohen_villegas_zagier()` for an
alternating series, and `accel_wynn_eps()` for more general series.
"""
module SeriesAccelerations


export accel_cohen_villegas_zagier
export accel_wynn_eps


################# Cohen, Villegas, Zagier 2000 algorithm #####################

@doc raw"""
    accel_cohen_villegas_zagier(ak)

Estimate the sum over the array `ak`. It is assumed without check that the
elements in the array have alternating sign and the series converges.

This algorithm is
[presented](https://people.mpim-bonn.mpg.de/zagier/files/exp-math-9/fulltext.pdf)
in [Cohen et al 2000]
(https://www.tandfonline.com/doi/abs/10.1080/10586458.2000.10504632).
"""
function accel_cohen_villegas_zagier(ak)
    T = eltype(ak)
    n = length(ak)
    d = (3 + √T(8))^n
    d = (d + 1/d) / 2
    b = -T(1)
    c = -d
    s = T(0)
    for k=0:n-1
        c = b - c
        s += (-1)^k * c * ak[k+1]
        b = (k + n) * (k - n) * b / ((k + 1/T(2)) * (k + 1))
    end
    return s / d
end


################ Wynn 1956 epsilon-algorithm #################################

function next_j!(c2, ϵ, n)
    c1 = (c2 & 1) + 1  # c1 is the other column (row, actually)
    for i=1:n
        ϵ[c1,i] = ϵ[c1,i+1] + 1 / (ϵ[c2,i+1] - ϵ[c2,i])
    end
    return c1
end


function accel_wynn_eps!(ϵ)
    N = size(ϵ, 2)
    c = 2
    for j=1:N-1
        c = next_j!(c, ϵ, N - j)
    end
    return ϵ[c,1]
end


@doc raw"""
    accel_wynn_eps(ak)

Estimate the sum over the array `ak`. It is assumed without check that the
series converges.

The algorithm used is [Wynn's
epsilon-algorithm](https://www.ams.org/journals/mcom/1956-10-054/S0025-5718-1956-0084056-6/S0025-5718-1956-0084056-6.pdf)
presented in [Wynn 1956](https://www.jstor.org/stable/2002183).
"""
function accel_wynn_eps(ak)
    ΔN = Int(iseven(length(ak)))
    N = length(ak) - ΔN
    ϵ = Array{eltype(ak)}(undef, 2, N)
    for n=1:N
        ϵ[1,n] = zero(eltype(ak))
        ϵ[2,n] = sum(ak[1:n+ΔN])
    end
    return accel_wynn_eps!(ϵ)
end


end


# vim: set sw=4 et sts=4 :
