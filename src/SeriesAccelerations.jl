@doc raw"""
    SeriesAccelerations

This modules provides an algorithm for calculating infinite series efficiently.
Currently, only one algorithm for an alternating series is provided, see
`accel_cohen_villegas_zagier()` for more details.
"""
module SeriesAccelerations


export accel_cohen_villegas_zagier


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



end


# vim: set sw=4 et sts=4 :
