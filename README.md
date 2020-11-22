# QuadOsc.jl

This Julia package exports exactly one function `quadosc()` that is used to
integrate oscillatory functions to infinity. The algorithm works by integrating
the integrand between successive zeros using [QuadGK][1] and then summing the
resulting alternating series with a series acceleration, [described][2] in
[Cohen et al 2000][3].

Given an oscillatory function `fn(x)`,
```julia
julia> using QuadOsc
julia> a = 0.0
julia> I, E = quadosc(fn, a, Inf, n->fnzeros(n))
```
integrates `fn(x)` from `a` to infinity, and `fnzeros(n)` is the `n`-th zero of
`fn(x)`. That is, `fn(fnzeros(n)) == 0` for integer `n`.

The series acceleration is available via
```julia
julia> ak = @. (-1)^(1:5) / (1:5)
julia> sum_ak = QuadOsc.accel_cohen_villegas_zagier(ak)
```

Tests can be run by loading the package, entering package mode with pressing
`]`, and calling `test` on the package:
```julia
julia> using QuadOsc
pkg> test QuadOsc
```
'nuff said?


[1]: https://github.com/JuliaMath/QuadGK.jl
[2]: https://people.mpim-bonn.mpg.de/zagier/files/exp-math-9/fulltext.pdf
[3]: https://doi.org/10.1080/10586458.2000.10504632
