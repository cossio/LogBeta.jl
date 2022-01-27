# LogBeta Julia package

![](https://github.com/cossio/LogBeta.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/cossio/LogBeta.jl/branch/master/graph/badge.svg?token=YHVMRQOGXE)](https://codecov.io/gh/cossio/LogBeta.jl)

This package doesn't export any symbols.
The main function is `LogBeta.logbeta`, which computes the log of the incomplete Beta function in a numerically stable way:

```
LogBeta.logbeta(a, b) = log(beta(a, b))
LogBeta.logbeta(a, b, x) = log(beta(a, b, x))
LogBeta.logbeta(a, b, x1, x2) = log(beta(a, b, x2) - beta(a, b, x1))
```

This package is not registered.
Install with:

```julia
using Pkg
Pkg.add(url="https://github.com/cossio/LogBeta.jl")
```