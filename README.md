# GWSC.jl

Gravitational-Wave Sensitivity Curves.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bingining.github.io/GWSC.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bingining.github.io/GWSC.jl/dev)

GWSC is the julia package to calculate and plot sensitivity curves for gravitational-wave detectors.

## Installation

```julia
julia> using Pkg;
julia> Pkg.add(PackageSpec(url="https://github.com/bingining/GWSC.jl.git"))
```


## Usage

```julia
julia> using GWSC
julia> ipta = PTA(NP=36, σRMS=1e2, TObs=20.);
julia> plotΩPI(ipta, plotΩeff=true, plotΩPILines=true, ΩPlotRange=(1e-15, 1e-8))
```

More complete examples can be found at https://github.com/bingining/GWSC.jl/tree/master/examples

## TODO List
- Add other detectors: LIGO, KAGRA, BBO, DECIGO ...
- Add more comments, explanations and examples
- Implement the Bayesian approch for PTA
- Probably write an article to explain the algorithms implemented in the package
- Do more tests
- Register as an offical julia package if possible

## References

* https://arxiv.org/abs/1310.5300
* https://arxiv.org/abs/1406.5199
* https://arxiv.org/abs/1803.01944
