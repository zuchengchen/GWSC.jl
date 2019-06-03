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
julia> plotΩPI(ipta, plotΩeffQ=true, plotΩPILinesQ=true, ΩPlotRange=(1e-15, 1e-8))
```
![](https://github.com/bingining/GWSC.jl/blob/master/test/pta.png)

More complete examples can be found at [examples](https://github.com/bingining/GWSC.jl/tree/master/examples)

## TODO List

- Add other detectors: BBO, DECIGO ...
- Add more comments, explanations and examples
- Implement the Bayesian approch for PTA
- Add h_c PI curve
- Probably write an article to explain the algorithms implemented in the package
- Do more tests
- Register as an offical julia package if possible

## References

* https://arxiv.org/abs/1310.5300
* https://arxiv.org/abs/1406.5199
* https://arxiv.org/abs/1803.01944

## Contributors

* **Zu-Cheng Chen** - [bingining](https://github.com/bingining/)

Please email the author with any bugs or requests. 

## License

This project is licensed under the GPL3 License - see the [LICENSE](LICENSE) file for details.
