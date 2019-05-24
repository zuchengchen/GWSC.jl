<<<<<<< HEAD
# GWSC

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bingining.github.io/GWSC.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bingining.github.io/GWSC.jl/dev)
=======
# GWSC.jl

Sensitivity curves for gravitational-wave detectors.

=======
# GWSC.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bingining.github.io/GWSC.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bingining.github.io/GWSC.jl/dev)

GWSC is the julia package to calculate and plot sensitivity curves for gravitational-wave detectors.

## Install
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/bingining/GWSC.jl.git"))

## Example

```julia
using GWSC

ipta = PTA(NP=36, σRMS=1e2, TObs=20.);
plotΩPI(ipta, plotΩeff=true, plotΩPILines=true, 
    ΩPlotRange=(1e-15, 1e-8))
```
## References
>>>>>>> 5559996083fa71abd7e610227d2da4329b6a2499
LISA sensitivity curve from https://arxiv.org/abs/1803.01944
