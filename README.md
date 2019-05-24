# GWSC.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bingining.github.io/GWSC.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bingining.github.io/GWSC.jl/dev)

GWSC is the julia package to calculate and plot sensitivity curves for gravitational-wave detectors.

## Example

```julia
using GWSC

ipta = PTA(NP=36, σRMS=1e2, TObs=20.);
plotΩPI(ipta, plotΩeff=true, plotΩPILines=true, 
    ΩPlotRange=(1e-15, 1e-8))
```
## References
LISA sensitivity curve from https://arxiv.org/abs/1803.01944
