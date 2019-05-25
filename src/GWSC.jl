__precompile__()

module GWSC

using DelimitedFiles:readdlm
using Interpolations
using Unitful, UnitfulAstro
using PyPlot
using PhysicalConstants.CODATA2014: c_0
using QuadGK
using Roots

export
    YEAR,
    DAY,
    c0,
    H0,
    Float,
    Detector,
    PTA,
    LISA,
    plotΩPI,
    plotCharacteristicStrain,
    plotSpectralDensity,
    ΩPI

include("utils.jl")
include("LISA/LISA.jl")
include("PTA/PTA.jl")
include("SNR.jl")
include("plot.jl")

end # module