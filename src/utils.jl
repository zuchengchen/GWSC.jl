using DelimitedFiles:readdlm
using Interpolations
using Unitful, UnitfulAstro
using PyPlot
using PhysicalConstants.CODATA2014: c_0
using QuadGK
using Roots
# using HypergeometricFunctions

abstract type Detector end

# spedd of light in vaccum (m/s)
const c0 = ustrip(float(c_0))

# Julian year in seconds
const YEAR = ustrip(uconvert(Unitful.s, 1*u"yr"))

# Day in seconds
const DAY = YEAR/365.25

# Hubble constant (s^-1)
H0 = ustrip(uconvert(u"s"^-1, 67.66 * u"km" * u"Mpc"^-1 * u"s"^-1))

# Astronomical unit (meters)
# const AU = ustrip(uconvert(Unitful.m, 1*u"AU"))

const Float = Float64