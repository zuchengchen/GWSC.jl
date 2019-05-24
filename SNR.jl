# power-law spectra of gravitational backgrounds
Ωgw(Ωβ, β, fRef) = f -> Ωβ * (f/fRef)^β


# """
# SNR with `1yr` observation.
# Note we have omit √2 for convenience.
# """
# function SNR(Ωgw::Function, Ωeff::Function, fRange::Array)
#     integral(f) = (Ωgw(f)/Ωeff(f))^2
#     fMin = fRange[1]
#     fMax = fRange[2]
#     result = quadgk(integral, fMin, fMax, rtol=1e-3)
#     result[1]^(1/2)
# end


function SNR(pta::PTA, β, Ωβ)
    # increase precision
    β, Ωβ = big(β), big(Ωβ)
    
    Δt = pta.Δt # (s)
    σ = pta.σRMS*1e-9 # (s)
    A = 16*π^4*Δt*σ^2/H0^2
    b = β-5
    
    fMin, fMax, fRef = pta.fMin, pta.fMax, pta.fRef 
    
    N = pta.NP
    
    T = pta.TObs * YEAR # (s)
    ζ2 = 1/48. # square of the average of Hellings and Downs factor
    
    # B*f^b/A
    Bfb(f) = Ωβ*(f/fRef)^β*f^(-5.)/A
    
    int(f) = f/b * (b + 1/(1 + Bfb(f)) 
        - (1+b)*_₂F₁(1, 1/b, 1+1/b, - Bfb(f)))
    
    snr2 = 2T*N*(N-1)/2*ζ2*(int(fMax) - int(fMin))
    snr2^(1/2)
end

function getΩβ(pta::PTA, β)
    find_zero(Ωβ -> SNR(pta, β, Ωβ) - pta.ρThSGWB, (1e-40, 1e40))
end


# """
# SNR with `1yr` observation.
# Note we have omit √2 for convenience.
# """
# function SNR(det::Detector, Ωβ, β; fRange=nothing)
    
#     Ωeff = det.Ωeff
    
#     fMin = det.fMin
#     fMax = det.fMax
    
# #     if fRange == nothing
# #         fRange = [det.fMin, det.fMax]
# #     end
    
# #     snr = SNR(Ωgw(Ωβ, β), Ωeff, fRange)
    
#     if typeof(det) == LISA
#         integral(f) = (Ωgw(Ωβ, β)(f)/Ωeff(f))^2
#         result = quadgk(integral, fMin, fMax, rtol=1e-3)
#         snr =result[1]^(1/2)
#         snr *= sqrt(det.TObs * YEAR)
        
#     elseif typeof(det) == PTA
#         integral = f -> big(Ωgw(Ωβ, β)(f)/Ωeff(f))^2
#         result = quadgk(integral, fMin, fMax, rtol=1e-3)
#         snr =result[1]^(1/2)
#         snr *= 0.85 * sqrt(2 * det.TObs * YEAR)
        
#     else
#         snr *= sqrt(2 * det.TObs * YEAR)
#     end        
# end

# function getΩβ(det::Detector, β) 
#     integral(f) = f^(2*β)/det.Ωeff(f)^2
    
#     fMin = det.fMin
#     fMax = det.fMax
    
#     # Eq.(29)
#     coeff1 = typeof(det)==LISA ? 1. : 2.
#     coeff2 = typeof(det)==PTA ? 1 : 1.
    
#     result = quadgk(integral, fMin, fMax, rtol=1e-3)
#     Ωβ = coeff2 * det.ρThSGWB * (coeff1 * det.TObs * YEAR * result[1])^(-1/2)
    
#     Ωβ
# end

# function ΩPI(det::Detector; fRange=nothing, nPoints=10^3, 
#         ρth=nothing, βs=-8:0.5:8)
    
#     if fRange != nothing
#         fMin = fRange[1]
#         fMax = fRange[2]
#     else
#         fMin, fMax = det.fPlotRange
#     end
    
#     logfMin = log10(fMin)
#     logfMax = log10(fMax)
    
#     fs = 10 .^ range(logfMin, logfMax, length=nPoints)
    
#     if ρth==nothing
#         ρth = det.ρThSGWB
#     end      
    
#     ΩgwLines = pmap(β -> Ωgw(getΩβ(det, β), β).(fs), βs)
    
#     ΩPIs = [maximum(hcat(ΩgwLines...)'[:,i]) for i in 1:nPoints]
    
#     fs, ΩPIs, ΩgwLines
# end

function ΩPI(det::PTA; fRange=nothing, nPoints=10^3, 
        ρth=nothing, βs=-20:1:20)
    
    if fRange != nothing
        fMin = fRange[1]
        fMax = fRange[2]
    else
        fMin, fMax = det.fPlotRange
    end
    
    logfMin = log10(fMin)
    logfMax = log10(fMax)
    
    fs = 10 .^ range(logfMin, logfMax, length=nPoints)
    fRef = det.fRef
    
    if ρth==nothing
        ρth = det.ρThSGWB
    end      
    
    ΩgwLines = map(β -> Ωgw(getΩβ(det, β), β, fRef).(fs), βs)
    
    ΩPIs = [maximum(hcat(ΩgwLines...)'[:,i]) for i in 1:nPoints]
    
    fs, ΩPIs, ΩgwLines
end