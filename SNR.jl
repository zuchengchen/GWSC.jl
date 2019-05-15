# power-law spectra of gravitational backgrounds
Ωgw(Ωβ, β) = f -> Ωβ * f^β


"""
SNR with `1yr` observation.
Note we have omit √2 for convenience.
"""
function SNR(Ωgw::Function, Ωeff::Function, fRange::Array)
    integral(f) = (Ωgw(f)/Ωeff(f))^2
    fMin = fRange[1]
    fMax = fRange[2]
    result = quadgk(integral, fMin, fMax, rtol=1e-3)
    result[1]^(1/2)
end


"""
SNR with `1yr` observation.
Note we have omit √2 for convenience.
"""
function SNR(det::Detector, Ωβ, β; fRange=nothing)
    
    Ωeff = det.Ωeff
    
    fMin = det.fMin
    fMax = det.fMax
    
#     if fRange == nothing
#         fRange = [det.fMin, det.fMax]
#     end
    
#     snr = SNR(Ωgw(Ωβ, β), Ωeff, fRange)
    
    if typeof(det) == LISA
        integral(f) = (Ωgw(Ωβ, β)(f)/Ωeff(f))^2
        result = quadgk(integral, fMin, fMax, rtol=1e-3)
        snr =result[1]^(1/2)
        snr *= sqrt(det.TObs * YEAR)
        
    elseif typeof(det) == PTA
        integral = f -> big(Ωgw(Ωβ, β)(f)/Ωeff(f))^2
        result = quadgk(integral, fMin, fMax, rtol=1e-3)
        snr =result[1]^(1/2)
        snr *= 0.85 * sqrt(2 * det.TObs * YEAR)
        
    else
        snr *= sqrt(2 * det.TObs * YEAR)
    end        
end

function getΩβ(det::Detector, β) 
    integral(f) = f^(2*β)/det.Ωeff(f)^2
    
    fMin = det.fMin
    fMax = det.fMax
    
    # Eq.(29)
    coeff1 = typeof(det)==LISA ? 1. : 2.
    coeff2 = typeof(det)==PTA ? 1/0.85 : 1.
    
    result = quadgk(integral, fMin, fMax, rtol=1e-3)
    Ωβ = coeff2 * det.ρThSGWB * (coeff1 * det.TObs * YEAR * result[1])^(-1/2)
    
    Ωβ
end

function ΩPI(det::Detector; fRange=nothing, nPoints=10^3, 
        ρth=nothing, βs=-8:0.1:8)
    
    if fRange != nothing
        fMin = fRange[1]
        fMax = fRange[2]
    else
        fMin, fMax = det.fPlotRange
    end
    
    logfMin = log10(fMin)
    logfMax = log10(fMax)
    
    fs = 10 .^ range(logfMin, logfMax, length=nPoints)
    
    if ρth==nothing
        ρth = det.ρThSGWB
    end      
    
    ΩgwLines = pmap(β -> Ωgw(getΩβ(det, β), β).(fs), βs)
    
    ΩPIs = [maximum(hcat(ΩgwLines...)'[:,i]) for i in 1:nPoints]
    
    fs, ΩPIs, ΩgwLines
end