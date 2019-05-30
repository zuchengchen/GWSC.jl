# power-law spectra of gravitational backgrounds
Ωgw(Ωβ, β, fRef) = f -> Ωβ * (f/fRef)^β

# TODO: add more explanation
function getΩβ(det::Detector, β) 
    
    ρTh = det.ρThSGWB
    
    if typeof(det)==PTA
        equation = Ωβ -> SNR(det, Ωgw(Ωβ, β, det.fRef)) - ρTh
        return find_zero(equation, (1e-500, 1e100))
    end
    
    integral(f) = (f/det.fRef)^(2*β)/det.Ωeff(f)^2
    
    fMin = det.fMin
    fMax = det.fMax
    
    # Eq.(36) & Eq.(29)
    result = quadgk(integral, fMin, fMax, rtol=1e-4)
    Ωβ = ρTh * (det.TObs*YEAR * result[1])^(-1/2)
end

# TODO: add more explanation
function ΩPI(det::Detector; 
        fRange=nothing, 
        nPoints=10^3, 
        ρth=nothing,
        file=nothing)
    
    if typeof(det)==PTA
        βs=[-150:20:-10; -10:2:0; 0:0.5:10]
    else
        βs=-10:0.5:10
    end
    
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
    
    if file != nothing
        backup(fs, ΩPIs, file)
    end
    
    fs, ΩPIs, ΩgwLines
end
