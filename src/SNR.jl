# power-law spectra of gravitational backgrounds
Ωgw(Ωβ, β, fRef) = f -> Ωβ * (f/fRef)^β

# TODO: add more explanation
function getΩβ(det::Detector, β) 
    
    if typeof(det)==PTA
        return find_zero(Ωβ -> SNR(det, β, Ωβ) - det.ρThSGWB, (1e-500, 1e100))
    end
    
    integral(f) = (f/det.fRef)^(2*β)/det.Ωeff(f)^2
    
    fMin = det.fMin
    fMax = det.fMax
    
    # Eq.(36) & Eq.(29)
    result = quadgk(integral, fMin, fMax, rtol=1e-3)
    Ωβ = det.ρThSGWB * (det.TObs*YEAR * result[1])^(-1/2)
end

# TODO: add more explanation
function ΩPI(det::Detector; 
        fRange=nothing, 
        nPoints=10^3, 
        ρth=nothing,
        file=nothing)
    
    if typeof(det)==PTA
        βs=[-200:20:-10; -10:2:0; 0:0.5:10]
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

# TODO: add more explanation
function plotΩPI(det::Detector; 
        figure_file=nothing, 
        plotΩPI=true, 
        plotΩeff=false, 
        plotΩPILines=false,
        fPlotRange=nothing,
        ΩPlotRange=nothing)
    
    fs, ΩPIResult, ΩPILines  = ΩPI(det)
    Ωeffs = det.Ωeff.(fs)
    
    if fPlotRange == nothing
        fPlotRange = det.fPlotRange
    end
    
    if ΩPlotRange == nothing
        ΩPlotRange = det.ΩPlotRange
    end
    
    xlim(fPlotRange)
    ylim(ΩPlotRange)    
    
    grid("on", which="both", linestyle="--")
    
    if plotΩPI != false
        loglog(fs, ΩPIResult, "black", label=L"Ω$_{\mathrm{PI}}$")        
    end
    
    if plotΩeff==true
        loglog(fs, Ωeffs, "b", label=L"Ω$_{\mathrm{eff}}$")
        if typeof(det)==PTA 
            vlines(det.fMin, det.Ωeff(det.fMin), ΩPlotRange[2], "b")
        end
    end
    
    if plotΩPILines==true 
        [loglog(fs, ΩPILines[i], "--", color="gray") for i in 1:length(ΩPILines)]
    end
    
    legend()
    show()
end