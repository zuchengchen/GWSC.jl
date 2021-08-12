# Set ASD file for different detectors
function ASDFile(name::String)
    if name == "LIGO_O1"
        ASDFile = "LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt"
    elseif name == "LIGO_O2"
        ASDFile = "LIGO-P1200087-v18-aLIGO_MID_LOW.txt"
    elseif name == "LIGO_O3"
        ASDFile = "LIGO-P1200087-v18-aLIGO_MID_HIGH.txt"
    elseif name == "LIGO_O5"
        ASDFile = "LIGO-P1200087-v18-aLIGO_LATE_HIGH.txt"
    elseif name == "KAGRA_Design"
        ASDFile = "LIGO-T1600593-v1-KAGRA_Design.txt"
    elseif name == "LIGO_Design"
        ASDFile = "LIGO-P1200087-v18-aLIGO_DESIGN.txt"
    elseif name == "ET"
#         ASDFile = "LIGO-P1600143-v18-ET_D.txt"
        ASDFile = "ET-0000A-18_ETDSensitivityCurveTxtFile.txt"
    elseif name == "CE"
        ASDFile = "LIGO-P1600143-v18-CE.txt"
    end
    ASDFile
end

# interpolate two arrays as a function using log-log scales internal
# and return to the normal scale
function interLogxLogyIn(xs, ys; out=Inf)
    logxs = log10.(xs)
    logys = log10.(ys)
    
    logInter0 = interpolate((logxs,), logys, Gridded(Linear()))
    logInter = extrapolate(logInter0, out)
        
    # convert back to normal scale
    f -> 10^logInter(log10(f))
end

# interpolate two arrays as a function using log-log scales internal
# and return to the normal scale
function interLogxIn(xs, ys; out=Inf)
    logxs = log10.(xs)
    
    logInter0 = interpolate((logxs,), ys, Gridded(Linear()))
    logInter = extrapolate(logInter0, out)
        
    # convert back to normal scale
    f -> logInter(log10(f))
end

struct LIGO <: Detector
    
    name::String # detector name
    TObs::Float # observation time (yr)
    NDet::Int # number of detectors
    
    R::Function # transfer function
    
    fMin::Float # minimum frequency (Hz) 
    fMax::Float # maximum frequency (Hz)
    fRef::Float # reference frequency (Hz)    
    fPlotRange::Tuple{Float,Float} # frequency range to show (Hz, Hz)
    ΩPlotRange::Tuple{Float,Float} # energy density range to show
    hPlotRange::Tuple{Float,Float} # strain range to show
    
    Pn::Function # Power Spectral Density without confusion (WC) noise (Hz^-1)
    Sn::Function # Sensitivity curve of LISA 
    # without the galactic binary confusion noise
    
    Ωn::Function # fractional energy density of noise
    Ωeff::Function # fractional effective energy density of noise
    
    ρThSGWB::Float # threshold SNR for deteting a gravitational background
    
    ASD
    
    function LIGO(;
            name="LIGO_O1",
            TObs=1.,
            NDet=2, 
            fRef=25., 
            fPlotRange::Tuple{Float,Float}=(1e1, 1e4), 
            ΩPlotRange::Tuple{Float,Float}=(1e-14, 1e-6), 
            hPlotRange::Tuple{Float,Float}=(1e-21, 1e-16),
            ρThSGWB=1.0)
        
        # read in the data interpolate ASD 
        ASD_file = joinpath(sensitivity_path, ASDFile(name))
        println("using $ASD_file")
        ASDData = readdlm(ASD_file) 
        
        fsASD = ASDData[:, 1] # frequency
        ASDs = ASDData[:, 2] # ASD        
        fMinASD, fMaxASD = fsASD[1], fsASD[end]
        fPlotRange = (fMinASD, fMaxASD)
               
        # convert back to normal scale
        ASD = interLogxLogyIn(fsASD, ASDs)
        
        Sn(f) = ASD(f)^2
        
        R(f) = 1 # This is true for LIGO-type interferometer
        
        Pn(f) = Sn(f) # because R(f)=1
        
        # interpolate the normalized overlap function
        overlap_file = joinpath(sensitivity_path, "H1L1_orf.dat")
        γData = readdlm(overlap_file)
        fsγ = γData[:, 1]
        γs = γData[:, 2]
        
        fMinγ, fMaxγ = fsγ[1], fsγ[end]        
        fMin = max(fMinASD, fMinγ)
        fMax = min(fMaxASD, fMaxγ)
        
        γ = interLogxIn(fsγ, γs)
        
        Γ(f) = 1/5 * γ(f)
        
        Ωn(f) = (2*π^2/3/H0^2) * f^3 * Sn(f) # eq.(3)
        Ωeff(f) = sqrt((2/NDet/(NDet-1))) * Ωn(f)/Γ(f)
        
        new(name, TObs, NDet, R, 
            fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange, 
            Pn, Sn, Ωn, Ωeff, ρThSGWB, ASD)
    end       
end


# TODO: add more explanation
function SNR(det::LIGO, Ωgw::Function)
    
    fMin, fMax, fRef = det.fMin, det.fMax, det.fRef 
    
    T = det.TObs * YEAR # (s)
    
    integral(f) = (Ωgw(f)/det.Ωeff(f))^2
    
    logfMin = log10(det.fMin)
    logfMax = log10(det.fMax)
    f1, f2, f3, f4, f5, f6 = 10 .^ range(logfMin, logfMax, length=6)
    
    int0 = quadgk(integral, f1, f2, f3, f4, f5, f6, rtol=1e-6)[1]
        
    snr = sqrt(2T*int0)
end