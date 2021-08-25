struct DECIGO <: SpaceDetector
    TObs::Float # observation time of LISA (yr)
    LArm::Float # Arm length of for every arm of LISA (meter)
    NC::Int # Number of data channels
    
    fStar::Float # transfer frequency (Hz)
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
    
    function DECIGO(;
            TObs::Float=4., 
            LArm::Float=3e9, 
            NC::Int=2, 
            fRef::Float=1e-3, 
            fPlotRange::Tuple{Float,Float}=(1e-5, 1e0), 
            ΩPlotRange::Tuple{Float,Float}=(1e-14, 1e-6), 
            hPlotRange::Tuple{Float,Float}=(1e-21, 1e-16),
            ρThSGWB=10.0)
        
        fStar = c0/(2*π*LArm)
        R(f) = NC * 3.0/20.0/(1.0 + 0.6*(f/fStar)^2)
        file = joinpath(sensitivity_path, "DECIGO.txt")
        data = readdlm(file) # read in the data
        fMin = data[1, 1]
        fMax = data[end, 1]
        function Sn(f)
            
            x = data[:, 1]
            y = data[:, 2]
    
            logx = log10.(x)
            logy = log10.(y)
    
            # interpolate the transfer function with log10-log10 scale
            logInter0 = interpolate((logx,), logy, Gridded(Linear()))
            logInter = extrapolate(logInter0, Inf)  
            # convert back to normal scale
            10^logInter(log10(f))
        end
        Pn(f) = Sn(f) * R(f)
        
        Ωn(f) = (2*π^2/3/H0^2) * f^3 * Sn(f)
        Ωeff(f) = (2/NC/(NC-1)) * Ωn(f)
        
        new(TObs, LArm, NC, fStar, R, 
            fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange, 
            Pn, Sn, Ωn, Ωeff, ρThSGWB)
    end       
end