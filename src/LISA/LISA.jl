# TODO: extend this module

"LISA modulation frequency"
fm = 3.168753575e-8  


""" 
    RLISA(fStar::Float=19.09e-3, NC::Int=2; approQ=false)

Calculate the transfer function (sky and polarization averaged).

`approQ=false`: whether to approximate the transfer function or not;

`fStar`: transfer frequency;

`NC`: number of data channels for LISA.
"""
function RLISA(fStar::Float=19.09e-3, NC::Int=2; approRQ=false)  
    
    if approRQ == true        
        # return an approximated version of transfer function
        return f::Float -> NC * 3.0/20.0/(1.0 + 0.6*(f/fStar)^2)
    end
    
    file = joinpath(sensitivity_path, "RLISA.txt")
    RData = readdlm() # read in the data
        
    f = RData[:, 1]*fStar # convert to frequency
    R = RData[:, 2]*NC # response gets improved by more data channels
    
    # convert f and R to log10 scale
    logf = log10.(f)
    logR = log10.(R)
    
    # interpolate the transfer function with log10-log10 scale
    logRInter0 = interpolate((logf,), logR, Gridded(Linear()))
    logRInter = extrapolate(logRInter0, Inf)  
    
    # convert back to normal scale
    f::Float -> 10^logRInter(log10(f))
end


"""
SnC(f::Float, TObs::Float=4.0, NC::Int=2)

An estimation of the galactic binary confusion noise for
Tobs = {0.5 yr, 1 yr, 2 yr, 4yr}.
Enter Tobs as a year or fraction of a year.

`f`: frequency (Hz);

`TObs`: observation time of LISA (yr);

`NC`: number of data channels for LISA.
"""
function getSnC(f::Float, TObs::Float=4.0, NC::Int=2)
    
    # Fix the parameters of the confusion noise fit
    if (TObs < 0.75)
        α = 0.133
        β = 243.
        κ = 482.
        γ = 917.
        fk = 2.58e-3  
    elseif (0.75 < TObs && TObs < 1.5)
        α = 0.171
        β = 292.
        κ = 1020.
        γ = 1680.
        fk = 2.15e-3 
    elseif (1.5 < TObs && TObs < 3.0) 
        α = 0.165
        β = 299.
        κ = 611.
        γ = 1340.
        fk = 1.73e-3 
    else
        α = 0.138
        β = -221.
        κ = 521.
        γ = 1680.
        fk = 1.13e-3 
    end
    
    A = 1.8e-44/NC
    
    # Eq.(14)
    A*f^(-7.0/3.) * exp(-f^α + β*f*sin(κ*f)) * (1.0 + tanh(γ*(fk-f)))
end




"""
    LISA type

`TObs::Float`: observation time of LISA (yr);

`LArm`: Arm length of for every arm of LISA (meter);

`NC`: number of data channels for LISA;

`approQ=false`: whether to approximate the transfer function or not;

`fStar::Float`: transfer frequency (Hz)

`R::Function`: transfer function

`Pn::Function`: Power Spectral Density without confusion (WC) noise

`Sn::Function`: Sensitivity curve of LISA 
    without the galactic binary confusion noise
    
`Pn_WC::Function`: Power Spectral Density with confusion (WC) noise
    
`Sn_WC::Function`: Sensitivity curve of LISA 
    with the galactic binary confusion noise
    
`SnC::Function`: galactic binary confusion noise

`Ωn::Function`: effective energy density of noise
"""
struct LISA <: Detector
    TObs::Float # observation time of LISA (yr)
    LArm::Float # Arm length of for every arm of LISA (meter)
    NC::Int # Number of data channels
    
    approRQ # whether to approximate the transfer function or not
    
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
    
    Pn_WC::Function # Power Spectral Density with confusion (WC) noise
    Sn_WC::Function # Sensitivity curve of LISA 
    # with the galactic binary confusion noise
    
    SnC::Function # galactic binary confusion noise
    
    Ωn::Function # fractional energy density of noise
    Ωeff::Function # fractional effective energy density of noise
    
    ρThSGWB::Float # threshold SNR for deteting a gravitational background
    
    function LISA(;
            TObs::Float=4., 
            LArm::Float=2.5e9, 
            NC::Int=2, 
            approRQ=false, 
            fMin=1e-5, 
            fMax=1e0,
            fRef::Float=1e-3, 
            fPlotRange::Tuple{Float,Float}=(1e-5, 1e0), 
            ΩPlotRange::Tuple{Float,Float}=(1e-14, 1e-6), 
            hPlotRange::Tuple{Float,Float}=(1e-21, 1e-16),
            ρThSGWB=10.0)
        
        fStar = c0/(2*π*LArm)
        R = RLISA(fStar, NC, approRQ=approRQ)
        
        # single-link optical metrology noise (1/Hz), Eq.(10)
        P_OMS(f::Float) = (1.5e-11)^2 * (1 + (2e-3/f)^4)
        
        # single test mass acceleration noise, Eq.(11) 
        P_acc(f::Float) = (3e-15)^2 * (1 + (0.4e-3/f)^2) * (1 + (f/(8e-3))^4) 
        
        # total noise in Michelson-style LISA data channel, Eq.(12) 
        Pn(f) = P_OMS(f)/LArm^2 + 2(1 + cos(f/fStar)^2)*P_acc(f)/(2pi*f)^4/LArm^2
        
        Sn(f) = fMin<f<fMax ? Pn(f)/R(f) : Inf
        
        SnC(f) = getSnC(f, TObs, NC)
        Sn_WC(f) = Sn(f) + SnC(f)
        Pn_WC(f) = Sn_WC(f)*R(f)         
        
        Ωn(f) = (2*π^2/3/H0^2) * f^3 * Sn_WC(f)
        Ωeff(f) = (2/NC/(NC-1)) * Ωn(f)
        
        new(TObs, LArm, NC, approRQ, fStar, R, 
            fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange, 
            Pn, Sn, Pn_WC, Sn_WC, SnC, Ωn, Ωeff, ρThSGWB)
    end       
end


# TODO: add more explanation
function SNR(det::LISA, Ωgw::Function)
    
    fMin, fMax, fRef = det.fMin, det.fMax, det.fRef 
    
    T = det.TObs * YEAR # (s)
    
    integral(f) = (Ωgw(f)/det.Ωeff(f))^2
    
    logfMin = log10(det.fMin)
    logfMax = log10(det.fMax)
    f1, f2, f3, f4, f5, f6 = 10 .^ collect(range(logfMin, logfMax, length=6))
    
    int0 = quadgk(integral, f1, f2, f3, f4, f5, f6, rtol=1e-7)[1]
        
    snr = sqrt(T*int0)
end