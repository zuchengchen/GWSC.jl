"""
Detector type for pulsar timing array.
"""
struct PTA <: Detector
    c::Float # cadence of observations (yr^-1)
    Δt::Float # 1/c (s)
    σRMS::Float # root-mean-square timing noise (ns)
    NP::Int # number of pulsars
    TObs::Float # observation time of PTA (yr)
    
    fMin::Float # minimum frequency (Hz) 
    fMax::Float # maximum frequency (Hz)
    fRef::Float # reference frequency (Hz)
    
    fPlotRange::Tuple{Float,Float} # frequency range to show (Hz, Hz)
    ΩPlotRange::Tuple{Float,Float} # energy density range to show
    hPlotRange::Tuple{Float,Float} # strain range to show
    
    R::Function # transfer function
    Pn::Function # detector power spectral density (Hz^-1)
    Sn::Function # strain noise power spectral density (Hz^-1)
    Seff::Function # effective strain noise power spectral density (Hz^-1)
    
    Ωn::Function # fractional energy density of noise
    Ωeff::Function # effective fractional energy density of noise
    
    ρThSGWB::Float # threshold SNR for deteting a gravitational wave background
    
    function PTA(;c::Float=20., σRMS::Float=1e2, NP::Int=20, TObs::Float=5., 
            fRef::Float=1e-8, fPlotRange=(1e-9, 1e-7), 
            ΩPlotRange=(1e-11, 1e-7), hPlotRange=(1e-17, 1e-12), ρThSGWB=1.0)
        
        # convert units
        Δt = YEAR/c # (s)
        σ = σRMS*1e-9 # (s)
        T = TObs*YEAR # (s)
        
        fMin = 1/T # (Hz)
        fMax = 0.5/Δt # Nyquist frequency (Hz)
        
        R(f) = 1/(12π^2*f^2)
        Pn(f) = f>=fMin&&f<fMax ? 2Δt*σ^2 : Inf # Eq.(40)
        Sn(f) = f>=fMin&&f<fMax ? 24π^2*f^2*Δt*σ^2 : Inf
        
        ζ2 = 1/48. # square of the average of Hellings and Downs factor
        
        Seff(f) = Sn(f) * (NP*(NP-1)/2 * ζ2)^(-1/2)
        
        Ωn(f) = f>=fMin&&f<fMax ? 16π^4*Δt*σ^2*f^5/H0^2 : Inf
        Ωeff(f) = 2π^2/(3H0^2) * f^3 * Seff(f)
        
        new(c, Δt, σRMS, NP, TObs, fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange,
            R, Pn, Sn, Seff, Ωn, Ωeff, ρThSGWB)
    end       
end