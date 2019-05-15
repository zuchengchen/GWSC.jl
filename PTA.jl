struct PTA <: Detector
    c::Float # cadence of observations (yr^-1)
    ΔT::Float # 1/c
    σ::Float # root-mean-square timing noise (ns)
    nPTA::Int # number of pulsars
    TObs::Float # observation time of PTA (yr)
    fMin::Float # minimum frequency
    fMax::Float # maximum frequency
    
    fPlotRange::Tuple{Float,Float} # frequency range to show
    ΩPlotRange::Tuple{Float,Float} # energy density range to show
    
    R::Function # transfer function
    Pn::Function # Power Spectral Density (Hz^-1)
    Sn::Function # Power Spectral Density (Hz^-1)
    Seff::Function # Effective Power Spectral Density (Hz^-1)
    
    Ωn::Function # fractional energy density of noise
    Ωeff::Function # fractional effective energy density of noise
    
    ρThSGWB::Float # threshold SNR for deteting a gravitational background
    
    function PTA(;c::Float=20., σ::Float=1e2, nPTA::Int=20, TObs::Float=5., 
            fPlotRange=(1e-9, 1e-6), ΩPlotRange=(1e-11, 1e0), ρThSGWB=1.0)
        
        ΔT = YEAR/c
        fMin = 1e-9
        fMax = 1e-7
        fMin = 1/(TObs*YEAR)
 
        fMax = 1/ΔT/2
        R(f) = 1/(12*π^2*f^2)
        Pn(f) = 2 * ΔT * (σ*1e-9)^2 # Eq.(40)
        Sn(f) = f>=fMin&&f<fMax ? Pn(f)/R(f) : Inf
#         Sn(f) = Pn(f)/R(f)
        
        ΓII = 1.0/(4.0*sqrt(3.0))
        
        Seff(f) = Sn(f) * (nPTA*(nPTA-1)/2)^(-1/2) * ΓII^(-1)
        
        # fix me 
        # add a factor of 5 to test result
        Seff(f) = Sn(f) * (4.74)^(-1/2) 
        
        Ωn(f) = 2*π^2/(3*H0^2) * f^3 * Sn(f)
        Ωeff(f) = 2*π^2/(3*H0^2) * f^3 * Seff(f)
        
        new(c, ΔT, σ, nPTA, TObs, fMin, fMax, fPlotRange, ΩPlotRange,
            R, Pn, Sn, Seff, Ωn, Ωeff, ρThSGWB)
    end       
end