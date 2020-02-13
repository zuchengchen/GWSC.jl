"""
    PTA <: Detector

Detector type for pulsar timing array (PTA).

# Fields

- `Δt::Float=14DAY`: (s); 1/Δt is the cadence of observations

- `σRMS::Float=1e2`: root-mean-square timing noise (ns)

- `NP::Int=20`: number of pulsars

- `TObs::Float=5.`: observation time of PTA (yr)

- `fMin::Float`: minimum frequency (Hz) 

- `fMax::Float`: maximum frequency (Hz)

- `fRef::Float=1e-8`: reference frequency (Hz)

- `fPlotRange::Tuple{Float,Float}=(1e-9, 1e-7)`: frequency range to show (Hz, Hz)

- `ΩPlotRange::Tuple{Float,Float}=(1e-11, 1e-7)`: energy density range to show

- `hPlotRange::Tuple{Float,Float}=(1e-17, 1e-12)`: strain range to show

- `R::Function`: transfer function

- `Pn::Function`: detector power spectral density (Hz^-1)

- `Sn::Function`: strain noise power spectral density (Hz^-1)

- `Seff::Function`: effective strain noise power spectral density (Hz^-1)

- `Ωn::Function`: fractional energy density of noise

- `Ωeff::Function`: effective fractional energy density of noise

- `ρThSGWB::Float=1.0`: threshold SNR for deteting a gravitational wave background

# Examples
```jldoctest
julia> pta = PTA(Δt=14*DAY, σRMS=1e2, NP=36, TObs=5.);
```

"""
struct PTA <: Detector
    
    Δt::Float # (s); 1/Δt is the cadence of observations    
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
    
    function PTA(;
            Δt::Float=14DAY, 
            σRMS::Float=1e2, 
            NP::Int=20, 
            TObs::Float=5., 
            fRef::Float=1e-8,  
            ΩPlotRange::Tuple{Float,Float}=(1e-11, 1e-5), 
            hPlotRange::Tuple{Float,Float}=(1e-17, 1e-12), 
            ρThSGWB::Float=1.0)
        
        # convert units
        σ = σRMS*1e-9 # (s)
        T = TObs*YEAR # (s)
        
        fMin = 1/T # (Hz)
        fMax = 0.5/Δt # Nyquist frequency (Hz)
        fPlotRange=(fMin*1e-1, fMax*1e1)
        
        R(f) = 1/(12π^2*f^2)
        Pn(f) = fMin ≤ f ≤ fMax ? 2Δt*σ^2 : Inf # Eq.(40)
        Sn(f) = fMin ≤ f ≤ fMax ? 24π^2*f^2*Δt*σ^2 : Inf
        
        ζ2 = 1/48. # square of the average of Hellings and Downs factor
        
        Seff(f) = Sn(f) * (NP*(NP-1)/2 * ζ2)^(-1/2)
        
        Ωn(f) = fMin ≤ f ≤ fMax ? 16π^4*Δt*σ^2*f^5/H0^2 : Inf
        Ωeff(f) = 2π^2/(3H0^2) * f^3 * Seff(f)
        
        new(Δt, σRMS, NP, TObs, fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange,
            R, Pn, Sn, Seff, Ωn, Ωeff, ρThSGWB)
    end       
end


# # TODO: add more explanation
# function SNR(pta::PTA, β, Ωβ)
#     β, Ωβ = big(β), big(Ωβ)
    
#     fMin, fMax, fRef = pta.fMin, pta.fMax, big(pta.fRef) 
    
#     N = pta.NP
    
#     T = pta.TObs * YEAR # (s)
#     ζ2 = 1/48. # square of the average of Hellings and Downs factor
    
#     # B*f^b
#     integral(f) = 1.0/(1+pta.Ωn(f)/Ωgw(Ωβ, β, fRef)(f))^2
#     println(integral(pta.fRef))
    
#     int0 = quadgk(integral, fMin, fMax, rtol=1e-4)[1]
    
#     snr2 = 2T*N*(N-1)/2*ζ2*int0
#     snr = sqrt(snr2)
# end


# TODO: add more explanation
function SNR(det::PTA, Ωgw::Function)
    
    fMin, fMax, fRef = det.fMin, det.fMax, det.fRef 
    
    N = det.NP
    
    T = det.TObs * YEAR # (s)
    ζ2 = 1/48. # square of the average of Hellings and Downs factor
    
    integral(f) = 1.0/(1+det.Ωn(f)/Ωgw(f))^2.0
    
    logfMin = log10(det.fMin)
    logfMax = log10(det.fMax)
    f1, f2, f3, f4, f5, f6 = 10 .^ range(logfMin, logfMax, length=6)
    
    int0 = quadgk(integral, f1, f2, f3, f4, f5, f6, rtol=1e-7)[1]
    
    snr2 = 2T*N*(N-1)/2*ζ2*int0
    snr = sqrt(snr2)
end