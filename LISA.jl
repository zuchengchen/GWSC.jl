# LISA modulation frequency
fm = 3.168753575e-8  

"""
fstar(LArm::Float=2.5e9)

Return the transfer frequency (Hz).

`LArm`: Arm length of for every arm of LISA (meter).
"""
fstar(LArm::Float=2.5e9) = c0/(2*π*LArm)


""" 
loadTransfer(transfer_file::String, fStar::Float=19.09e-3, 
        NC::Int=2)

Load the data file containing the numerically calculate 
transfer function (sky and polarization averaged)
and return an interpolated/approximated transfer function.

`transfer_file`: name of the data file containing transfer function;

`fStar`: transfer frequency;

`NC`: number of data channels for LISA.
"""
function loadTransfer(transfer_file::String, fStar::Float=19.09e-3, 
        NC::Int=2)  
    RData = try    # try to read in the data file
         readdlm(transfer_file) # read in the data
    catch
        # If file isn't successfully read in, use approximate transfer function
        print("Warning: Could not find transfer function file ")
        printstyled("$(transfer_file)!\n", color=:red)
        print("\t Approximation will be used...")
        
        # return an approximated version of transfer function
        return f::Float -> NC * 3.0/20.0/(1.0 + 0.6*(f/fStar)^2)
    end
    
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

# single-link optical metrology noise (1/Hz), Eq.(10)
P_OMS(f::Float) = (1.5e-11)^2 * (1. + (2.0e-3/f)^4) 

# single test mass acceleration noise, Eq.(11)
P_acc(f::Float) = (3.0e-15)^2 * (1. + (0.4e-3/f)^2) * (1. + (f/(8.0e-3))^4) 

# total noise in Michelson-style LISA data channel, Eq.(12)
getPn(f::Float, LArm::Float=2.5e9) = P_OMS(f)/LArm^2 + 2.0*(1. + cos(f/fstar(LArm))^2) * P_acc(f)/(2.0*pi*f)^4/LArm^2 

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
    A*f^(-7.0/3.) * exp(-f^α + β*f*sin(κ*f)) * (1.0 + tanh(γ*(fk - f)))
end


""" 
getSn(f::Float, transfer_file::String="RLISA.txt", 
        LArm::Float=2.5e9, NC::Int=2)  

Sensitivity curve of LISA without the galactic binary confusion noise.

`f`: frequency (Hz);

`transfer_file`: name of the data file containing transfer function;

`LArm`: Arm length of for every arm of LISA (meter);

`NC`: number of data channels for LISA.
"""
function getSn(f::Float, transfer_file::String="RLISA.txt", 
        LArm::Float=2.5e9, NC::Int=2)  
    
    R = loadTransfer(transfer_file, fstar(LArm), NC)
    getPn(f, LArm)/R(f) 
end


""" 
getSn_WC(f::Float, transfer_file::String="RLISA.txt", 
    LArm::Float=2.5e9, TObs::Float=4.0, NC::Int=2)

Calculate the sensitivity curve of LISA accounting for
the galactic binary confusion noise.

`f`: frequency (Hz);

`transfer_file`: name of the data file containing transfer function;

`LArm`: Arm length of for every arm of LISA (meter);

`TObs`: observation time of LISA (yr);

`NC`: number of data channels for LISA.
"""
getSn_WC(f::Float, transfer_file::String="RLISA.txt", 
    LArm::Float=2.5e9, TObs::Float=4.0, NC::Int=2) =
getSn(f, transfer_file, LArm, NC) + getSnC(f, TObs, NC)


""" 
getPn_WC(f::Float, transfer_file::String="RLISA.txt", 
    LArm::Float=2.5e9, TObs::Float=4.0, NC::Int=2)

Calculate Power Spectral Density with confusion (WC) noise estimate.

`f`: frequency (Hz);

`transfer_file`: name of the data file containing transfer function;

`LArm`: Arm length of for every arm of LISA (meter);

`TObs`: observation time of LISA (yr);

`NC`: number of data channels for LISA.
"""
function getPn_WC(f::Float, transfer_file::String="RLISA.txt", 
    LArm::Float=2.5e9, TObs::Float=4.0, NC::Int=2)
    
    R = loadTransfer(transfer_file, fstar(LArm), NC)
    getPn(f, LArm) + getSnC(f, TObs, NC)*R(f)
end


"""
    LISA type

`TObs::Float`: observation time of LISA (yr);

`LArm`: Arm length of for every arm of LISA (meter);

`NC`: number of data channels for LISA;

`transfer_file`: name of the data file containing transfer function;

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
    
    # name of the data file containing transfer function
    transfer_file::String 
    
    fStar::Float # transfer frequency (Hz)
    R::Function # transfer function
    
    fMin::Float # minimum frequency
    fMax::Float # maximum frequency
    
    fPlotRange::Tuple{Float,Float} # frequency range to show
    ΩPlotRange::Tuple{Float,Float} # energy density range to show
    
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
    
    function LISA(;TObs::Float=4., LArm::Float=2.5e9, NC::Int=2, 
            transfer_file="RLISA.txt", fMin=1e-5, fMax=1e0,
            fPlotRange=(1e-5, 1e0), ΩPlotRange=(1e-14, 1e-6), ρThSGWB=10.0)
        
        fStar = fstar(LArm)
        R = loadTransfer(transfer_file, fStar, NC)
        Pn(f) = getPn(f, LArm)
        Sn(f) = getSn(f, transfer_file, LArm, NC)
        
        Pn_WC(f) = getPn_WC(f, transfer_file, LArm, TObs, NC)        
        Sn_WC(f) = getSn_WC(f, transfer_file, LArm, TObs, NC)
        
        SnC(f) = getSnC(f, TObs, NC) 
        
        Ωn(f) = (2*π^2/3/H0^2) * f^3 * Sn_WC(f)
        Ωeff(f) = (2/NC/(NC-1)) * Ωn(f)
        
        new(TObs, LArm, NC, transfer_file, fStar, R, 
            fMin, fMax, fPlotRange, ΩPlotRange, Pn, Sn, 
            Pn_WC, Sn_WC, SnC, Ωn, Ωeff, ρThSGWB)
    end       
end