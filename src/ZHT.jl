struct ZHT <: Detector

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

    function ZHT(;
        name="ZHT",
        TObs=5.0,
        NDet=2,
        fRef=25.0,
        fPlotRange::Tuple{Float,Float}=(5e0, 1e4),
        ΩPlotRange::Tuple{Float,Float}=(1e-14, 1e-6),
        hPlotRange::Tuple{Float,Float}=(1e-21, 1e-16),
        ρThSGWB=1.0)

        # read in the data interpolate ASD 
        ASD_file = joinpath(sensitivity_path, ASDFile(name))
        #         println("using $ASD_file")
        # ASD_file = "data_zht_psd.csv"
        ASDData = readdlm(ASD_file, ',', Float64)
        ASDData = ASDData[sortperm(ASDData[:, 1]), :]

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
        # overlap_file = joinpath(sensitivity_path, "H1L1_orf.dat")
        # γData = readdlm(overlap_file)
        # fsγ = γData[:, 1]
        # γs = γData[:, 2]

        # fMinγ, fMaxγ = fsγ[1], fsγ[end]
        # fMin = max(fMinASD, fMinγ)
        # fMax = min(fMaxASD, fMaxγ)

        # γ = interLogxIn(fsγ, γs)

        # Γ(f) = 1 / 5 * γ(f)
        fMin, fMax = fMinASD, fMaxASD
        Γ(f) = 1 / 5

        Ωn(f) = (2 * π^2 / 3 / H0^2) * f^3 * Sn(f) # eq.(3)
        Ωeff(f) = sqrt((2 / NDet / (NDet - 1))) * Ωn(f) / Γ(f)

        new(name, TObs, NDet, R,
            fMin, fMax, fRef, fPlotRange, ΩPlotRange, hPlotRange,
            Pn, Sn, Ωn, Ωeff, ρThSGWB, ASD)
    end
end

function SNR(det::ZHT, Ωgw::Function)

    fMin, fMax, fRef = det.fMin, det.fMax, det.fRef

    T = det.TObs * YEAR # (s)

    integral(f) = (Ωgw(f) / det.Ωeff(f))^2

    logfMin = log10(det.fMin)
    logfMax = log10(det.fMax)
    f1, f2, f3, f4, f5, f6 = 10 .^ range(logfMin, logfMax, length=6)

    int0 = quadgk(integral, f1, f2, f3, f4, f5, f6, rtol=1e-6)[1]

    snr = sqrt(2T * int0)
end