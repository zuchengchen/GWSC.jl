""" 
plotCharacteristicStrain(lisa::LISA; figure_file=undef, 
        plotSn_WC=true, plotSn=false, plotPn=false)

Plot the sensitivity curve of characteristic strain.    
If `figure_file` is provided, the figure will be saved.

`lisa::LISA`: lisa is a LISA type;

`figure_file=undef`: file name of the figue to be saved to;
        
`plotSn_WC=true`: plot `√f*Sn_WC` if it is true;

`plotSn=false`: plot `√f*Sn` if it is true;

`plotPn=false`: plot `√f*Pn` if it is true;

`plotPn_WC=false`: plot `√f*Pn_WC` if it is true;
"""
function plotCharacteristicStrain(lisa::LISA; file=nothing, 
        plotSn_WCQ=true, plotSnQ=false, plotPnQ=false, plotPn_WCQ=false)
    
    f = 10 .^ range(-5, stop=0, length=10^3)
    
    grid("on", which="both", linestyle="--")

    xlabel("Frequency (Hz)") 
    ylabel("Characteristic Strain") 
    
    xlim(1.0e-5, 1.0e0)
    ylim(3.0e-22, 1.0e-15)
    
    if plotSn_WCQ != false
        Sn_WC = lisa.Sn_WC.(f)
        fSn_WC = sqrt.(f .* Sn_WC)
        loglog(f, fSn_WC, label=L"(f Sₙ)$^{1/2}$")  
    end
    
    if plotSnQ == true
        Sn = lisa.Sn.(f)
        fSn = sqrt.(f .* Sn)
        loglog(f, fSn, label=L"(f Sₙ)$^{1/2}$")
    end
    
    if plotPnQ == true
        Pn = lisa.Pn.(f)
        fPn = sqrt.(f .* Pn)
        loglog(f, fPn, label=L"(f Pₙ)$^{1/2}$")
    end
    
    if plotPn_WCQ == true
        Pn_WC = lisa.Pn_WC.(f)
        fPn_WC = sqrt.(Pn_WC)
        loglog(f, fPn_WC, label=L"(f Pₙ)$^{1/2}$")
    end    
    
    legend()
    show()
    
    # save figure
    tight_layout() # save room for the labels
    if (file != nothing)
        savefig(file)
    end
end

""" 
plotSpectralDensity(lisa::LISA; figure_file=undef, 
        plotSn_WC=true, plotSn=false, plotPn=false, plotPn_WC=false)

Plot the sensitivity curve of characteristic strain.    
If `figure_file` is provided, the figure will be saved.

`lisa::LISA`: lisa is a LISA type;

`figure_file=undef`: file name of the figue to be saved to;
        
`plotSn_WC=true`: plot `√Sn_WC` if it is true;

`plotSn=false`: plot `√Sn` if it is true;

`plotPn=false`: plot `√Pn` if it is true;

`plotPn_WC=false`: plot `√Pn_WC` if it is true.
"""
function plotSpectralDensity(lisa::LISA; file=nothing, 
        plotSn_WCQ=true, plotSnQ=false, plotPnQ=false, plotPn_WCQ=false, plotSnCQ=false)
    
    f = 10 .^ range(-5, stop=0, length=10^3)
    
    grid("on", which="both", linestyle="--")

    xlabel("Frequency (Hz)")
    ylabel(L"Spectral Density (Hz$^{-1/2}$)")
    
    xlim(1.0e-5, 1.0e0)
    ylim(3.0e-21, 1.0e-14)
        
    if plotSn_WCQ != false
        Sn_WC = lisa.Sn_WC.(f)
        fSn_WC = sqrt.(Sn_WC)
        loglog(f, fSn_WC, label=L"Sₙ$^{1/2}$") 
    end
    
    if plotSnQ == true
        Sn = lisa.Sn.(f)
        fSn = sqrt.(Sn)
        loglog(f, fSn, label=L"Sₙ$^{1/2}$")
    end
    
    if plotPnQ == true
        Pn = lisa.Pn.(f)
        fPn = sqrt.(Pn)
        loglog(f, fPn, label=L"Pₙ$^{1/2}$")
    end
    
    if plotPn_WCQ == true
        Pn_WC = lisa.Pn_WC.(f)
        fPn_WC = sqrt.(Pn_WC)
        loglog(f, fPn_WC, label=L"Pₙ$^{1/2}$")
    end
       
    if plotSnCQ == true
        SnC = lisa.SnC.(f)
        fSnC = sqrt.(SnC)
        loglog(f, fSnC, label=L"S$_\mathrm{c}^{1/2}$")
    end
    
    legend()
    show()
    
    # save figure
    tight_layout() # save room for the labels
    if (file != nothing)
        savefig(file)
    end
end

