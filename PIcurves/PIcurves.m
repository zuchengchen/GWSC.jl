function [f, fplot, df, Omega, S, h, Omega_eff, S_eff, h_eff, Omega_n, S_n, ...
  h_n, web] = PIcurves(noisefile, transferfile, orffile, T, fmin, fmax, ...
  fplot_min, fplot_max, params)
%
% by Eric Thrane and Joe Romano
% contact: eric.thrane@ligo.org, joseph.romano@ligo.org
% Calculates and Power-law Integrated (PI) curves
%
% INPUTS:
%
%   noisefile    = path to file containing ifo or pulsar timing noise psd
%   transferfile = path to file containing ifo or pulsar timing transfer
%                  function
%   orffile      = path to file containing the overlap reduction function 
%                  for a pair of ifos or a list of (ra,dec) for a set
%                  of pulsars comprising a pulsar timing array
%   (each file has two columns: the first column is the frequency in Hz;
%   the second colum are the values of the noise psd, transfer function, etc.)
%
%   NOTE: for ifos, the transfer function and overlap reduction are assumed
%   to be normalized to unity in the limit f->0 for colocated and coaligned
%   detectors.  For a pulsar timing array, the transfer function is
%   R(f)=1/(2pi f^2)*(1/3) and the Hellings-Downs factors can be computed from
%   the pulsars file.
%
%   T                     = observation time in years
%   [fmin fmax]           = start and stop frequencies for detector noise
%   [fplot_min fplot_max] = start and stop frequencies for plotting PI curves
%   params                = struct describing details of plot.  Subfields are:
%                           ifo_pair: true or false
%                           pta_network: true or false
%                           (Must be one or the other!)
%                           rho: expected squared SNR
%                           beta: opening angle for ifo; required if ifo_pair=1
% OUTPUTS:
%   f             = array of frequencies (Hz)
%   fplot         = array of plot frequencies (Hz)
%   df            = frequency spacing (Hz)
%   Omega, S, h   = power-law integrated sensitivity curves
%   Omega_eff, S_eff, h_eff = multi-detector sensitivity curves
%   Omega_n, S_n, h_n = single detector sensitivity curves: 
%   web           = a struct with subfields: Omega, S, and h.  Each one
%                   contains the full set of power-law curves used in the
%                   construction of the power-law integrated curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some constants
HubbleConstant = 3.2408e-18;
H100 = 0.68 * HubbleConstant; % use h0 from the Planck result.
yr = 3.15569e7; % s

% convert T to seconds
T = T*yr;

% load and parse files
temp = load(noisefile);
f_noise = temp(:,1);
P_n = temp(:,2);
  
temp = load(transferfile);
f_transfer = temp(:,1);
transfer = temp(:,2);

%--------CHECK PARAM VALUES----------------------------------------------------
% define parameter defaults beginning with ifo_pair
if ~isfield(params, 'ifo_pair')
  params.ifo_pair = false;
end
% pta_network
if ~isfield(params, 'pta_network')
  params.pta_network = false;
end
% make sure pta_network of ifo_pair is defined
if params.ifo_pair==false & params.pta_network==false
  error('Set one of the following fields to true: params.ifo_pair or params.pta_network.');
end
% ...but not both
if params.ifo_pair==true & params.pta_network==true
  error('Set ONLY one of the following fields to true: params.ifo_pair or params.pta_network.');
end
% if ifo_pair==true, make sure that beta is defined
if params.ifo_pair
  if ~isfield(params, 'beta')
    error('ifo_pair requires params.beta.');
  end
end
% make sure rho is defined
if ~isfield(params, 'rho')
  fprintf('params.rho is not defined; setting rho=1.\n');
  params.rho = 1;
end
%------------------------------------------------------------------------------

% for an ifo pair, the orf is assumed to be normalized to unity as f->0
% for a pta network, the orf file contains a list of pulsars so do nothing for that
if params.ifo_pair
  temp = load(orffile);
  f_gamma = temp(:,1);
  gamma = temp(:,2);
  %  apply normalization factor to the orf for an ifo pair
  if params.ifo_pair
    Gamma = (sin(params.beta)^2/5)*gamma;
    transfer = (sin(params.beta)^2/5)*transfer;
    fprintf('normalizing Gamma_IJ(f) and R(f) for ifos\n');
  end
end

% chose an appropriate bin for interpolation
% use fmax as a guide
if fmax>10
  % use typical LIGO bin width
  df = 0.25;
elseif fmin < 1/(0.5*yr)
  % pulsar timing regime
  df = 1/T;
else
  % scale based on fmax
  df = fmax/400;
end

% make sure that df>=1/T
if df<1/T
  fprintf('Warning: setting df=1/T.\n');
  df = 1/T;
end
fprintf('Assuming a bin width of %1.1e Hz\n', df);
f = [fmin:df:fmax]';

% interpolate relevant functions
P_n = interp1(f_noise, P_n, f);
transfer = interp1(f_transfer, transfer, f);
if params.ifo_pair
  Gamma = interp1(f_gamma, Gamma, f);
end

% calculate single-detector strain power S_n, energy density Omega_n,
% and strain amplitude h_n
S_n = P_n./transfer;
Omega_n = (2 * pi^2)/(3 * H100^2) * f.^3 .* S_n;
h_n = (f .* S_n).^(1/2);

% calculate multi-detector strain power Seff
if params.ifo_pair
  % S_eff(f) = sqrt(P_I*P_J)/abs(Gamma_IJ), but assume P_I = P
  S_eff = P_n./abs(Gamma);

elseif params.pta_network
  % S_eff = [sum_I sum_J>I Gamma_{IJ}^2 /P_I P_J]^{-1/2}, but assume P_I = P
  % so that S_eff = S_n / sum_I sum_J>I zeta_IJ^2 
  [Npairs_eff, zetaIJ, IJ, thetaIJ] = sumHDfactors2(orffile);
  fprintf('Effective number of pulsar pairs = %f\n',  Npairs_eff)
  S_eff = S_n/sqrt(Npairs_eff);
  
else
  error('unknown type of detector!')

end

% convert Seff to effective energy density and effective strain amplitude 
Omega_eff = (2 * pi^2)/(3 * H100^2) * f.^3 .* S_eff;
h_eff = (f .* S_eff).^(1/2);

% define power law indices
beta = -8:1:8;

% calculate power-law Omega_beta that yields desired squared snr value
rho = params.rho;

for ii=1:length(beta)
  b = beta(ii);

  integral = 0;
  for jj=1:length(f)
    integral = integral + df * f(jj)^(2*b) / Omega_eff(jj)^2;
  end 

  OmegaByFrefBeta(ii) = (rho/sqrt(2*T)) / sqrt(integral);

end

% calculate power law integrated curves for plot range frequencies
fplot = [fplot_min:df:fplot_max];

for ii=1:length(beta)
  b = beta(ii);

  % calculate power law omega
  Omega(ii,:) = OmegaByFrefBeta(ii) * fplot.^b;

  % convert power-law Omega to S and h 
  S(ii,:) = (3 * H100^2)/(2 * pi^2) * Omega(ii,:)./(fplot.^3);
  h(ii,:) = (fplot .* S(ii,:)).^(1/2);
end

% save struct containing "webs" of power-law curves
web.Omega = Omega;
web.S = S;
web.h = h;

% power-law integrated curves are the locus of these power-law curves
Omega = max(Omega, [], 1);
S = max(S, [], 1);
h = max(h, [], 1);

return

