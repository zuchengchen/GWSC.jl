function pta_noisepower(fmin, fmax, cadence, sigma)
%
% calculate timing noise psd (1-sided) for pulsar timing
%
% Pn(f) = 2 \Delta t sigma (units %
% then write results to file pta_noisepower.dat 
%
% fmin - min frequency (Hz)
% fmax - max frequency (Hz)
% cadence - number of observations (yr^-1)
% sigma - rms timing noise (ns)
%
% typically, fmin = 1e-9, fmax = 1e-7, cadence = 20, sigma = 100
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
yr = 3.15569e7; % s
ns = 1e-9; % s

% convert sigma to seconds
sigma = sigma*ns;

% discrete frequencies
N = 5000;
f = linspace(fmin, fmax, N);

% calculate Delta t (=1/cadence) in seconds
Delta_t = yr/cadence;

% calculate 1-sided noise psd
P = 2 * Delta_t * (sigma^2) * ones(1,N);

% write data to file
fid = fopen('pta_noisepower.dat','w');
for ii=1:N
  fprintf(fid, '%g\t%g\n', f(ii), P(ii));
end

fclose(fid);

return

