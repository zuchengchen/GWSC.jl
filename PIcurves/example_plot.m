function example_plot
% function example_plot

% define parameters
% noise curve
noisefile = 'dat/ZERO_DET_high_P_psd.txt';
% transfer function
transferfile = 'dat/LIGO_transfer.dat';
% overlap reduction function
orffile = 'dat/H1L1_orf.dat';
% observation time in years
T = 1;
% analysis band in Hz
fmin = 10;
fmax = 2000;
% plotting band in Hz
fplot_min = fmin;
fplot_max = fmax;
% specify whether this is an ifo_pair or a pta_network
params.ifo_pair = true;
% since it is an ifo_pair, define the ifo opening angle
params.beta = pi/2;
% define the SNR squared
params.rho = 1;

% calculate power-law integrated curve
[f, fplot, df, Omega, S, h, Omega_eff, S_eff, h_eff, Omega_n, S_n, h_n, ...
  web] = PIcurves(noisefile, transferfile, orffile, T, fmin, fmax, ...
  fplot_min, fplot_max, params);

% plot curve
figure;
loglog(f, web.Omega, 'k');
hold on;
loglog(f, Omega, 'linewidth', 2);
axis([fmin fmax 1e-10 1e-4]);
grid on;
xlabel('f (Hz)');
ylabel('\Omega(f)');
print('-dpng', 'example_locus');

return
