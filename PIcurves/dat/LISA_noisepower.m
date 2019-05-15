function LISA_noisepower(fmin, fmax)
%
% calculate noise psd (1-sided) for standard LISA
% using data from cornish-larson 
%
% fmin - min frequency (Hz)
% fmax - max frequency (Hz)
%
% typically, fmin = 1e-4, fmax = 1e-1 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% discrete frequencies
N = 5000;
f = linspace(fmin, fmax, N);

% calculate 1-sided noise psd
%P = 9.24e-40*(0.001./f).^4 + 1.6e-41; % from cornish-larson 2001

% parameters from table II of crowder-cornish
Spos = 4e-22; % (\tilde delta x)^2
Saccel = 9e-30; % (\tilde delta a)^2 
L = 5e9; % arm length (m)

% calculate 1-sided noise psd (see cornish-larson, eqs. 49,50)
% (factor of 4 times Spos is from cutler-harms)
P = Spos/L^2 + Saccel*(2./(L*(2*pi*f).^2)).^2; 

% write data to file
fid = fopen('LISA_noisepower.dat','w');
for ii=1:N
  fprintf(fid, '%g\t%g\n', f(ii), P(ii));
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare against Shane Larson's strain power curve
% from http://www.srl.caltech.edu/~shane/sensitivity/MakeCurve.html
temp = load('lisa_psd_larson.dat');
fII = temp(:,1);
SII = temp(:,2);

% load normalized transfer function
temp = load('LISA_transfer.dat');
f_R = temp(:,1);
R = temp(:,2); 
R = (sin(pi/3)^2/5) * R; % R = Gamma 
R = interp1(f_R, R, f);
S = P./R;

% make a plot
figure;
loglog(f, P, 'b');
hold on;
loglog(f, S, 'g');
hold on;
loglog(fII, SII, 'r--');
grid on;
xlim([1e-4 1e-1])
xlabel('f (Hz)');
ylabel('P(f), S(f) Hz^{-1}');

print('-depsc2', 'LISA_comparison.eps');

return
