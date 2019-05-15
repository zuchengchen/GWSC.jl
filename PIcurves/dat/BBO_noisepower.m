function BBO_noisepower(fmin, fmax)
%
% calculate noise psd (1-sided) for standard BBO 
% using data from crowder-cornish 2005, "Beyond
% LISA: Exploring Future GW Missions
%
% fmin - min frequency (Hz)
% fmax - max frequency (Hz)
%
% typically, fmin = 1e-3, fmax = 1e2 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;   % speed of light (m/s)
 
% discrete frequencies
N = 16*4000;
f = linspace(fmin, fmax, N);

% parameter from table II of crowder-cornish
Spos = 2e-34; % (\tilde delta x)^2
Saccel = 9e-34; % (\tilde delta a)^2 
L = 5e7; % arm length (m)
fstar = c/(2*pi*L);

% calculate 1-sided noise psd (see cornish-larson, eqs. 49,50)
% (factor of 4 times Spos is from cutler-harms)
P = 4*Spos/L^2 + Saccel*(2./(L*(2*pi*f).^2)).^2; 

% write data to file
fid = fopen('BBO_noisepower.dat','w');
for ii=1:N
  fprintf(fid, '%g\t%g\n', f(ii), P(ii));
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare against crowder-cornish

% load normalized transfer function
temp = load('BBO_transfer.dat');
f_R = temp(:,1);
R = temp(:,2); 
R = (sin(pi/3)^2/5) * R; % R = Gamma 
R = interp1(f_R, R, f);
S = P./R;

% make plots
figure(1)
loglog(f, R)
grid on
xlabel('f (Hz)')
ylabel('R(f)')

figure(2)
loglog(f, P, 'b');
hold on;
loglog(f, S, 'g');
hold on;
grid on;
xlim([1e-3 1e2])
xlabel('f (Hz)');
ylabel('P(f), S(f) Hz^{-1}');

figure(3)
loglog(f, (f.*S).^(1/2), 'b');
grid on;
xlim([1e-3 1e2])
xlabel('f (Hz)');
ylabel('h_c(f)');
print('-depsc2', 'BBO_comparison.eps');

return
