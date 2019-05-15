function pta_transfer(fmin, fmax)
%
% calculate single pulsar timing transfer function, 
%
% R(f) = GammaII(f) = 1/(2 pi f)^2 1/3
%
% the write results to file pta_transfer.dat
%
% fmin - min frequency (Hz) 
% fmax - max frequency (Hz) 
%
% typically, fmin = 1e-9, fmax = 1e-7
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate single pulsar timing transfer function
N = 5000;
f = linspace(fmin, fmax, N);
R = 1./(12*pi^2*f.^2);

% write data to file
fid = fopen('pta_transfer.dat','w');
for ii=1:N
  fprintf(fid, '%g\t%g\n', f(ii), R(ii));
end

fclose(fid);

return
