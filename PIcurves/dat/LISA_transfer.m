function LISA_transfer()
%
% calculate normalized LISA transfer function
%
% R(f) = gammaII(f) 
%
% then write result to file LISA_transfer.dat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data 
load('overlap_default_default_exact.mat');

% rescale dimensionless frequencies by c/L;
c = 299792458;   % speed of light (m/s)
L = 5e9; % LISA armlength (m)
f = f*c/L;

R = real(orf);

% write data to file
fid = fopen('LISA_transfer.dat','w');
for ii=1:length(f)
  fprintf(fid, '%g\t%g\n', f(ii), R(ii));
end

fclose(fid);

return
