function BBO_transfer()
%
% calculate normalized BBO transfer function
% (transfer function for a mini-LISA with L=5e7 m arms)
%
% R(f) = gammaII(f) 
%
% then write result to file BBO_transfer.dat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data 
% (first use overlapExamples.m, example=4)
load('overlap_default_default_exact_broadband.mat');

% rescale dimensionless frequencies by c/L;
c = 299792458;   % speed of light (m/s)
L = 5e7; % LISA armlength (m)
f = f*c/L;

R = real(orf);

% write data to file
fid = fopen('BBO_transfer.dat','w');
for ii=1:length(f)
  fprintf(fid, '%g\t%g\n', f(ii), R(ii));
end

fclose(fid);

return
