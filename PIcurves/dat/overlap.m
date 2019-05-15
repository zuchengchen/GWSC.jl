function orf = overlap(f, det1, det2, beta, method)
%
% calculate overlap reduction function
%
% see overlapScript.m for examples
% 
% you can use load('whatever filename was', '-mat') to 
% subsequently reload the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;   % speed of light (m/s)

% get detector geometry information
[r1, u1, v1, T1] = getdetectorNew(det1);
[r2, u2, v2, T2] = getdetectorNew(det2);
deltaX = r1-r2;

% get angles from healpix pixelization of the sky
%[theta,phi]=healpix2ang('pixelization_3072.dat');
[theta,phi]=healpix2ang('pixelization_12288.dat');
%[theta,phi]=healpix2ang('pixelization_49152.dat');
Npix = length(theta);
dArea = 4*pi/Npix;

% construct overlap reduction function
orf = zeros(1,length(f));

for ii = 1:1:Npix
  fprintf('working on %d of %d\n', ii, Npix);

  [F1p, F1c] = FpFc(f, theta(ii), phi(ii), 0, u1, v1, T1, method, 'id');
  [F2p, F2c] = FpFc(f, theta(ii), phi(ii), 0, u2, v2, T2, method, 'id');

  H = (1/(sin(beta)^2)) * (5/(8*pi)) * (F1p.*conj(F2p) + F1c.*conj(F2c));

  R = Ry(theta(ii)) * Rz(phi(ii));
  RdeltaX = R * deltaX;
  nDotDeltaX = RdeltaX(3);

  phaseFac = 2*pi*f*nDotDeltaX/c;

  orf   = orf + H.*exp(i*phaseFac)*dArea; 

end

% plot overlap reduction function
figure(1)
plot(f, real(orf), 'b');
xlabel('f (Hz)');
ylabel('\gamma(f)');
grid on
filename = ['overlap_' det1 '_' det2 '_' method '_linear.eps'];
%print('-depsc', filename);

figure(2)
semilogx(f, real(orf), 'b');
xlabel('f (Hz)');
ylabel('\gamma(f)');
grid on
filename = ['overlap_' det1 '_' det2 '_' method '_log.eps'];
%print('-depsc', filename);

% save data in a .mat file
filename = ['overlap_' det1 '_' det2 '_' method '.mat'];
save(filename, '-mat');

return
