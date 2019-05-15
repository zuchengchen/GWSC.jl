function overlapScript(example)
%
% script for running overlap.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
c = 299792458;   % speed of light (m/s)

switch example
  case 1
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'H1';
    det2 = 'L1';
    beta = pi/2;
    method = 'lw';

  case 2
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'H1';
    det2 = 'L1';
    beta = pi/2;
    method = 'exact';

  case 3;
    N = 8192;
    det1 = 'default';
    det2 = 'default';
    beta = pi/2;
    T = 1;
    fsr = 1/(2*T);
    flow = 1e-3*fsr;
    fhigh = 10*fsr;
    f = linspace(flow, fhigh, N);
    method = 'exact';

  case 4;
    N = 32768;
    det1 = 'default';
    det2 = 'default';
    beta = pi/2;
    T = 1;
    fsr = 1/(2*T);
    flow = 1e-5*fsr;
    fhigh = 1e5*fsr;
    f = logspace(log10(flow), log10(fhigh), N);
    method = 'exact';

  case 5;
    N = 8192;
    det1 = 'BBO1';
    det2 = 'BBO2';
    beta = pi/3;
    T = 5e7/c;
    fsr = 1/(2*T);
    flow = 1e-3*fsr;
    fhigh = 1e1*fsr;
    f = logspace(log10(flow), log10(fhigh), N);
    method = 'exact';

  case 6
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'H1';
    det2 = 'V1';
    beta = pi/2;
    method = 'lw';

  case 7
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'L1';
    det2 = 'V1';
    beta = pi/2;
    method = 'lw';

  case 8
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'H1';
    det2 = 'T1';
    beta = pi/2;
    method = 'lw';

  case 9
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'L1';
    det2 = 'T1';
    beta = pi/2;
    method = 'lw';

  case 10
    N = 4096;
    f = linspace(1, 2048, N);
    det1 = 'T1';
    det2 = 'V1';
    beta = pi/2;
    method = 'lw';

end

overlap(f, det1, det2, beta, method);

return

