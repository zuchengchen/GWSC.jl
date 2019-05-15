% commands for producing the gammaII plot
% (first use overlapExamples.m, example=3)

clear all
close all
load('overlap_default_default_exact.mat', '-mat')
fsr = f(end)/10;
%semilogx(f/fsr, real(orf), 'b')
loglog(f/fsr, real(orf), 'b')
grid on
xlim([1e-2 10])
xlabel('2fL/c', 'fontsize', 24);
ylabel('\gamma_I_I(f)', 'fontsize', 24);

print -depsc2 gammaII.eps
 
