% commands for producing the BBO overlap plot

clear all
close all
data = load('BBO1BBO2_orf.dat');
f = data(:,1);
orf = data(:,2);

semilogx(f, real(orf), 'b')
grid on
xlim([1e-2 10])
xlabel('f (Hz)', 'fontsize', 24);
ylabel('\gamma(f)', 'fontsize', 24);

print -depsc2 BBOoverlap.eps
 
