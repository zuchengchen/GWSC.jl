% calculate and plot network overlap reduction function 
% for ground-based interferometers located at H1, L1, T1, V1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

N = 4;
Npairs = N*(N-1)/2;

% assumes that overlapExamples.m was used to generate .mat files
load('overlap_H1_L1_lw.mat', '-mat');
gamma(1,:) = real(orf);
load('overlap_H1_T1_lw.mat', '-mat');
gamma(2,:) = real(orf);
load('overlap_H1_V1_lw.mat', '-mat');
gamma(3,:) = real(orf);
load('overlap_L1_T1_lw.mat', '-mat');
gamma(4,:) = real(orf);
load('overlap_L1_V1_lw.mat', '-mat');
gamma(5,:) = real(orf);
load('overlap_T1_V1_lw.mat', '-mat');
gamma(6,:) = real(orf);

% make plot
figure(1)
semilogx(f, gamma(1,:), f, gamma(2,:), f, gamma(3,:), f, gamma(4,:), f, gamma(5,:), f, gamma(6,:));
xlabel('f (Hz)','fontsize',24);
ylabel('\gamma(f)','fontsize',24)
xlim([1 2e3])
set(gca,'XTick',logspace(0,3,4))
grid on
legend('H1L1', 'H1K1', 'H1V1', 'L1K1', 'L1V1', 'K1V1')

% print to file
print('-depsc2', 'all_orf.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate network overlap reduction function for S_eff(f)
% assuming P_nI(f) = P_n(f) for all detectors

temp = zeros(1,length(f));
for ii=1:Npairs
  temp = temp + gamma(ii,:).^2;  
end
gamma_network = temp.^(1/2);

% make plot
figure(2)
semilogx(f, gamma_network)
xlabel('f (Hz)','fontsize',24);
%ylabel('\gamma_{network}(f)','fontsize',24)
ylabel('R_{eff}(f)','fontsize',24)
xlim([1 2e3])
set(gca,'XTick',logspace(0,3,4))
grid on

% print to file
print('-depsc2', 'network_orf.eps');

% write to file
fid = fopen('network_orf.dat', 'w');
for ii=1:length(f)
  fprintf(fid, '%g\t%g\n', f(ii), gamma_network(ii));
end

fclose(fid);

return

