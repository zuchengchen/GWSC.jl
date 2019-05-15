function write_orf(filename, det1, det2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assumes that overlapExamples.m was used to generate .mat file
load(filename, '-mat');
gamma = real(orf);

% make plot
figure(1)
semilogx(f, gamma, 'b')
xlabel('f (Hz)','fontsize',24);
ylabel('\gamma(f)','fontsize',24)
set(gca,'XTick',logspace(-2,1,4))
grid on
titlestr = [det1 '-' det2];
title(titlestr)

% print to file
outputfilename = [det1 det2 '_orf.eps'];
print('-depsc2', outputfilename);

% write data to file
outputfilename = [det1 det2 '_orf.dat'];
fid = fopen(outputfilename, 'w');
for ii=1:length(f)
  fprintf(fid, '%g\t%g\n', f(ii), gamma(ii));
end

fclose(fid);

return

