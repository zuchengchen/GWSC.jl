function [theta,phi] = healpix2ang(filename)

data = load(filename);
theta = data(:,1);
phi = data(:,2);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of pixels
fprintf('Number of pixels = %d\n',length(theta));

% plot points
plot(phi, theta, '*');
xlim([0 2*pi])
ylim([0 pi])
xlabel('phi')
ylabel('theta')

