% addpath('/Volumes/Music/Dropbox/Documents/Matlab/jlab')
addpath('/Users/jearly/Dropbox/Documents/Matlab/jlab')
addpath('../GLOceanKit/Matlab/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Drifter Spectra
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = '/Volumes/Data/QGPlusSlab/TurbulenceExperimentNonStiff';
file = 'QGDampedSlabTrajectories.mat';
figure_out = sprintf('%s/DrifterSpectra.png',folder);

load(sprintf('%s/%s',folder,file));

numDrifters = size(xpos1,2);
days = t/86400;
dt = t(2)-t(1);

cv1 = (diff(xpos1,1,1) + sqrt(-1)*diff(ypos1,1,1))/dt;
cv2 = (diff(xpos2,1,1) + sqrt(-1)*diff(ypos2,1,1))/dt;

[psi,lambda]=sleptap(size(cv1,1),3);
[omega,spp,snn,spn]=mspec(dt,cv1,psi);

% convert from radians/second to cycles/day
f=omega*86400/(2*pi);

omega_floats = [ flip(-omega(2:end),1); omega];
f_floats = [ flip(-f(2:end),1); f];
S_floats = cat(1, flip(snn(2:end,:),1), spp);

S_mean = mean(S_floats,2);
S_std = std(S_floats,0,2);
inertialRange = find(f_floats < -0.63 & f_floats > -1.0);
inertial_speed = sqrt( (1/2/pi)*trapz(omega_floats(inertialRange),S_mean(inertialRange)));
inertial_std1 = sqrt( (1/2/pi)*trapz(omega_floats(inertialRange),S_mean(inertialRange) + S_std(inertialRange)));
inertial_std2 = sqrt( (1/2/pi)*trapz(omega_floats(inertialRange),S_mean(inertialRange) - S_std(inertialRange)));


% a == b is true for a complex time series
Z = cv1(:,1);
[F,SPP,SNN]=mspec(dt,Z,[]);
df = F(2)-F(1);
a = (1/2/pi)*(sum(SPP)+sum(SNN(2:end-1)))*df;
b = std(Z).^2 + mean(Z).*conj(mean(Z));

% c == d is true, for some real time series X
X = imag(Z);
[F,S]=mspec(X,[]);
S(1)=0.5*S(1);
c = (1/2/pi)*sum(S,1)*(F(2)-F(1));
d = std(X).^2 + mean(X).^2;
g = mean(X.*X);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Slab Layer Spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the winds
file = 'winds.nc';
u_wind = ncread(file, 'u');
v_wind = ncread(file, 'v');
time_wind = ncread(file, 't');
index_range = find(time_wind<max(t));
u_wind = u_wind(index_range);
v_wind = v_wind(index_range);
time_wind = time_wind(index_range);

% Mooring location
latitude = 24;

depth = 50;
slab_damp = 4;

[t, u, v] = OBLModel_DampedSlab( time_wind/86400, u_wind, v_wind, depth, latitude, slab_damp );

cv_wind = u + sqrt(-1)*v;
dt = time_wind(2)-time_wind(1);

[psi,lambda]=sleptap(size(cv_wind,1),3);
[omega,spp,snn,spn]=mspec(dt,cv_wind,psi);

f=omega*86400/(2*pi);

f_slab = [ flip(-f(2:end),1); f];
S_slab = [ flip(snn(2:end),1); spp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units', 'points', 'Position', [50 50 1000 400])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



plot(f_floats,S_floats)
ylog
hold on
h_slab = plot(f_slab,S_slab, 'k', 'LineWidth', 2);
set(gca,'FontSize', 16)
title('Drifter Spectra from the QG+Damped Slab Model', 'FontSize', 28.0, 'FontName', 'Helvetica')
xlabel('frequency (cycles per day)', 'FontSize', 20.0, 'FontName', 'Helvetica');
ylabel('power (m^2/s)', 'FontSize', 20.0, 'FontName', 'Helvetica');
xlim([-3 3])
ylim([2e-2 3e4])
legend(h_slab,{'Damped Slab Model'}, 'FontSize', 20.0, 'FontName', 'Helvetica')

export_fig(figure_out, '-r150')