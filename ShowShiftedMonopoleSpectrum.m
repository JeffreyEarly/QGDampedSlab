%addpath('/Volumes/Music/Dropbox/Documents/Matlab/jlab')
addpath('/Users/jearly/Dropbox/Documents/Matlab/jlab')
addpath('../GLOceanKit/Matlab/')

load('/Volumes/Data/QGPlusSlab/MonopoleExperimentWithNoWindsLongDamp_+/QGDampedSlabTrajectories_Monopole.mat');
%load('/Volumes/Music/Dropbox/Shared/Lilly-Sykulski-Early/MonopoleExperiment/QGDampedSlabTrajectories_Monopole.mat');



numDrifters = size(xpos1,2);
days = t/86400;
dt = t(2)-t(1);

cv1 = (diff(xpos1,1,1) + sqrt(-1)*diff(ypos1,1,1))/dt;
cv2 = (diff(xpos2,1,1) + sqrt(-1)*diff(ypos2,1,1))/dt;

% find a drifter with a large shift
% [val, idx]=min(rv1(1,:));
% idx = 136;
timeRange = find(days>0 & days<5);

% cv = cv1(timeRange,idx);
cv = cv1(timeRange,:);

[psi,lambda]=sleptap(size(cv,1),1);
[omega,spp,snn,spn]=mspec(dt,cv,psi);

% convert from radians/second to cycles/day
f=omega*86400/(2*pi);

f0 = corfreq(24)/3600;
% zeta = rv1(:,idx);
% sigma2 = strain_n1(:,idx).^2 + strain_s1(:,idx).^2;
% omega_shift = sqrt( (f0 + zeta/2).^2 - sigma2/4 );
% 
% figure
% plot(days, omega_shift/f0)
% 
% figure
% plot(days, omega_shift*86400/(2*pi))

% file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';
% u1 = squeeze(ncread(file, 'u-1', [135 128 1], [1 1 max(timeRange)], [1 1 1]));
% v1 = squeeze(ncread(file, 'v-1', [135 128 1], [1 1 max(timeRange)], [1 1 1]));
% cv_mooring = u1 + sqrt(-1)*v1;
% [omega_mooring,spp_mooring,snn_mooring,spn_mooring]=mspec(dt,cv_mooring,psi);
% 
% eta2 = squeeze(ncread(file, 'eta-2', [1 1 1], [Inf Inf 1], [1 1 1]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Units', 'points', 'Position', [50 50 1000 400])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double the zero frequency for plotting purposes
snn(1,:)=2*snn(1,:);
spp(1,:)=2*spp(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(f,spp),ylog
hold on
plot(-f,snn)
vlines([f0*86400/(2*pi)])

% % Read the winds
% file = 'winds.nc';
% u_wind = ncread(file, 'u');
% v_wind = ncread(file, 'v');
% time_wind = ncread(file, 't');
% index_range = find(time_wind<t(length(cv)));
% u_wind = u_wind(index_range);
% v_wind = v_wind(index_range);
% time_wind = time_wind(index_range);
% 
% % Mooring location
% latitude = 24;
% 
% depth = 50;
% slab_damp = 4;
% 
% [t, u, v] = OBLModel_DampedSlab( time_wind/86400, u_wind, v_wind, depth, latitude, slab_damp );
% 
% cv_wind = u + sqrt(-1)*v;
% dt = time_wind(2)-time_wind(1);
% 
% [psi,lambda]=sleptap(size(cv_wind,1),1);
% [omega,spp,snn,spn]=mspec(dt,cv_wind,psi);
% 
% f=omega*86400/(2*pi);
% plot(f,spp, 'k', 'LineWidth', 2)
% plot(-f,snn, 'k', 'LineWidth', 2)
% 
% xlim([-3 3])
% ylim([2e-2 3e4])