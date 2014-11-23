addpath('/Volumes/Music/Dropbox/Documents/Matlab/jlab')
addpath('../GLOceanKit/Matlab/')


load('/Volumes/Music/Model_Output/QGDampedSlabTrajectories_Monopole.mat');

numDrifters = size(xpos1,2);
days = t/86400;
dt = t(2)-t(1);

cv1 = (diff(xpos1,1,1) + sqrt(-1)*diff(ypos1,1,1))/dt;
cv2 = (diff(xpos2,1,1) + sqrt(-1)*diff(ypos2,1,1))/dt;

[psi,lambda]=sleptap(size(cv1,1),3);
[omega,spp,snn,spn]=mspec(dt,cv1,psi);

% convert from radians/second to cycles/day
f=omega*86400/(2*pi);

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
plot(f,spp, 'k', 'LineWidth', 2)
plot(-f,snn, 'k', 'LineWidth', 2)

xlim([-3 3])
ylim([2e-2 3e4])