addpath('../GLOceanKit/Matlab/')

%addpath('/Users/jearly/Dropbox/Documents/Matlab/jlab')
%load('/Users/jearly/Dropbox/Shared/Lilly-Sykulski-Early/MonopoleExperiment/QGDampedSlabTrajectories_Monopole.mat');

addpath('/Volumes/Music/Dropbox/Documents/Matlab/jlab')
%load('/Volumes/Music/Dropbox/Shared/Lilly-Sykulski-Early/MonopoleExperiment/QGDampedSlabTrajectories_Monopole.mat');
file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
latitude = ncreadatt(file, '/', 'latitude');
dRho1 = ncreadatt(file, '/', 'dRho1');
dRho2 = ncreadatt(file, '/', 'dRho2');
g1 = 9.81*dRho1;
g2 = 9.81*dRho2;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );

dt = t(2)-t(1);
timeIndex = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
psi2 = -(g2/f0)*squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
[strain_s_eul2, strain_n_eul2, rv_eul2, k, l] = FieldsFromStreamFunction( latitude, x, y, psi2, 'strain_s', 'strain_n', 'rv', 'k', 'l');

%%%%%%%%%%%%%%%%%%%%%
%
% Vorticity Plot
%
%%%%%%%%%%%%%%%%%%%%%%

theFigure = figure('Position', [50 50 1000 1000]);
theFigure.PaperPositionMode = 'auto';
theFigure.Color = 'white';

rvPlot = subplot(2,1,1);
theSSH = pcolor(x, y, rv_eul2/f0);
theSSH.EdgeColor = 'none';
axis(rvPlot, 'equal', 'tight');
rvPlot.Title.String = 'Upper layer vorticity';
rvPlot.XTick = [];
rvPlot.YTick = [];
colormap(rvPlot,gray(1024))
cb = colorbar('eastoutside');

xidx = 130;
yidx = 128;
hold on
scatter(x(xidx),y(yidx),'filled')
timeRange = find(t>=50*86400 & t<100*86400);

u1 = squeeze(ncread(file, 'u-1', [xidx yidx min(timeRange)], [1 1 length(timeRange)], [1 1 1]));
v1 = squeeze(ncread(file, 'v-1', [xidx yidx min(timeRange)], [1 1 length(timeRange)], [1 1 1]));
cv_mooring = u1 + sqrt(-1)*v1;
[psi,lambda]=sleptap(size(cv_mooring,1),1);
[omega_mooring,spp_mooring,snn_mooring,spn_mooring]=mspec(dt,cv_mooring,psi);
f=omega_mooring*86400/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double the zero frequency for plotting purposes
snn_mooring(1,:)=2*snn_mooring(1,:);
spp_mooring(1,:)=2*spp_mooring(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2);
plot(f,spp_mooring,'b'),ylog
hold on
plot(-f,snn_mooring,'b')


% Read the winds
file = 'winds.nc';
u_wind = ncread(file, 'u');
v_wind = ncread(file, 'v');
time_wind = ncread(file, 't');
index_range = find(time_wind<max(t(timeRange)) & time_wind>min(t(timeRange)) );
index_range = min(index_range):6:max(index_range);
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

vlines(-f0*86400/(2*pi))

xlim([-3 3])
