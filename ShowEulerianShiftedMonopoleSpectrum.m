addpath('../GLOceanKit/Matlab/')

addpath('/Users/jearly/Dropbox/Documents/Matlab/jlab')
load('/Users/jearly/Dropbox/Shared/Lilly-Sykulski-Early/MonopoleExperiment/QGDampedSlabTrajectories_Monopole.mat');
file = '/Volumes/Data/QGPlusSlab/MonopoleExperiment/QGDampedSlab_Monopole.nc';

%addpath('/Volumes/Music/Dropbox/Documents/Matlab/jlab')
%load('/Volumes/Music/Dropbox/Shared/Lilly-Sykulski-Early/MonopoleExperiment/QGDampedSlabTrajectories_Monopole.mat');
%file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';

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

timeRange = find(t>=250*86400 & t<300*86400);
%xidx = [10, 50, 110, 115, 120, 125, 128];
xidx = [25, 110, 125, 130];
yidx = 128*ones(size(xidx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
psi2 = -(g2/f0)*squeeze(ncread(file, 'eta-2', [1 1 min(timeRange)], [Inf Inf 1], [1 1 1]));
[strain_s_eul2, strain_n_eul2, rv_eul2, k, l] = FieldsFromStreamFunction( latitude, x, y, psi2, 'strain_s', 'strain_n', 'rv', 'k', 'l');

u1 = squeeze(ncread(file, 'u-1', [1 1 min(timeRange)], [Inf Inf 1], [1 1 1]));
v1 = squeeze(ncread(file, 'u-1', [1 1 min(timeRange)], [Inf Inf 1], [1 1 1]));

u1FD = fft2(u1);
v1FD = fft2(v1);
[K, L] = meshgrid(k,l);

rv_eul1 = double(ifft2((K.*v1FD - L.*u1FD)*2*pi*sqrt(-1),'symmetric'));

labels = cell(length(xidx),1);
for iIndex=1:length(xidx)
    labels{iIndex} = sprintf('%.2f f0', rv_eul1(xidx(iIndex),yidx(iIndex))/f0);
end

%%%%%%%%%%%%%%%%%%%%%
%
% Vorticity Plot (lower layer)
%
%%%%%%%%%%%%%%%%%%%%%%

theFigure = figure('Position', [50 50 800 1300]);
theFigure.PaperPositionMode = 'auto';
theFigure.Color = 'white';

rvPlot = subplot(3,1,1);
subrange = (length(y)/4):(3*length(y)/4);
theSSH = pcolor(x, y(subrange), rv_eul2(subrange,:)/f0);
theSSH.EdgeColor = 'none';
axis(rvPlot, 'equal', 'tight');
rvPlot.Title.String = 'Lower layer vorticity';
rvPlot.XTick = [];
rvPlot.YTick = [];
colormap(rvPlot,gray(1024))
cb = colorbar('eastoutside');

hold on
for iIndex=1:length(xidx)
    scatter(x(xidx(iIndex)),y(yidx(iIndex)),7*7,rvPlot.ColorOrder(mod(iIndex,7)+1,:),'filled')
end

%%%%%%%%%%%%%%%%%%%%%
%
% Vorticity Plot (upper layer)
%
%%%%%%%%%%%%%%%%%%%%%%

rvPlot2 = subplot(3,1,2);
subrange = (length(y)/4):(3*length(y)/4);
theSSH2 = pcolor(x, y(subrange), rv_eul1(subrange,:)/f0);
theSSH2.EdgeColor = 'none';
axis(rvPlot2, 'equal', 'tight');
rvPlot2.Title.String = 'Upper layer vorticity';
rvPlot2.XTick = [];
rvPlot2.YTick = [];
colormap(rvPlot2,gray(1024))
cb = colorbar('eastoutside');

hold on
for iIndex=1:length(xidx)
    scatter(x(xidx(iIndex)),y(yidx(iIndex)),7*7,rvPlot.ColorOrder(mod(iIndex,7)+1,:),'filled')
end

%%%%%%%%%%%%%%%%%%%%%
%
% Eulerian spectrum plot
%
%%%%%%%%%%%%%%%%%%%%%%


if ~exist('cv_mooring','var') || ~isequal(size(cv_mooring) , [length(timeRange) length(xidx)])
    cv_mooring = zeros(length(timeRange),length(xidx));
    for iIndex=1:length(xidx) 
        u1 = squeeze(ncread(file, 'u-1', [xidx(iIndex) yidx(iIndex) min(timeRange)], [1 1 length(timeRange)], [1 1 1]));
        v1 = squeeze(ncread(file, 'v-1', [xidx(iIndex) yidx(iIndex) min(timeRange)], [1 1 length(timeRange)], [1 1 1]));
        cv_mooring(:,iIndex) = u1 + sqrt(-1)*v1;
    end
end
[psi,lambda]=sleptap(size(cv_mooring,1),1);
[omega_mooring,spp_mooring,snn_mooring,spn_mooring]=mspec(dt,cv_mooring,psi);
f=omega_mooring*86400/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double the zero frequency for plotting purposes
snn_mooring(1,:)=2*snn_mooring(1,:);
spp_mooring(1,:)=2*spp_mooring(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,1,3);
plot(f,spp_mooring),ylog
hold on
plot(-f,snn_mooring)

vlines(-f0*86400/(2*pi))

xlim([-2 2])
ylim([1e-2 1e6])
legend(labels,'FontSize',12)
