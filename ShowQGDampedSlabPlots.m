day = 300;

%file = '/Users/jearly/Desktop/QGDampedSlab.nc';
file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';
file = '/Volumes/Music/Model_Output/TurbulenceExperimentLongDampNonStiff/QGDampedSlab.nc';
%output = '/Volumes/Music/Model_Output/QGDampedSlabTrajectories_Monopole.mat';
t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');

t_days = t/86400;
timeIndex = find( t_days <= day, 1, 'last');

% u = squeeze(ncread(file, 'u', [1 1 1], [1 1 timeIndex], [1 1 1]));
% v = squeeze(ncread(file, 'v', [1 1 1], [1 1 timeIndex], [1 1 1]));
% 
% tau_x = squeeze(ncread(file, 'tau_x'));
% tau_y = squeeze(ncread(file, 'tau_y'));

stride = 4;

xFloat = ncread(file, 'x-float');
yFloat = ncread(file, 'y-float');
xPosition1 = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));
yPosition1 = squeeze(ncread(file, 'y-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));

% Reshape to [time, float]
xpos1 = (reshape(xPosition1, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';
ypos1 = (reshape(yPosition1, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';

xPosition2 = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));
yPosition2 = squeeze(ncread(file, 'y-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));

% Reshape to [time, float]
xpos2 = (reshape(xPosition2, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';
ypos2 = (reshape(yPosition2, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';

eta1 = squeeze(ncread(file, 'eta-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
eta2 = squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));

%save(output, 't', 'xpos1', 'ypos1', 'xpos2', 'ypos2')

theFigure = figure('Position', [50 50 1000 1000]);
theFigure.PaperPositionMode = 'auto';
theFigure.Color = 'white';


%%%%%%%%%%%%%%%%%%%%%
%
% eta1 Plot
%
%%%%%%%%%%%%%%%%%%%%%%

sshPlot = subplot(2,2,1);
theSSH = pcolor(x, y, eta1);
theSSH.EdgeColor = 'none';
axis(sshPlot, 'equal', 'tight');
sshPlot.Title.String = 'Upper layer height';
sshPlot.XTick = [];
sshPlot.YTick = [];
colormap(sshPlot,gray(1024))

%%%%%%%%%%%%%%%%%%%%%
%
% eta2 Plot
%
%%%%%%%%%%%%%%%%%%%%%%

ssh2Plot = subplot(2,2,2);
theSSH = pcolor(x, y, eta2);
theSSH.EdgeColor = 'none';
axis(ssh2Plot, 'equal', 'tight');
ssh2Plot.Title.String = 'Lower layer height';
ssh2Plot.XTick = [];
ssh2Plot.YTick = [];
colormap(ssh2Plot,gray(1024))

%%%%%%%%%%%%%%%%%%%%%
%
% Upper layer float Plot
%
%%%%%%%%%%%%%%%%%%%%%%

ulFloatPlot = subplot(2,2,3);
plot(xpos1, ypos1)
ulFloatPlot.Title.String = 'Upper layer floats';
xlim([min(min(xpos1)) max(max(xpos1))])
ylim([min(min(ypos1)) max(max(ypos1))])

%%%%%%%%%%%%%%%%%%%%%
%
% Lower layer float Plot
%
%%%%%%%%%%%%%%%%%%%%%%

llFloatPlot = subplot(2,2,4);
plot(xpos2, ypos2)
llFloatPlot.Title.String = 'Lower layer floats';
xlim([min(min(xpos2)) max(max(xpos2))])
ylim([min(min(ypos2)) max(max(ypos2))])

% figure
% plot(xPosition2,yPosition2)
% 
% figure
% plot(xPosition1,yPosition1)
% dt = t(3)-t(2);
% u1 = diff(xPosition1)/dt;
% v1 = diff(yPosition1)/dt;

% t_days = t/86140;

% figure,
% plot(t_days(2:end), u1, 'b')
% hold on
% plot(t_days(2:end), v1, 'r')
% xlim([0 40])

% figure,
% plot(t_days, u, 'b')
% hold on
% plot(t_days, v, 'r')
% xlim([0 40])

% h = double(squeeze(ncread(file, 'eta1', [1 1 70], [Inf Inf 1], [1 1 1])));
% figure, pcolor(x,y,h), shading flat

% figure
% plot(t_days,tau_x, 'b')
% hold on
% plot(t_days, tau_y, 'r')
% xlim([0 40])