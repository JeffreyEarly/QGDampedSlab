addpath('/Users/jearly/Dropbox/Documents/Matlab/jlab')
addpath('../GLOceanKit/Matlab/')

folder = '/Volumes/Data/QGPlusSlab/TurbulenceExperimentNonStiff';
file = 'QGDampedSlab.nc';
file = sprintf('%s/%s',folder,file);
framesFolder = sprintf('%s/DisplacementFrames',folder);

if exist(framesFolder, 'dir') == 0
	mkdir(framesFolder);
end


day = 300;
t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');

t_days = t/86400;

timeIndex = find( t_days-1 <= day, 1, 'last');

eta1 = squeeze(ncread(file, 'eta-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
eta2 = squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));

eta1ColorLimits = 0.8*[-max(max(abs(eta1))) max(max(abs(eta1)))];
eta2ColorLimits = 0.8*[-max(max(abs(eta2))) max(max(abs(eta2)))];

theFigure = figure('Position', [50 50 1000 470]);
theFigure.PaperPositionMode = 'auto';
theFigure.Color = 'white';

for timeIndex=545:length(t_days)

    eta1 = squeeze(ncread(file, 'eta-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
    eta2 = squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
    
    if (max(max(abs(eta1))) < 0.003*eta1ColorLimits(2))
        eta1 = 1000*eta1;
        units = 'millimeters';
    elseif (max(max(abs(eta1))) < 0.03*eta1ColorLimits(2))
        eta1 = 100*eta1;
        units = 'centimeters';
    elseif (max(max(abs(eta1))) < 0.3*eta1ColorLimits(2))
        eta1 = 10*eta1;
        units = 'decimeters';
    else
        units = 'meters';
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %
    % eta1 Plot
    %
    %%%%%%%%%%%%%%%%%%%%%%

    sshPlot = subplot(1,2,1);
    theSSH = pcolor(x, y, eta1);
    theSSH.EdgeColor = 'none';
    axis(sshPlot, 'equal', 'tight');
    sshPlot.Title.String = 'Upper layer height';
    sshPlot.Title.FontName = 'Helvetica';
    sshPlot.Title.FontSize = 28;
    sshPlot.XTick = [];
    sshPlot.YTick = [];
    sshPlot.CLim = eta1ColorLimits;
    colormap(sshPlot,parula(1024))

    sshBar = colorbar('westoutside');
    sshBar.Label.String = sprintf('Displacement (%s)',units);
    sshBar.Label.FontName = 'Helvetica';
    sshBar.Label.FontSize = 16;

    %%%%%%%%%%%%%%%%%%%%%
    %
    % eta2 Plot
    %
    %%%%%%%%%%%%%%%%%%%%%%

    ssh2Plot = subplot(1,2,2);
    theSSH = pcolor(x, y, eta2);
    theSSH.EdgeColor = 'none';
    axis(ssh2Plot, 'equal', 'tight');
    ssh2Plot.Title.String = 'Lower layer height';
    ssh2Plot.Title.FontName = 'Helvetica';
    ssh2Plot.Title.FontSize = 28;
    ssh2Plot.XTick = [];
    ssh2Plot.YTick = [];
    ssh2Plot.CLim = eta2ColorLimits;
    colormap(ssh2Plot,parula(1024))

    ssh2Bar = colorbar('eastoutside');
    ssh2Bar.Label.String = 'Displacement (meters)';
    ssh2Bar.Label.FontName = 'Helvetica';
    ssh2Bar.Label.FontSize = 16;

    text( 1.5e5, -4.5e5,  sprintf('Day %d @ %2d:00',floor(t(timeIndex)/86400), round(mod(t(timeIndex),86400)/3600)), 'fontsize', 16, 'FontName', 'Helvetica', 'BackgroundColor', 'white' )
    packfig(1,2)
    
    output = sprintf('%s/t_%03d', framesFolder,timeIndex-1);
    export_fig(output,'-r300')
end