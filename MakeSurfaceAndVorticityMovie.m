addpath('../GLOceanKit/Matlab/')

file = '/Volumes/Data/QGPlusSlab/MonopoleExperiment/QGDampedSlab_Monopole.nc';
FramesFolder ='/Volumes/Data/QGPlusSlab/MonopoleExperiment/SurfaceVorticityFrames';

file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';
FramesFolder ='/Volumes/Music/Model_Output/MonopoleExperiment/SurfaceVorticityFrames';

file = '/Volumes/RadiativeTransfer/QGDampedSlab/MonopoleExperiment/QGDampedSlab_Monopole.nc';
FramesFolder ='/Volumes/RadiativeTransfer/QGDampedSlab/MonopoleExperiment/SurfaceVorticityFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');
xFloat = ncread(file, 'x-float');
yFloat = ncread(file, 'y-float');
latitude = ncreadatt(file, '/', 'latitude');
dRho1 = ncreadatt(file, '/', 'dRho1');
dRho2 = ncreadatt(file, '/', 'dRho2');
g1 = 9.81*dRho1;
g2 = 9.81*dRho2;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );
stride = 1;
layer1floatSize = 10;
layer2floatSize = 15;

xlength = max(x)-min(x) + x(2)-x(1);
ylength = max(y)-min(y) + y(2)-y(1);

[X,Y] = meshgrid(x,y);

xPosition1Initial = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
xPosition2Initial = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
xPosition1Initial = (reshape(xPosition1Initial, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
xPosition2Initial = (reshape(xPosition2Initial, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';

psi2 = -(g2/f0)*squeeze(ncread(file, 'eta-2', [1 1 length(t)], [Inf Inf 1], [1 1 1]));
[rv_eul2, k, l] = FieldsFromStreamFunction( latitude, x, y, psi2, 'rv', 'k', 'l');
rv_eul2 = rv_eul2/f0;
maxRV2 = max(max(rv_eul2));
minRV2 = min(min(rv_eul2));
[K, L] = meshgrid(k,l);
u1 = squeeze(ncread(file, 'u-1', [1 1 length(t)], [Inf Inf 1], [1 1 1]));
v1 = squeeze(ncread(file, 'v-1', [1 1 length(t)], [Inf Inf 1], [1 1 1]));
u1FD = fft2(u1);
v1FD = fft2(v1);
rv_eul1 = double(ifft2((K.*v1FD - L.*u1FD)*2*pi*sqrt(-1),'symmetric'));
rv_eul1 = rv_eul1/f0;
maxRV1 = max(max(rv_eul1));
minRV1 = min(min(rv_eul1));

% Set them to the same scale? I prefer they be slightly saturated, rather
% than span the color space.
minRV2 = min(minRV2,minRV1);
minRV1 = minRV2;
maxRV2 = max(maxRV2,maxRV1);
maxRV1 = maxRV2;

rv1offset = 1.001*(maxRV1 - minRV1);
cbmin = minRV2;
cbmax = rv1offset + maxRV2;

    
colorposition = [xPosition2Initial, (max(max(xPosition2Initial))-min(min(xPosition2Initial)))+xPosition1Initial];
%graymap = [linspace(.35,.85,128)', linspace(.35,.85,128)', linspace(.35,.85,128)'];
layer2map = gray(128);
layer1map = parula(128);
combomap = [layer2map;layer1map];
floatSize = [ layer2floatSize*layer2floatSize*ones(size(xPosition2Initial)), layer1floatSize*layer1floatSize*ones(size(xPosition1Initial))];

fig = figure('Position', [50 50 1080 1080]);
fig.PaperPositionMode = 'auto';
fig.Color = 'w';

%for timeIndex=1:2:length(t)
    for timeIndex = 6001:6001;

    mainPlot = subplot(1,10,1:8);
    
    xPosition1 = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
    yPosition1 = squeeze(ncread(file, 'y-position-layer-1', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));

    % Reshape to [time, float]
    xpos1 = (reshape(xPosition1, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    ypos1 = (reshape(yPosition1, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    xpos1 = mod( xpos1-min(x), xlength ) + min(x);
    ypos1 = mod( ypos1-min(y), ylength ) + min(y);
    
    % Determine the relative vorticity of the upper layer
    u1 = squeeze(ncread(file, 'u-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
    v1 = squeeze(ncread(file, 'v-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
    u1FD = fft2(u1);
    v1FD = fft2(v1);
    rv_eul1 = double(ifft2((K.*v1FD - L.*u1FD)*2*pi*sqrt(-1),'symmetric'))/f0;
    rv1 = interp2( X, Y, rv_eul1, xpos1, ypos1 );
    
    xPosition2 = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
    yPosition2 = squeeze(ncread(file, 'y-position-layer-2', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));

    % Reshape to [time, float]
    xpos2 = (reshape(xPosition2, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    ypos2 = (reshape(yPosition2, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    xpos2 = mod( xpos2-min(x), xlength ) + min(x);
    ypos2 = mod( ypos2-min(y), ylength ) + min(y);
    
    psi2 = -(g2/f0)*squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
    [rv_eul2] = FieldsFromStreamFunction( latitude, x, y, psi2, 'rv');
    rv2 = interp2( X, Y, rv_eul2/f0, xpos2, ypos2 );
    
   % default color map is only 128 shades---we need more!
	xpos = [xpos2, xpos1];
    ypos = [ypos2, ypos1];
   
    
    colorposition = [rv2, rv1offset + rv1];
    %mesh([xpos;xpos],[ypos;ypos],[colorposition;colorposition],'mesh','column','marker','.','MarkerSize', floatSize), view(2)
    scatter(xpos, ypos, floatSize, colorposition,'filled','MarkerEdgeColor',[.0 .0 .0],'LineWidth',0.25);
    
    colormap(combomap)
    caxis([cbmin cbmax])
    
    % make the axes look better
	set( gca, 'TickDir', 'out');
	set( gca, 'Linewidth', 1.0);
	axis equal tight
	
	% get rid of the xticks because we're going to put a colorbar below with the same info.
% 	set( gca, 'xtick', [])
%     set( gca, 'xlabel', [])
%     set( gca, 'ytick', [])
%     set( gca, 'ylabel', [])
    axis off
    
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
    
	% label everything
	title( sprintf('Floats advected by a Quasigeostrophic eddy with wind'), 'fontsize', 28, 'FontName', 'Helvetica' );
    text( 1e5, -4.7e5,  sprintf('Day %d @ %2d:00',floor(t(timeIndex)/86400), round(mod(t(timeIndex),86400)/3600)), 'fontsize', 28, 'FontName', 'Helvetica', 'BackgroundColor', 'white' )
    
    cb1Plot = subplot(1,10,9);
    ycb = linspace(minRV2, maxRV2, 128)';
    xcb = (1:10)';
    cb1 = repmat(ycb,[1,10]);
    pcolor(xcb,ycb,cb1)
    shading interp
    cb1Plot.XTick = [];
    caxis([cbmin cbmax])
    for iIndex=1:length(cb1Plot.YTick)
        cb1Plot.YTickLabel{iIndex}=sprintf('%.2f f_0',cb1Plot.YTick(iIndex));
    end
    
    cb2Plot = subplot(1,10,10);
    ycb = linspace(rv1offset + minRV1, rv1offset + maxRV1, 128)';
    xcb = (1:10)';
    cb2 = repmat(ycb,[1,10]);
    pcolor(xcb,ycb,cb2)
    shading interp
    cb2Plot.XTick = [];
    cb2Plot.YTick = [];
    caxis([cbmin cbmax])
%     for iIndex=1:length(cb2Plot.YTick)
%         cb2Plot.YTickLabel{iIndex}=sprintf('%.2f f_0',cb2Plot.YTick(iIndex)-rv1offset);
%     end
    
    d = fig.PaperPosition;
    fig.PaperSize = [d(3) d(4)];

%     ti = mainPlot.TightInset;
%     mainPlot.Position=[-.12 ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)];
%     cb1Plot.Position = [mainPlot.Position(1)+mainPlot.Position(2) ti(2) cb1Plot.Position(3) 1-ti(4)-ti(2)];

   
   
	% write everything out	
	output = sprintf('%s/Hour_%05d', FramesFolder,timeIndex-1);
	%print(fig, '-dpsc2', output)
    print(fig, '-dpng', '-r150', output)

end