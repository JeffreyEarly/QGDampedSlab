file = '/Volumes/Data/QGPlusSlab/MonopoleExperiment/QGDampedSlab_Monopole.nc';
FramesFolder ='/Volumes/Data/QGPlusSlab/MonopoleExperiment/FloatFrames';

file = '/Volumes/Music/Model_Output/MonopoleExperiment/QGDampedSlab_Monopole.nc';
FramesFolder ='/Volumes/Music/Model_Output/MonopoleExperiment/FloatFrames';

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

stride = 2;
layer1floatSize = 10;
layer2floatSize = 15;

[X,Y] = meshgrid(x,y);

xPosition1Initial = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
xPosition2Initial = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
xPosition1Initial = (reshape(xPosition1Initial, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
xPosition2Initial = (reshape(xPosition2Initial, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';

colorposition = [xPosition2Initial, (max(max(xPosition2Initial))-min(min(xPosition2Initial)))+xPosition1Initial];
graymap = [linspace(.35,.85,128)', linspace(.35,.85,128)', linspace(.35,.85,128)'];
combomap = [graymap;jet(128)];
floatSize = [ layer2floatSize*layer2floatSize*ones(size(xPosition2Initial)), layer1floatSize*layer1floatSize*ones(size(xPosition1Initial))];

fig = figure('Position', [50 50 1080 1080]);
fig.PaperPositionMode = 'auto';
fig.Color = 'w';

%for timeIndex=1:2:length(t)
    for timeIndex = 1:1;

%     subplot(1,2,1)
    
    xPosition1 = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
    yPosition1 = squeeze(ncread(file, 'y-position-layer-1', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));

    % Reshape to [time, float]
    xpos1 = (reshape(xPosition1, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    ypos1 = (reshape(yPosition1, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';

    xPosition2 = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));
    yPosition2 = squeeze(ncread(file, 'y-position-layer-2', [ceil(stride/2) ceil(stride/2) timeIndex], [length(yFloat)/stride length(xFloat)/stride 1], [stride stride 1]));

    % Reshape to [time, float]
    xpos2 = (reshape(xPosition2, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    ypos2 = (reshape(yPosition2, [length(yFloat)*length(xFloat)/(stride*stride), 1]))';
    
   % default color map is only 128 shades---we need more!
	xpos = [xpos2, xpos1];
    ypos = [ypos2, ypos1];
   
    %mesh([xpos;xpos],[ypos;ypos],[colorposition;colorposition],'mesh','column','marker','.','MarkerSize', floatSize), view(2)
    scatter(xpos, ypos, floatSize, colorposition,'filled','MarkerEdgeColor',[.0 .0 .0]);
    
    colormap(combomap)
    
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
    text( 2e5, -4.7e5,  sprintf('Day %d @ %2d:00',floor(t(timeIndex)/86400), round(mod(t(timeIndex),86400)/3600)), 'fontsize', 28, 'FontName', 'Helvetica', 'BackgroundColor', 'white' )
%     subplot(1,2,2)
%     
%     colormap default
%     
%     eta1 = squeeze(ncread(file, 'eta-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
% 
%     [nx, ny, nz] = surfnorm(X,Y,eta1);
%     b = reshape([nx ny nz], length(x),length(y),3);
% 
%     eta1plot = surf(X,Y,eta1,'VertexNormals',b,'EdgeColor','none');
%     eta1axes = gca;
% 
%     lighting gouraud
%     camlight
%     %grid off
%     %axis off
%     % alpha(eta1plot,0.9)
% 
%     eta1axes.XTickLabel = [];
%     eta1axes.YTickLabel = [];
%     eta1axes.XTick = [];
%     eta1axes.YTick = [];
% 
%     xlim([min(x) max(x)])
%     ylim([min(y) max(y)])
%     zlim([-28 28])
%     caxis([-5 5])
% 
%     view(28,14)
    
   d = fig.PaperPosition;
   fig.PaperSize = [d(3) d(4)];
        
	% write everything out	
	output = sprintf('%s/Hour_%05d', FramesFolder,timeIndex-1);
	print(fig, '-dpsc2', output)

end