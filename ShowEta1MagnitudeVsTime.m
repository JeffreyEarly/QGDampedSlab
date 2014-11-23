file = '/Volumes/Music/Model_Output/QGDampedSlab_Monopole.nc';
addpath('../GLOceanKit/Matlab/')

[t,f0] = FieldsFromTurbulenceFileWithNamedSSH(file, 1, 'eta-1', 't', 'f0');

maxEta = zeros(size(t));

for timeIndex=1:length(t)
    maxEta(timeIndex) = max(max(abs(squeeze(ncread(file, 'eta-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1])))));
end

t_inertial = t*f0/(2*pi);

figure
etaPlot = plot(t_inertial,maxEta);
title('Upper layer height (m) vs. inertial periods');
xlabel('time (inertial periods)')
ylabel('maximum displacement (m)')
xlim([min(t_inertial) max(t_inertial)])