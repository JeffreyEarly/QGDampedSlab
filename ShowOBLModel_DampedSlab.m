addpath('../GLOceanKit/Matlab/')

% Read the winds
file = 'winds.nc';
u_wind = ncread(file, 'u');
v_wind = ncread(file, 'v');
time_wind = ncread(file, 't')/86400;

% Mooring location
latitude = 24;

depth = 50;
slab_damp = 4;

[t, u, v] = OBLModel_DampedSlab( time_wind, u_wind, v_wind, depth, latitude, slab_damp );

[t_stress, tau] = StressFromWindVector( time_wind*86400, u_wind, v_wind);

days_wind = time_wind - time_wind(1);
days_current = t - t(1);

figure,
plot(time_wind, u, 'b')
hold on
plot(time_wind, v, 'r')
xlim([0 40])

figure,
plot(time_wind, real(tau), 'b')
hold on
plot(time_wind, imag(tau), 'r')
xlim([0 40])

figure
subplot(2,1,1)
plot( days_wind, sqrt(u_wind.*u_wind + v_wind.*v_wind) )
xlabel('days')
ylabel('speed')
xlim([days_wind(1) days_wind(end)])

subplot(2,1,2)
plot( days_current, sqrt(u.*u + v.*v) )
xlabel('days')
ylabel('speed')
xlim([days_wind(1) days_wind(end)])