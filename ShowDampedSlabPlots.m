file = '/Users/jearly/Desktop/WindForcedFPlane2.nc';
t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');

u = squeeze(ncread(file, 'u', [1 1 1], [1 1 Inf], [1 1 1]));
v = squeeze(ncread(file, 'v', [1 1 1], [1 1 Inf], [1 1 1]));

figure,
plot(t, u, 'b')
hold on
plot(t, v, 'r')
xlim([0 40])

t = t*86400;
t = t(1:1205);
u = u(1:1205);
v = v(1:1205);
save('/Users/jearly/Desktop/dampedslab.mat','u', 'v', 't')

% h = double(squeeze(ncread(file, 'h', [1 1 70], [Inf Inf 1], [1 1 1])));
% figure, pcolor(x,y,h), shading flat