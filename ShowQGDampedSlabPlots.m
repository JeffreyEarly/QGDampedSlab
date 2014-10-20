file = '/Users/jearly/Desktop/QGDampedSlab.nc';
t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');

u = squeeze(ncread(file, 'u', [1 1 1], [1 1 Inf], [1 1 1]));
v = squeeze(ncread(file, 'v', [1 1 1], [1 1 Inf], [1 1 1]));

tau_x = squeeze(ncread(file, 'tau_x'));
tau_y = squeeze(ncread(file, 'tau_y'));

figure,
plot(t, u, 'b')
hold on
plot(t, v, 'r')

h = double(squeeze(ncread(file, 'eta1', [1 1 70], [Inf Inf 1], [1 1 1])));
figure, pcolor(x,y,h), shading flat