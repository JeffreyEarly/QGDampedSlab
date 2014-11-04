file = '/Users/jearly/Desktop/QGDampedSlab.nc';
t = ncread(file, 'time');
x = ncread(file, 'x');
y = ncread(file, 'y');

u = squeeze(ncread(file, 'u', [1 1 1], [1 1 Inf], [1 1 1]));
v = squeeze(ncread(file, 'v', [1 1 1], [1 1 Inf], [1 1 1]));

tau_x = squeeze(ncread(file, 'tau_x'));
tau_y = squeeze(ncread(file, 'tau_y'));

xPosition1 = squeeze(ncread(file, 'x-position-layer-1', [1 1 1], [1 1 Inf], [1 1 1]));
yPosition1 = squeeze(ncread(file, 'y-position-layer-1', [1 1 1], [1 1 Inf], [1 1 1]));

xPosition2 = squeeze(ncread(file, 'x-position-layer-2', [1 1 1], [1 1 Inf], [1 1 1]));
yPosition2 = squeeze(ncread(file, 'y-position-layer-2', [1 1 1], [1 1 Inf], [1 1 1]));

eta2 = squeeze(ncread(file, 'eta2', [1 1 800], [Inf Inf 1], [1 1 1]));
figure
pcolor(x,y,eta2), axis equal tight, shading interp

figure
plot(xPosition2,yPosition2)

figure
plot(xPosition1,yPosition1)
dt = t(3)-t(2);
u1 = diff(xPosition1)/dt;
v1 = diff(yPosition1)/dt;

t_days = t/86140;

figure,
plot(t_days(2:end), u1, 'b')
hold on
plot(t_days(2:end), v1, 'r')
xlim([0 40])

figure,
plot(t_days, u, 'b')
hold on
plot(t_days, v, 'r')
xlim([0 40])

% h = double(squeeze(ncread(file, 'eta1', [1 1 70], [Inf Inf 1], [1 1 1])));
% figure, pcolor(x,y,h), shading flat

figure
plot(t_days,tau_x, 'b')
hold on
plot(t_days, tau_y, 'r')
xlim([0 40])