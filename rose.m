n = 800;
a = 1.995653;
b = 1.27689;
c = 8;
r = linspace(0,1,n);
th = linspace(-2,20*pi,n);

[R,TH] = ndgrid(r,th);
petalNum = 3.6;
x = 1 - 0.5 * ((5/4) * (1 - mod(petalNum*TH,2*pi)/pi).^2 - 0.25) .^ 2;
phi = (pi/2) * exp(-TH/(c*pi));
y = a * (R.^2) .* (b*R-1).^2 .* sin(phi);
R2 = x .* (R.*sin(phi) + y.*cos(phi));

X = R2 .* sin(TH);
Y = R2 .* cos(TH);
Z = x .* (R.*cos(phi) - y.*sin(phi));

red_map = linspace(1,0.25,10)';
red_map(:,2) = 0;
red_map(:,3) = 0;

figure(100)
clf
surf(X,Y,Z,'LineStyle','none')
view([-40.50 42.00])
colormap(red_map)