% clear;
close all;

sim_name = 'twocell';

base_dir = pwd;
data_dir = ['~/Documents/tables/les/', sim_name];

v0 = 100; % This is the characteristic velocity, which scales the non-dimensional model results.  
% Common values for the two-cell vortex are 150 and 225 m/s.
% Suggested value for the suction vortex case is 250 m/s.
g = 9.8;

% LES file name
LES_name = 'LES_mean_1_6_fnum1.dat';

num_times = 10;

% Spatial dimensions are scaled by a factor of v0^2/9.8, so increasing the
% non-dimensional velocity scale increases the domain size.  Time also
% increases with increasing velocity scale by a factor of v0/g.

%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the LES grid %
%%%%%%%%%%%%%%%%%%%%%%%%

% cd(data_dir)
fid = fopen([data_dir, '/fort.10_2'], 'r', 'ieee-le');
fread(fid, 1, 'int'); % junk
Nx = fread(fid, 1, 'int');
Ny = fread(fid, 1, 'int');
Nz = fread(fid, 1, 'int');
fread(fid, 2, 'real*4'); % junk
Xc = fread(fid, Nx*Ny*Nz, 'real*4');
fread(fid, 2, 'real*4'); % junk
Yc = fread(fid, Nx*Ny*Nz, 'real*4');
fread(fid, 2, 'real*4'); % junk
Zc = fread(fid, Nx*Ny*Nz, 'real*4');

% Dimensionalize variables
Xc = squeeze(Xc(1:Nx));
Yc = Xc;
Zc = squeeze(Zc(1+(Nx-1)*(Ny-1) : (Nx-1)*(Ny-1) : (Nz)*(Nx-1)*(Ny-1)+1));
Xc = Xc * v0^2 / 9.80;
Yc = Yc * v0^2 / 9.80;
Zc = Zc * v0^2 / 9.80;

% Original values of Xm, Ym, Zm are the grid cell edges.  So, we need the
% position of the grid cell centers.
% xc = (Xc(2:max(size(Xc)))+Xc(1:max(size(Xc))-1))/2;
% yc(1,1:Ny-1) = (Yc(2:max(size(Yc)))+Yc(1:max(size(Yc))-1))/2;
% zc(1,1,1:Nz-1) = (Zc(2:max(size(Zc)))+Zc(1:max(size(Zc))-1))/2;
xc = Xc(1:max(size(Xc)));
yc = permute(Yc(1:max(size(Yc))), [2 1]);
zc = permute(Zc(1:max(size(Zc))), [3 2 1]);
Xm = repmat(xc, [1 Ny Nz]);
Ym = repmat(yc, [Nx 1  Nz]);
Zm = repmat(zc, [Nx Ny 1]);
fclose(fid);

% LES data is only stored in a subdomain where the tornado is located
ix1 = 1;
ix2 = Nx;
iy1 = 1;
iy2 = Ny;
iz1 = 1;
iz2 = Nz;
% ix1 = 15;
% ix2 = Nx - 16; % minimum and maximum (x,y,z) values of indices saved
% iy1 = 15;
% iy2 = Nx - 16;
% iz1 = 1;
% iz2 = 51;

Xmf = Xm(ix1:ix2, iy1:iy2, iz1:iz2); % Reconstruct grid only for saved data pts
Ymf = Ym(ix1:ix2, iy1:iy2, iz1:iz2);
Zmf = Zm(ix1:ix2, iy1:iy2, iz1:iz2);

r = sqrt(Xmf.^2 + Ymf.^2);
phi = atan(Xmf ./ Ymf); % Azimuth angles w.r.t TOR center (domain center)
% adjust angle values for each quadrant
phi(Ymf < 0) = phi(Ymf < 0) + pi; % Lower 2 quadrants
phi(Xmf < 0 & Ymf >= 0) = phi(Xmf < 0 & Ymf >= 0) + 2*pi; % Upper left quadrant

%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the LES file %
%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([data_dir, '/', LES_name], 'r', 'ieee-le');
fread(fid, 1, 'int'); % junk

time = zeros(1,num_times);
ustore = zeros(size(Xmf,1), size(Xmf,2), size(Xmf,3), num_times);
vstore = ustore;
wstore = ustore;
pstore = ustore;
tkestore = ustore;
urstore = ustore;
vtstore = ustore;
Lstore = ustore;

for kdx = 1:num_times
    time(kdx) = fread(fid, 1, 'real*4') * v0 / g;
    fread(fid, 2, 'int'); % junk
    u = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0;
    fread(fid, 2, 'int'); % junk
    v = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0;
    fread(fid, 2, 'int'); % junk
    w = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0;
    fread(fid, 2, 'int'); % junk
    p = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0^2;
    fread(fid, 2, 'int'); % junk
    tke = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0^2;
    fread(fid, 2, 'int'); % junk
    
    u = reshape(u, size(Xmf,1), size(Xmf,2), size(Xmf,3));
    v = reshape(v, size(Xmf,1), size(Xmf,2), size(Xmf,3));
    w = reshape(w, size(Xmf,1), size(Xmf,2), size(Xmf,3));
    p = reshape(p, size(Xmf,1), size(Xmf,2), size(Xmf,3));
    tke = reshape(tke, size(Xmf,1), size(Xmf,2), size(Xmf,3));
    

    ur = u.*sin(phi) + v.*cos(phi);
    vt = v.*sin(phi) - u.*cos(phi);
    L = vt .* r;

    ustore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = u; 
    vstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = v;
    wstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = w;
    pstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = p;
    tkestore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = tke;
    urstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = ur;
    vtstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = vt;
    Lstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = L;
    %clear u v w p tke ur vt L
end
fclose(fid);
cd(base_dir)

lev = 61; % z index to plot
tme = 1; % time index to plot

figure(1)
subplot(2,2,1)
pcolor(Xmf(:,:,lev), Ymf(:,:,lev), double(ustore(:,:,lev,tme)))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title('U (m/s)')

subplot(2,2,2)
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), squeeze(double(vstore(:,:,lev,tme))))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title('V (m/s)')

subplot(2,2,3)
tmp = squeeze(double(wstore(:,:,lev,tme)));
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title('W (m/s)')


subplot(2,2,4)
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), squeeze(double(tkestore(:,:,lev,tme))))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title('TKE (m^2/s^2)')



figure(2)

subplot(1,3,1)
tmp = squeeze(double(urstore(:,:,lev,tme)));
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp(:,:,1))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title(['U_r (m/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)

subplot(1,3,2)
tmp = squeeze(double(vtstore(:,:,lev,tme)));
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp(:,:,1))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title(['V_t (m/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)

subplot(2,2,3)
tmp = squeeze(double(wstore(:,:,lev,tme)));
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title(['W (m/s)at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)

subplot(1,3,3)
tmp = squeeze(double(Lstore(:,:,lev,tme)));
pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp(:,:,1))
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title(['L (m^2/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)

