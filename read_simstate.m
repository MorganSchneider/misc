% Some upgrades needed: Need to count number of rain drops and record their
% velocities

% clear
% close all


base_dir = '/Users/schneider/Documents/';
dir_loc = [base_dir 'sims']; % SimRadar output directory
sim_dir = uigetdir(dir_loc); % Location of IQ files
filename = blib('choosefile', sim_dir, '*.simstate');


filt_deb_ht = true; % Filters based on debris height
tol = 1; % Factor multiplying 0.5º to determine the range of elevation angles 
% where debris are recorded (e.g., 0.5º*1 collects debris within +/- 0.5º)
collect_vel_stats = true;

% Name of file input
fname = filename(1:end-9);
% Read Simstate file
sdat = simradarstate([fname '.simstate']);
% Read IQ file
simrad = simradariq([fname '.iq']);
% IMPORTANT!!! If you update SimRadar again (from git.arrc), go into simradarstate.m and
% change 63*1024 to 64*1024 at the bottom

x = sdat.pos(1,:);
y = sdat.pos(2,:);
z = sdat.pos(3,:);
us = sdat.vel(1,:);
vs = sdat.vel(2,:);
ws = sdat.vel(3,:);
ori_alp = sdat.ori(1,:);
ori_bet = sdat.ori(2,:);
ori_gam = sdat.ori(3,:);

% Find rain points (they all have the same orientation)
bet_rain = mode(ori_bet);
alp_rain = mode(ori_alp);
gam_rain = mode(ori_gam);

yn = ori_bet == bet_rain & ori_alp == alp_rain & ori_gam == gam_rain;
[rain_ind] = find(yn);
deb_ind = find(~yn);

% I'm guessing rcs is shh, shv, svh, svv, should verify
rcs_hh = sdat.rcs(1,:);
rcs_vv = sdat.rcs(2,:);
rcs_vh = sdat.rcs(3,:);
rcs_hv = sdat.rcs(4,:);

%%%%%%%%%%%%%%%%%%%%%
% Import radar data %
%%%%%%%%%%%%%%%%%%%%%
% Sector width
sec = abs(simrad.params.scan_end - simrad.params.scan_start);
np_per_deg = 1.0 / simrad.params.scan_delta;

% Number of pulses per ray
np = round(0.5 * np_per_deg);

fprintf('Using %d pulses per ray (%.1f deg)\n', np, np / np_per_deg);

% Total number of samples in the file
ns = numel(simrad.az_deg);

% Maximum number of rays resulted
nray = floor(ns / np);

% Number of pulses per sweep
pps = (simrad.params.scan_end - simrad.params.scan_start) / simrad.params.scan_delta;
if strcmp(simrad.params.scan_mode, 'PPI') == true
    nrs = floor(pps / np);
else
    nrs = floor(pps / np);
end

% Number of sweeps in the file
nsweep = floor(nray / nrs);

%% Gather the scan attributes and simrada
az_rad = deg2rad(simrad.az_deg(1 : np : np * nrs));
el_rad = deg2rad(simrad.el_deg(1 : np : np * nrs));
r = (0 : simrad.params.range_count - 1) * simrad.params.range_delta + simrad.params.range_start;
el = deg2rad(simrad.el_deg(1));

iqh = reshape(simrad.iqh(:, 1:nsweep * nrs * np), [simrad.params.range_count, np, nrs, nsweep]);
iqh = permute(iqh, [1 3 2 4]);
iqv = reshape(simrad.iqv(:, 1:nsweep * nrs * np), [simrad.params.range_count, np, nrs, nsweep]);
iqv = permute(iqv, [1 3 2 4]);

% Some moment products
sh = real(mean(iqh .* conj(iqh), 3));
sv = real(mean(iqv .* conj(iqv), 3));
v = -simrad.params.va / pi * angle(sum(iqh(:, :, 2:end, :) .* conj(iqh(:, :, 1:end-1, :)), 3));
mh = repmat(mean(iqh, 3), [1 1 np 1]);
mv = repmat(mean(iqv, 3), [1 1 np 1]);
sh_ac = mean((iqh - mh) .* conj(iqh - mh), 3);
sv_ac = mean((iqv - mv) .* conj(iqv - mv), 3);

rhohv = abs(mean((iqh - mh) .* conj(iqv - mv), 3)) ./ sqrt(sh_ac .* sv_ac);
deltadp = angle(mean((iqh - mh) .* conj(iqv - mv), 3) ./ sqrt(sh_ac .* sv_ac)) * 180 / pi;
% Signal in dB, zdr in dB
sh = 10 * log10(squeeze(sh));
sv = 10 * log10(squeeze(sv));
v = squeeze(v);
rhohv = squeeze(rhohv);
deltadp = squeeze(deltadp);
zdr = sh - sv;

% Corrections factor to normalize tx power, Gt, Gr, lambda, etc.
zcor =  -10 * log10(simrad.params.tx_power_watt) - simrad.params.antenna_gain_dbi;

% Range correction for z like
rcor = 10 * log10(r(:) .^ 2);
rcor = repmat(rcor, [1 numel(az_rad) nsweep]);

% Now we apply the range correction factor and z correction factor
zh = sh + rcor + zcor;
zv = sv + rcor + zcor;

% Grid data for plotting
[az_mat, r_mat] = meshgrid(az_rad, r);
xx = r_mat .* sin(az_mat) * cos(el);
yy = r_mat .* cos(az_mat) * cos(el);
zz = r_mat * sin(el);

% Plot debris locations in 3D (highlight those moving upward)
figure(1)
scatter3(x(deb_ind), y(deb_ind), z(deb_ind), '.r')
hold on
ind = ws(deb_ind) > 0;
scatter3(x(deb_ind(ind)), y(deb_ind(ind)), z(deb_ind(ind)), ws(deb_ind(ind)), 'b')
xlabel('X(m)', 'FontSize', 14)
ylabel('Y(m)', 'FontSize', 14)
zlabel('Z(m)', 'FontSize', 14)
hold off
view(3)
legend('Debris', 'Debris in updraft', 'Location', 'northeast')



if ~exist('WS','var')
    F = scatteredInterpolant(x(:),y(:),z(:),ws(:),'natural');
    [xq,yq,zq] = ndgrid(-430:10:430, 1500:10:2500, 0:10:850);
    WS = F(xq,yq,zq);
end

figure(8)
axis tight manual
scatter3(x(deb_ind), y(deb_ind), z(deb_ind), '.r')
hold on
ind = ws(deb_ind) > 0;
scatter3(x(deb_ind(ind)), y(deb_ind(ind)), z(deb_ind(ind)), ws(deb_ind(ind)), 'ok')
[f,v] = isosurface(xq, yq, zq, WS, 10);
patch('Vertices', v, 'Faces', f, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold off
xlabel('X(m)', 'FontSize', 14)
ylabel('Y(m)', 'FontSize', 14)
zlabel('Z(m)', 'FontSize', 14)
xlim([-400 400])
ylim([1600 2400])
zlim([10 700])
view(3)
legend('Debris', 'Debris in updraft', 'Location', 'northeast')


% Plot histogram of debris
figure(2)
[dnum,hgt] = hist(z(deb_ind), 5:10:995);
% semilogy(hgt,dnum)
% xlabel('Height (m)', 'FontSize', 14)
% ylabel('Number of debris', 'FontSize', 14)
semilogx(dnum, hgt, 'k', 'LineWidth', 1)
xlabel('Number of debris', 'FontSize', 14)
ylabel('Height (m)', 'FontSize', 14)
grid on

% Get coordinates of debris
x_deb = x(deb_ind);
y_deb = y(deb_ind);
z_deb = z(deb_ind);
rcs_hh_deb = rcs_hh(deb_ind);
% Compute mean beam height and beam heights +/- 0.5º*tol
bh_mean = nanmean(r) * sin(el);
bh_min = nanmean(r) * sin(el - 0.5*pi/180*tol);
bh_max = nanmean(r) * sin(el + 0.5*pi/180*tol);

% Filter by beam heights
if filt_deb_ht
    x_deb = x_deb(z_deb >= bh_min & z_deb <= bh_max);
    y_deb = y_deb(z_deb >= bh_min & z_deb <= bh_max);
    z_deb = z_deb(z_deb >= bh_min & z_deb <= bh_max);
    u_deb = us(z_deb >= bh_min & z_deb <= bh_max);
    v_deb = vs(z_deb >= bh_min & z_deb <= bh_max);
    w_deb = ws(z_deb >= bh_min & z_deb <= bh_max);
end

deb_cnt = zeros(size(xx));
rcs_hh_mean = zeros(size(xx));
rcs_hh_max = rcs_hh_mean;
% Count the number of debris in each gate, based on which gate is closest to the debris
% Loop through all debris
if collect_vel_stats
    for idx = 1:max(size(x_deb))
       % Compute distance between the debris and gate in xy plane
       rad_deb = sqrt((x_deb(idx)-xx).^2 + (y_deb(idx)-yy).^2);
       % Find the minimum distance and corresponding gate
       [row, col] = find(min(rad_deb(:)) == rad_deb);
       % Add the debris to the gate's count
       deb_cnt(row(1),col(1)) = deb_cnt(row(1),col(1)) + 1;
       % Collect the debris velocities by gate
       u_save{row(1),col(1)}(deb_cnt(row(1),col(1))) = u_deb(idx);
       v_save{row(1),col(1)}(deb_cnt(row(1),col(1))) = v_deb(idx);
       w_save{row(1),col(1)}(deb_cnt(row(1),col(1))) = w_deb(idx);
       % Calculate the mean RCS in the gate as well as the maximum RCS (not very
       % important for now)
       rcs_hh_mean(row(1),col(1)) = rcs_hh(idx) + rcs_hh_mean(row(1),col(1));
       rcs_hh_max(row(1),col(1)) = max(abs(rcs_hh_max(row(1),col(1))), abs(rcs_hh_deb(idx)));
    end
    rcs_hh_mean = rcs_hh_mean ./ deb_cnt;
else
    for idx = 1:max(size(x_deb))
       % Compute distance between the debris and gate in xy plane
       rad_deb = sqrt((x_deb(idx)-xx).^2 + (y_deb(idx)-yy).^2);
       % Find the minimum distance and corresponding gate
       [row,col] = find(min(rad_deb(:)) == rad_deb);
       % Add the debris to the gate's count
       deb_cnt(row(1),col(1)) = deb_cnt(row(1),col(1)) + 1;
       % Calculate the mean RCS in the gate as well as the maximum RCS (not very
       % important for now)
       rcs_hh_mean(row(1),col(1)) = rcs_hh(idx) + rcs_hh_mean(row(1),col(1));
       rcs_hh_max(row(1),col(1)) = max(abs(rcs_hh_max(row(1),col(1))), abs(rcs_hh_deb(idx)));
    end
    rcs_hh_mean = rcs_hh_mean ./ deb_cnt;
end

% The points just outside of the domain have really high numbers since the
% sector scan is not square
mask = zeros(size(zh));
mask(2:size(xx,1), 2:size(xx,2)) = 1;
deb_cnt = deb_cnt .* mask; % Remove border points


figure(3)
pcolor(xx, yy, deb_cnt)
shading flat
colorbar
xlabel('X(m)')
ylabel('Y(m)')
title('Debris Count in scan')

% Scatter plot of Zh vs debris count
figure(4)
subplot(1,2,1)
scatter(zh(:), deb_cnt(:))
xlabel('Zh')
ylabel('Debris Count')
subplot(1,2,2)
[N,c] = hist3([zh(:), deb_cnt(:)], 'Ctrs', {-7.5:2.5:57.5 0:2:24});
[cx,cy] = meshgrid(c{1}, c{2});
N(N==0) = nan;
pcolor(cx, cy, N')
colorbar
shading flat
xlim([-10 60])
ylim([0 25])
xlabel('Zh')
ylabel('Debris count')

% Scatter plot of Rhohv vs debris count
figure(5)
subplot(1,2,1)
scatter(rhohv(:), deb_cnt(:))
xlabel('Rhohv')
ylabel('Debris Count')
subplot(1,2,2)
[N,c] = hist3([rhohv(:), deb_cnt(:)], 'Ctrs', {0:0.05:1 0:2:24});
[cx,cy] = meshgrid(c{1}, c{2});
N(N==0) = nan;
pcolor(cx, cy, N')
colorbar
shading flat
xlim([0 1])
ylim([0 25])
xlabel('Rhohv')
ylabel('Debris count')

% Max RCS - Mean RCS
figure(6)
pcolor(xx, yy, rcs_hh_max - rcs_hh_mean)
figure(gcf)
shading flat
colorbar
xlabel('Rhohv')
ylabel('Max RCS - Mean RCS')

% Same as Figure 6 except only where there are 3 or more debris
tmp = rcs_hh_max - rcs_hh_mean;
ind_pos = deb_cnt > 2;
figure(7)
scatter(rhohv(ind_pos), tmp(ind_pos))
xlabel('Rhohv')
ylabel('Max RCS - Mean RCS')


