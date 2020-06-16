% A convenient script to browse and proecess an iq data file from the
% simulator
%
% Boon Leng Cheong
% Advanced Radar Research Center
% University of Oklahoma
% 2/9/2016
%
clear;
close all;

base_dir = '/Users/schneider/Documents/';
% fig_dir = [base_dir 'imgs'];

deb_titl = input('Enter folder name:','s');
dir_loc = [base_dir deb_titl];
% addpath('/Users/schneider/Documents/sims/');

if ~exist('skipload', 'var')
    if exist('~/Downloads/density', 'dir')
        % filename = boonlib('choosefile', '~/Downloads/density', '*.iq');
        filename = boonlib('choosefile', dirlist, '*.iq');
    else
        filename = boonlib('choosefile', dir_loc, '*.iq');
    end
    if ~isempty(filename)
        dat = simradariq(filename);
    else
        return;
    end
%     statename = [filename(1: max(size(filename)) - 3) '.simstate'];
%     cd(dir_loc)
%         if exist(statename, 'file')
%             deb = simradarstate(filename);
%         end
%     cd(base_dir)
end
 
smooth_plots = true;
load_LES = true; 
single_time = false;
if single_time
    st_idx = input('Time step: ');
    smooth_plots = false;
end
dx = 20;
g = 9.8;
wplot_int_pos = [0 10 20 30 40];
wplot_int_neg = [-20 -10];
rho_clims = [0 1];
zdr_clims = [-5 5];


%Set LES directory
% if load_LES
LES_names = {'suctvort', 'suctvort_large', 'moore', 'trans'};
fprintf(['Choose LES simulation','\n'])
for idx = 1:max(size(LES_names))
    fprintf([num2str(idx) ': ' LES_names{idx} '\n'])
end
LES_sim_cho = input('Choose LES simulation: ');
if LES_sim_cho == 1
    LES_dir = '/Users/schneider/Documents/tables/suctvort/';
    v0 = 100;
    sim_name = 'Suctvort';
elseif LES_sim_cho == 2
    LES_dir = '/Users/schneider/Documents/tables/suctvort_large/';
    v0 = 100;
    sim_name = 'Suctvort Large';
elseif LES_sim_cho == 3
    LES_dir = '/Users/schneider/Documents/tables/moore/';
    v0 = 100;
    sim_name = 'Moore';
    dx = 40;
    wplot_int_pos = [0 10 20 30 40];
    wplot_int_neg = [-10 -10];
elseif LES_sim_cho == 4
    LES_dir = '/Users/schneider/Documents/tables/trans/';
    v0 = 100;
    sim_name = 'Translating';
    dx = 20; 
end
% end


class_save = true;
mom_dir = [base_dir '/moms/' deb_titl '/'];
if ~exist(mom_dir, 'dir')
    mkdir(mom_dir);
end
fig_dir = [base_dir '/imgs/' deb_titl '/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

if class_save
    deb_flag = input('With debris (1: Yes, 0: No');
    if deb_flag
        deb_title = 'D';
    else
        deb_title = 'ND';
    end
end

% Sector width
sec = abs(dat.params.scan_end - dat.params.scan_start);
np_per_deg = 1.0 / dat.params.scan_delta;

% Number of pulses per ray
np = round(0.5 * np_per_deg);

fprintf('Using %d pulses per ray (%.1f deg)\n', np, np / np_per_deg);

% Total number of samples in the file
ns = numel(dat.az_deg);

% Maximum number of rays resulted
nray = floor(ns / np);

% Number of pulses per sweep
pps = (dat.params.scan_end - dat.params.scan_start) / dat.params.scan_delta;
if strcmp(dat.params.scan_mode, 'PPI') == true
    nrs = floor(pps / np);
else
    nrs = floor(pps / np);
end

% Number of sweeps in the file
nsweep = floor(nray / nrs);
sweep_time = zeros(1, nsweep);
for ldx = 1:nsweep
    sweep_time(ldx) = nanmean(dat.scan_time((ldx-1) * round(pps) + (1:pps))); 
end

%% Gather the scan attributes and data
az_rad = deg2rad(dat.az_deg(1 : np : np * nrs));
el_rad = deg2rad(dat.el_deg(1 : np : np * nrs));
r = (0 : dat.params.range_count - 1) * dat.params.range_delta + dat.params.range_start;
el = deg2rad(dat.el_deg(1));

iqh = reshape(dat.iqh(:, 1:nsweep * nrs * np), [dat.params.range_count, np, nrs, nsweep]);
iqh = permute(iqh, [1 3 2 4]);
iqv = reshape(dat.iqv(:, 1:nsweep * nrs * np), [dat.params.range_count, np, nrs, nsweep]);
iqv = permute(iqv, [1 3 2 4]);

% Some moment products
sh = real(mean(iqh .* conj(iqh), 3));
sv = real(mean(iqv .* conj(iqv), 3));
vh = -dat.params.va / pi * angle(sum(iqh(:, :, 2:end, :) .* conj(iqh(:, :, 1:end-1, :)), 3));
vv = -dat.params.va / pi * angle(sum(iqv(:, :, 2:end, :) .* conj(iqv(:, :, 1:end-1, :)), 3));
mh = repmat(mean(iqh, 3), [1 1 np 1]);
mv = repmat(mean(iqv, 3), [1 1 np 1]);
sh_ac = mean((iqh - mh) .* conj(iqh - mh), 3);
sv_ac = mean((iqv - mv) .* conj(iqv - mv), 3);

% rhohv = abs(mean(iqh .* conj(iqv), 3)) ./ sqrt(sh_ac .* sv_ac);  % Assume signals are zero-mean
rhohv = abs(mean((iqh - mh) .* conj(iqv - mv), 3)) ./ sqrt(sh_ac .* sv_ac);

% Signal in dB, zdr in dB
sh = 10 * log10(squeeze(sh));
sv = 10 * log10(squeeze(sv));
vh = squeeze(vh);
rhohv = squeeze(rhohv);
zdr = sh - sv;

% Corrections factor to normalize tx power, Gt, Gr, lambda, etc.
% zcor =  -10 * log10(dat.params.tx_power_watt) - dat.params.antenna_gain_dbi - 10 * log10(dat.params.body_per_cell / 10000) + 50;

% zcor =  -10 * log10(dat.params.tx_power_watt) - dat.params.antenna_gain_dbi + 50; % paper draft
% zcor =  -10 * log10(dat.params.tx_power_watt) - dat.params.antenna_gain_dbi + 54.48;
zcor =  -10 * log10(dat.params.tx_power_watt) - dat.params.antenna_gain_dbi;

% Range correction for z like
% rcor = 10 * log10((1.0e-3 * r(:)) .^ 2);  % paper draft
rcor = 10 * log10(r(:) .^ 2);
rcor = repmat(rcor, [1 numel(az_rad) nsweep]);

% Now we apply the range correction factor and z correction factor
zh = sh + rcor + zcor;
zv = sv + rcor + zcor;


%% Moment plots

% zdr_ind = boonlib('nwsd2ind', zdr);
% rhohv_ind = boonlib('nwsr2ind', rhohv);

figure(1)
if size(zh, 3) > 1
    zhs = zeros(size(zh,1), size(zh,2), size(zh,3));
    vs = zeros(size(vh,1), size(vh,2), size(vh,3));
    rhos = zeros(size(rhohv,1), size(rhohv,2), size(rhohv,3));
    zdrs = zeros(size(zdr,1), size(zdr,2), size(zdr,3));
    for idx = 1:size(zh, 3)
        zhs(:,:,idx) = filter2(ones(3,3) / 9, zh(:,:,idx), 'same');
        vs(:,:,idx) = filter2(ones(3,3) / 9, vh(:,:,idx), 'same');
        rhos(:,:,idx) = filter2(ones(3,3) / 9, rhohv(:,:,idx), 'same');
        zdrs(:,:,idx) = filter2(ones(3,3) / 9, zdr(:,:,idx), 'same');
    end
else
    zhs = filter2(ones(3,3) / 9, zh, 'same');
    vs = filter2(ones(3,3) / 9, vh, 'same');
    rhos = filter2(ones(3,3) / 9, rhohv, 'same');
    zdrs = filter2(ones(3,3) / 9, zdr, 'same');
end

clf

%Get first time from LES
if load_LES
     %Load LES model grid
    % [Xm, Ym, Zm] = read_LES_grid(LES_dir, base_dir, v0);
    read_LES;
    %Figure out LES file times needed based on scan time
    LES_int = 2e-4 * 625 * v0 / g;
    if nsweep == 1
        first_time = sweep_time;
    else
        first_time = sweep_time(1);
    end
    fnum = floor(nanmin(first_time / LES_int) / 10) + 1;
    fnum_orig = -1;
    if ~exist('st_idx', 'var')
        st_idx = 1;
        end_idx = size(zh, 3);
    end
end

for idx = st_idx:end_idx
    if load_LES
        fnum = floor(nanmin(sweep_time(idx) / LES_int) / 10) + 1;
        fprintf(['fnum is: ' num2str(fnum) '\n']);
        if fnum ~= fnum_orig
%             umean = zeros(size(Xm,1), size(Xm,2), size(Xm,3));
%             vmean = zeros(size(Xm,1), size(Xm,2), size(Xm,3));
%             wmean = zeros(size(Xm,1), size(Xm,2), size(Xm,3));
%             tkemean = zeros(size(Xm,1), size(Xm,2), size(Xm,3));
            %Load LES model data for the selected times
            %[umean, vmean, wmean, tkemean, LES_time] = read_LES(Xm, LES_dir, base_dir, v0, fnum);
            umean = ustore;
            vmean = vstore;
            wmean = wstore;
            tkemean = tkestore;
            LES_time = time;
%                 umean(:, :, :, (tidx-1) * 10 + [1:10]) = umn;
%                 vmean(:, :, :, (tidx-1) * 10 + [1:10]) = vmn;
%                 wmean(:, :, :, (tidx-1) * 10 + [1:10]) = wmn;
%                 tkemean(:, :, :, (tidx-1) * 10 + [1:10]) = tkemn;
%                 LES_time((tidx-1) * 10 + [1:10]) = LES_tme;

            %Interpolate for plotting
            beam_height = nanmean(r) * sin(el);
            [c, zlev] = min(abs(squeeze(Zm(1,1,:)) - beam_height));
        %     if strcmp(deb_titl, 'trans_z20')
                xctr = 0;
                yctr = nanmean(r);
        %     else
        %         xctr = 0;
        %         yctr = 0;
        %     end
            if ~exist('xx2', 'var')
                xx2 = ceil(nanmin(Xmf(:)) * dx) / dx : dx : floor(nanmax(Xmf(:)) * dx) / dx;
                yy2 = xx2;
                [XX2, YY2] = meshgrid(xx2, yy2);
                xtmp = double(Xmf(:,:,zlev));
                ytmp = double(Ymf(:,:,zlev));
                grid_set = false;
            end

            clear utmp vtmp wtmp tketmp;
            
            uint = zeros(size(XX2,1), size(XX2,2), size(umean,4));
            vint = zeros(size(XX2,1), size(XX2,2), size(umean,4));
            wint = zeros(size(XX2,1), size(XX2,2), size(umean,4));
            tkeint = zeros(size(XX2,1), size(XX2,2), size(umean,4));
            for ldx = 1:size(umean,4)
                utmp = double(umean(:,:,zlev,ldx));
                vtmp = double(vmean(:,:,zlev,ldx));
                wtmp = double(wmean(:,:,zlev,ldx));
                tketmp = double(tkemean(:,:,zlev,ldx));
                F.u = scatteredInterpolant(xtmp(:), ytmp(:), utmp(:), 'linear');
                F.v = scatteredInterpolant(xtmp(:), ytmp(:), vtmp(:), 'linear');
                F.w = scatteredInterpolant(xtmp(:), ytmp(:), wtmp(:), 'linear');
                F.tke = scatteredInterpolant(xtmp(:), ytmp(:), tketmp(:), 'linear');
                uint(:, :, ldx) = F.u(XX2, YY2); %(1:size(XX2,1), 1:size(XX2,2), ldx)
                vint(:, :, ldx) = F.v(XX2, YY2);
                wint(:, :, ldx) = F.w(XX2, YY2);
                tkeint(:, :, ldx) = F.tke(XX2, YY2);
                clear F;
            end
            if ~exist('XX2p', 'var')
                XX2p = XX2;
                YY2p = YY2 + yctr;
                grid_set = true;
            end
        end
        %Get LES_ind;
        [c, LES_ind] = min(abs(LES_time - sweep_time(idx))); % unused: c
        fnum_orig = fnum;
        fprintf(['LES index is: ' num2str(LES_ind) '\n']);
    end
    if strcmp(dat.params.scan_mode, 'PPI')
        % Pad one more extra at the end for pcolor()
        % ap = cat(2, az_rad, az_rad(end) + (az_rad(2) - az_rad(1)));

        % rp = cat(2, r, r(end) + (r(2) - r(1)));
        % sp = s;
        % sp = cat(1, sp, sp(end, :));
        % sp = cat(2, sp, sp(:, end));
        clf;
        [az_mat, r_mat] = meshgrid(az_rad, r);
        % [az_mat, r_mat] = meshgrid(ap, rp);

        xx = r_mat .* sin(az_mat) * cos(el);
        yy = r_mat .* cos(az_mat) * cos(el);
        zz = r_mat * sin(el);

        if single_time
            figure(1)
            feval('boonlib', 'bsizewin', 1, [700 700])
        else
            ha = subplot(2, 2, 1);
        end
        hs = pcolor(xx, yy, zh(:, :, idx));
        set(gca, 'DataAspect', [1 1 1])
        caxis([0 80])
        colormap(ha, boonlib('zmap'))
        shading flat
        colorbar
        set(gca, 'YDir', 'Normal')
        title('Z - Reflectivity (dBZ)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2]=contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        if single_time
            figure(2)
            feval('boonlib', 'bsizewin', 2, [700 700])
        else
            ha(2) = subplot(2, 2, 2);
        end
        hs(2) = pcolor(xx, yy, vh(:, :, idx));
        caxis([-1 1] * dat.params.va)
        colormap(ha(2), boonlib('carbmap'))
        shading flat
        colorbar
        title('V - Velocity (m/s)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        if single_time
            figure(3)
            feval('boonlib', 'bsizewin', 3, [700 700])
        else
            ha(3) = subplot(2, 2, 3);
        end
        hs(3) = pcolor(xx, yy, zdr(:, :, idx));
    %     hs(3) = pcolor(xx, yy, zdr_ind(:, :, 1));
    %     colormap(ha(3), boonlib('nwsdmap'))
        colormap(ha(3), boonlib('zmap'))
        colorbar
        caxis(zdr_clims)
        shading flat
    %     boonlib('nwsdbar');
        title('D - Differential Reflectivity (dB)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        if single_time
            figure(4)
            feval('boonlib', 'bsizewin', 4, [700 700])
        else
            ha(4) = subplot(2, 2, 4);
        end
    %     hs(4) = pcolor(xx, yy, rhohv_ind(:, :, 1));
        hs(4) = pcolor(xx, yy, real(rhohv(:,:,idx)));
        % colormap(ha(4), boonlib('rhomap'))
        % colorbar
        caxis(rho_clims)
        shading flat
        colormap(ha(4), blib('nwsrmap'))
        blib('nwsrbar');
        title('R - RhoHV')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2,h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        set(ha, 'DataAspect', [1 1 1])

    elseif strcmp(dat.params.scan_mode, 'RHI')

        [el_mat, r_mat] = meshgrid(el_rad, r);

        xx = r_mat .* cos(el_mat);
        zz = r_mat .* sin(el_mat);

        ha = subplot(2, 2, 1);
        hs = pcolor(xx, zz, zh(:, :, idx));
        caxis([0 80])
        colormap(ha, boonlib('zmap'))
        colorbar
        shading flat
        set(gca, 'YDir', 'Normal')
        title('Z - Reflectivity (dBZ)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        ha(2) = subplot(2, 2, 2);
        hs(2) = pcolor(xx, zz, vh(:, :, idx));
        caxis([-1 1] * dat.params.va)
        colormap(ha(2), boonlib('carbmap'))
        colorbar
        shading flat
        title('V - Velocity (m/s)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end


        ha(3) = subplot(2, 2, 3);
    %     hs(3) = pcolor(xx, zz, zdr_ind(:, :, 1));
        hs(3) = pcolor(xx, zz, zdr(:, :, idx));
        set(gca, 'DataAspect', [1 1 1])
    %     colormap(ha(3), boonlib('nwsdmap'))
    %     boonlib('nwsdbar');
        colormap(ha(3), boonlib('zmap'))
        colorbar
        shading flat
        caxis(zdr_clims)
        title('D - Differential Reflectivity (dB)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2,h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        rhohv(zh < 10) = NaN;
        ha(4) = subplot(2, 2, 4);
    %     hs(4) = pcolor(xx, zz, rhohv_ind(:, :, 1));
        hs(4) = pcolor(xx, zz, real(rhohv(:,:,idx)));
        % colormap(ha(4), boonlib('rhomap'))
        % colorbar
        caxis(rho_clims)
        shading flat
        colormap(ha(4), blib('nwsrmap'))
        blib('nwsrbar');
        title('R - RhoHV')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k'); % unused: c
            [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
            hold off
        end

        set(ha, 'DataAspect', [1 1 1])
    else

        fprintf('Not a scan mode I know how to plot.\n');

    end

    if single_time
        var_name = {'Z', 'v_{r}', 'Z_{DR}', '\rho_{HV}'};
        for plotdx = 1:4
            axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01])
        %         title_str = filename(max(size(dir_loc)) + 1 : max(size(filename)));
    %         if strcmp(title_str(1), '/')
    %             title_str = title_str(2:max(size(title_str)));
    %         end
            if el * 180 / pi < 10
                title_str = ['t=' num2str(roundn(sweep_time(idx),-1)) 's-E0' num2str(roundn(el*180/pi,-1))];
            else
                title_str = ['t=' num2str(roundn(sweep_time(idx),-1)) 's-E' num2str(roundn(el*180/pi,-1))];
            end
            tstr = sprintf('%s', title_str);
%             title_str=filename(max(size(dir_loc))+1:max(size(filename)));
%             if strcmp(title_str(1), '/')
%                 title_str = title_str(2:max(size(title_str)));
%             end
%             if dat.debris_counts(3) > 0
%                 tstr = sprintf(tstr, dat.params.body_per_cell, dat.debris_counts(1:3));
%                 fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2)) '_' num2str(dat.debris_counts(3))];
%             elseif dat.debris_counts(2) > 0
%                 tstr = sprintf('%s (D = %.1f, w = %d, d = [%d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:2));
%                 fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2))];
%             else
%                 tstr = sprintf('%s (D = %.1f, w = %d, no debris)', tstr, dat.params.body_per_cell, dat.debris_counts(1));
%                 fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1))];
%             end
            title(tstr, 'FontSize', 14);
            axis off

            boonlib('bsizewin', gcf, [1400 700])
            set(gcf, 'PaperPositionMode', 'Auto')
            cd(fig_dir)
                print(['singleplot_' num2str(plotdx) '_' num2str(idx) '.fig']);
            cd(base_dir)
        end
    else
        axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01])
    %         title_str=filename(max(size(dir_loc))+1:max(size(filename)));
%         if(strcmp(title_str(1),'/'))
%             title_str=title_str(2:max(size(title_str)));
%         end
        if el * 180 / pi < 10
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E0' num2str(roundn(el*180/pi,-1))];
        else
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E' num2str(roundn(el*180/pi,-1))];
        end
        tstr = sprintf('%s', title_str);
        title_str = filename(max(size(dir_loc))+1:max(size(filename)));
        if strcmp(title_str(1), '/')
            title_str = title_str(2:max(size(title_str)));
        end
        if dat.debris_counts(3) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d %d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:3));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2)) '_' num2str(dat.debris_counts(3))];
        elseif dat.debris_counts(2) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:2));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2))];
        else
            tstr = sprintf('%s (D = %.1f, w = %d, no debris)', tstr, dat.params.body_per_cell, dat.debris_counts(1));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1))];
        end
        title(tstr, 'FontSize', 14);
        axis off

        boonlib('bsizewin', gcf, [1400 700])
        set(gcf, 'PaperPositionMode', 'Auto')
        cd(fig_dir)
            print([fig_name '_' num2str(idx) '.png'],'-dpng','-r0')
        cd(base_dir)
    end

    %Smoothed figure
    if smooth_plots
        if strcmp(dat.params.scan_mode, 'PPI')
            % Pad one more extra at the end for pcolor()
            % ap = cat(2, az_rad, az_rad(end) + (az_rad(2) - az_rad(1)));
            % rp = cat(2, r, r(end) + (r(2) - r(1)));
            % sp = s;
            % sp = cat(1, sp, sp(end, :));
            % sp = cat(2, sp, sp(:, end));
%             if load_LES
%             [c,LES_ind]=min(abs(LES_time-sweep_time(idx)));
%             fprintf(['LES index is: ' num2str(LES_ind) '\n']);
%             end
            figure(2)
            clf
            [az_mat, r_mat] = meshgrid(az_rad, r);
            % [az_mat, r_mat] = meshgrid(ap, rp);

            xx = r_mat .* sin(az_mat) * cos(el);
            yy = r_mat .* cos(az_mat) * cos(el);
            zz = r_mat * sin(el);

            ha = subplot(2, 2, 1);
            hs = pcolor(xx, yy, zhs(:, :, idx));
            set(gca, 'DataAspect', [1 1 1])
            caxis([0 80])
            colormap(ha, boonlib('zmap'))
            shading flat
            colorbar
            set(gca, 'YDir', 'Normal')
            title('Z - Reflectivity (dBZ)')
            if(load_LES)
                hold on
                [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k');
                [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
                hold off
            end

            ha(2) = subplot(2, 2, 2);
            hs(2) = pcolor(xx, yy, vs(:, :, idx));
            caxis([-1 1] * dat.params.va)
            colormap(ha(2), boonlib('carbmap'))
            shading flat
            colorbar
            title('V - Velocity (m/s)')
            if load_LES
                hold on
                [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k');
                [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
                hold off
            end

            ha(3) = subplot(2, 2, 3);
            hs(3) = pcolor(xx, yy, zdrs(:, :, idx));
        %     hs(3) = pcolor(xx, yy, zdr_ind(:, :, 1));
        %     colormap(ha(3), boonlib('nwsdmap'))
            colormap(ha(3), boonlib('zmap'))
            colorbar
            caxis(zdr_clims)
            shading flat
        %     boonlib('nwsdbar');
            title('D - Differential Reflectivity (dB)')
            if load_LES
                hold on
                [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k');
                [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
                hold off
            end

            ha(4) = subplot(2, 2, 4);
        %     hs(4) = pcolor(xx, yy, rhohv_ind(:, :, 1));
            hs(4) = pcolor(xx, yy, real(rhos(:,:,idx)));
            % colormap(ha(4), boonlib('rhomap'))
            % colorbar
            caxis(rho_clims)
            shading flat
            colormap(ha(4), blib('nwsrmap'))
            blib('nwsrbar');
            title('R - RhoHV')
            if load_LES
                hold on
                [c, h] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_pos, '-k');
                [c2, h2] = contour(XX2p, YY2p, wint(:,:,LES_ind), wplot_int_neg, '--k');
                hold off
            end

            set(ha, 'DataAspect', [1 1 1])
        end

        axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01])
    %         title_str=filename(max(size(dir_loc))+1:max(size(filename)));
    %         if(strcmp(title_str(1),'/'))
    %             title_str=title_str(2:max(size(title_str)));
    %         end
        if el * 180 / pi < 10
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E0' num2str(roundn(el*180/pi,-1))];
        else
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E' num2str(roundn(el*180/pi,-1))];
        end
        tstr = sprintf('%s', title_str);
        title_str = filename(max(size(dir_loc))+1:max(size(filename)));
        if strcmp(title_str(1), '/')
             title_str = title_str(2:max(size(title_str)));
        end
        if dat.debris_counts(3) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d %d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:3));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2)) '_' num2str(dat.debris_counts(3))];
        elseif dat.debris_counts(2) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:2));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2))];
        else
            tstr = sprintf('%s (D = %.1f, w = %d, no debris)', tstr, dat.params.body_per_cell, dat.debris_counts(1));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1))];
        end
        title(tstr, 'FontSize', 14);
        axis off

        boonlib('bsizewin', gcf, [1400 700])
        set(gcf, 'PaperPositionMode', 'Auto')
        cd(fig_dir)
            print(['smooth_' fig_name '_' num2str(idx) '.png'],'-dpng')
        cd(base_dir)
    end
end

%% 
% if size(sh, 3) > 1
%     nloop = 1;
%     for k = 1:nloop
%         for ii = 1:size(sh, 3)
%             set(hs(1), 'CData', zh(:, :, ii))
%             set(hs(2), 'CData', v(:, :, ii))
% %             set(hs(3), 'CData', zdr_ind(:, :, ii))
% %             set(hs(4), 'CData', rhohv_ind(:, :, ii))
%             set(hs(3), 'CData', zdr(:, :, ii))
%             set(hs(4), 'CData', rhohv(:, :, ii))
%             % pause(0.15)
%             cd(fig_dir)
%             print([fig_name '_' num2str(ii) '.png'],'-dpng', '-r0');
%             cd(base_dir)
%         end
%         pause(0.75)
%     end
% end

if class_save
    param = dat.params;
    iqh = dat.iqh;
    iqv = dat.iqv;
    az_deg_iq = dat.az_deg;
    el_deg_iq = dat.el_deg;
    scan_time_iq = dat.scan_time;
    moms.zh = zh;
    moms.zv = zv;
    moms.zdr = zdr;
    moms.rhohv = rhohv;
    moms.vh = vh;
    moms.vv = vv;
    moms.az_rad = az_rad;
    moms.el_rad = el_rad;
    moms.range = r_mat;
    moms.xx = xx;
    moms.yy = yy;
    moms.zz = zz;
    cd(mom_dir)
    save(['rs_' num2str(roundn(el*180/pi,-1)) '.mat'], 'param', 'iqh', 'iqv', 'az_deg_iq', 'el_deg_iq',...
        'scan_time_iq', 'moms');
    cd(base_dir)
end
%%
% figure(2)
% clf
% subplot(3, 1, 1)
% cplot(squeeze(iqh(5, 12, :, 1)))
% subplot(3, 1, 2)
% cplot(squeeze(iqv(5, 12, :, 1)))
% subplot(3, 1, 3)
% cplot(squeeze(iqh(5, 12, :, 1)) - squeeze(iqv(5, 12, :, 1)))

% for idx=1:size(zh,3)
%     if strcmp(dat.params.scan_mode, 'PPI')
%         % Pad one more extra at the end for pcolor()
%         % ap = cat(2, az_rad, az_rad(end) + (az_rad(2) - az_rad(1)));
%         % rp = cat(2, r, r(end) + (r(2) - r(1)));
%         % sp = s;
%         % sp = cat(1, sp, sp(end, :));
%         % sp = cat(2, sp, sp(:, end));
%         if(load_LES)
%         [c,LES_ind]=min(abs(LES_time-sweep_time(idx)));
%         fprintf(['LES index is: ' num2str(LES_ind) '\n']);
%         end
%         figure(2)
%         clf
%         [az_mat, r_mat] = meshgrid(az_rad, r);
%         % [az_mat, r_mat] = meshgrid(ap, rp);
% 
%         xx = r_mat .* sin(az_mat) * cos(el);
%         yy = r_mat .* cos(az_mat) * cos(el);
%         zz = r_mat * sin(el);
% 
%         ha = subplot(2, 2, 1);
%         hs = pcolor(xx, yy, zhs(:, :, idx));
%         set(gca, 'DataAspect', [1 1 1])
%         caxis([0 80])
%         colormap(ha, boonlib('zmap'))
%         shading flat
%         colorbar
%         set(gca, 'YDir', 'Normal')
%         title('Z - Reflectivity (dBZ)')
%         if(load_LES)
%             hold on
%             [c,h]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_pos,'-k');
%             [c2,h2]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_neg,'--k');
%             hold off
%         end
% 
%         ha(2) = subplot(2, 2, 2);
%         hs(2) = pcolor(xx, yy, vs(:, :, idx));
%         caxis([-1 1] * dat.params.va)
%         colormap(ha(2), boonlib('carbmap'))
%         shading flat
%         colorbar
%         title('V - Velocity (m/s)')
%         if(load_LES)
%             hold on
%             [c,h]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_pos,'-k');
%             [c2,h2]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_neg,'--k');
%             hold off
%         end
% 
%         ha(3) = subplot(2, 2, 3);
%         hs(3) = pcolor(xx, yy, zdrs(:, :, idx));
%     %     hs(3) = pcolor(xx, yy, zdr_ind(:, :, 1));
%     %     colormap(ha(3), boonlib('nwsdmap'))
%         colormap(ha(3),boonlib('zmap'))
%         colorbar
%         caxis(zdr_clims)
%         shading flat
%     %     boonlib('nwsdbar');
%         title('D - Differential Reflectivity (dB)')
%         if(load_LES)
%             hold on
%             [c,h]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_pos,'-k');
%             [c2,h2]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_neg,'--k');
%             hold off
%         end
% 
%         ha(4) = subplot(2, 2, 4);
%     %     hs(4) = pcolor(xx, yy, rhohv_ind(:, :, 1));
%         hs(4) = pcolor(xx, yy, real(rhos(:,:,idx)));
%         colormap(ha(4),boonlib('rhomap'))
%         colorbar
%         caxis(rho_clims)
%         shading flat
%     %     colormap(ha(4), boonlib('nwsrmap'))
%     %     boonlib('nwsrbar');
%         title('R - RhoHV')
%         if(load_LES)
%             hold on
%             [c,h]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_pos,'-k');
%             [c2,h2]=contour(XX2p,YY2p,wint(:,:,LES_ind),wplot_int_neg,'--k');
%             hold off
%         end
% 
%         set(ha, 'DataAspect', [1 1 1])
%     end
% 
%     axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01])
% %         title_str=filename(max(size(dir_loc))+1:max(size(filename)));
% %         if(strcmp(title_str(1),'/'))
% %             title_str=title_str(2:max(size(title_str)));
% %         end
%     if(el*180/pi<10)
%         title_str=[sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E0' num2str(roundn(el*180/pi,-1))];
%     else
%         title_str=[sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E' num2str(roundn(el*180/pi,-1))];
%     end
%     tstr = sprintf('%s', title_str);
%     title_str=filename(max(size(dir_loc))+1:max(size(filename)));
%     if(strcmp(title_str(1),'/'))
%          title_str=title_str(2:max(size(title_str)));
%     end
%     if dat.debris_counts(3) > 0
%         tstr = sprintf('%s (D = %.1f, w = %d, d = [%d %d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:3));
%         fig_name=[title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2)) '_' num2str(dat.debris_counts(3))];
%     elseif dat.debris_counts(2) > 0
%         tstr = sprintf('%s (D = %.1f, w = %d, d = [%d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:2));
%         fig_name=[title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2))];
%     else
%         tstr = sprintf('%s (D = %.1f, w = %d, no debris)', tstr, dat.params.body_per_cell, dat.debris_counts(1));
%         fig_name=[title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1))];
%     end
%     title(tstr, 'FontSize', 14);
%     axis off
% 
%     boonlib('bsizewin', gcf, [1400 700])
%     set(gcf,'PaperPositionMode','Auto')
%     cd(fig_dir)
%         print(['smooth_' fig_name '_' num2str(idx) '.png'],'-dpng')
%     cd(base_dir)
% end

figure(3)
clf;
if load_LES
    [c, LES_ind_st] = nanmin(abs(nanmin(sweep_time) - LES_time));
    [c, LES_ind_fn] = nanmin(abs(nanmax(sweep_time) - LES_time));
end
if size(zh,3) > 1
    if strcmp(dat.params.scan_mode, 'PPI')
        % Pad one more extra at the end for pcolor()
        % ap = cat(2, az_rad, az_rad(end) + (az_rad(2) - az_rad(1)));
        % rp = cat(2, r, r(end) + (r(2) - r(1)));
        % sp = s;
        % sp = cat(1, sp, sp(end, :));
        % sp = cat(2, sp, sp(:, end));

        figure(3)
        [az_mat, r_mat] = meshgrid(az_rad, r);
        % [az_mat, r_mat] = meshgrid(ap, rp);

        xx = r_mat .* sin(az_mat) * cos(el);
        yy = r_mat .* cos(az_mat) * cos(el);
        zz = r_mat * sin(el);

        ha = subplot(2, 2, 1);
        hs = pcolor(xx, yy, nanmean(zh,3));
        set(gca, 'DataAspect', [1 1 1])
        caxis([0 80])
        colormap(ha, boonlib('zmap'))
        shading flat
        colorbar
        set(gca, 'YDir', 'Normal')
        title('Z - Reflectivity (dBZ)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_pos, '-k');
            [c2, h2] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_neg, '--k');
            hold off
        end

        ha(2) = subplot(2, 2, 2);
        hs(2) = pcolor(xx, yy,nanmean(vs,3));
        caxis([-1 1] * dat.params.va)
        colormap(ha(2), boonlib('carbmap'))
        shading flat
        colorbar
        title('V - Velocity (m/s)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_pos, '-k');
            [c2, h2] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_neg, '--k');
            hold off
        end

        ha(3) = subplot(2, 2, 3);
        hs(3) = pcolor(xx, yy, nanmean(zdr,3));
    %     hs(3) = pcolor(xx, yy, zdr_ind(:, :, 1));
    %     colormap(ha(3), boonlib('nwsdmap'))
        colormap(ha(3), boonlib('zmap'))
        colorbar
        caxis(zdr_clims)
        shading flat
    %     boonlib('nwsdbar');
        title('D - Differential Reflectivity (dB)')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_pos, '-k');
            [c2, h2] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_neg, '--k');
            hold off
        end

        ha(4) = subplot(2, 2, 4);
    %     hs(4) = pcolor(xx, yy, rhohv_ind(:, :, 1));
        hs(4) = pcolor(xx, yy, real(nanmean(rhohv,3)));
        % colormap(ha(4), boonlib('rhomap'))
        % colorbar
        caxis(rho_clims)
        shading flat
        colormap(ha(4), blib('nwsrmap'))
        blib('nwsrbar');
        title('R - RhoHV')
        if load_LES
            hold on
            [c, h] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_pos, '-k');
            [c2, h2] = contour(XX2p, YY2p, nanmean(wint(:,:,LES_ind_st:LES_ind_fn), 3), wplot_int_neg, '--k');
            hold off
        end

        set(ha, 'DataAspect', [1 1 1])

        axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01])
%         title_str=filename(max(size(dir_loc))+1:max(size(filename)));
%         if(strcmp(title_str(1),'/'))
%             title_str=title_str(2:max(size(title_str)));
%         end
        if el * 180 / pi < 10
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E0' num2str(roundn(el*180/pi,-1))];
        else
            title_str = [sim_name ' t' num2str(roundn(sweep_time(idx),-1)) 's-E' num2str(roundn(el*180/pi,-1))];
        end
        tstr = sprintf('%s', title_str);
        title_str = filename(max(size(dir_loc))+1:max(size(filename)));
        if strcmp(title_str(1), '/')
            title_str = title_str(2:max(size(title_str)));
        end
        if dat.debris_counts(3) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d %d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:3));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2)) '_' num2str(dat.debris_counts(3))];
        elseif dat.debris_counts(2) > 0
            tstr = sprintf('%s (D = %.1f, w = %d, d = [%d])', tstr, dat.params.body_per_cell, dat.debris_counts(1:2));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1)) '_' num2str(dat.debris_counts(2))];
        else
            tstr = sprintf('%s (D = %.1f, w = %d, no debris)', tstr, dat.params.body_per_cell, dat.debris_counts(1));
            fig_name = [title_str '_' num2str(dat.params.body_per_cell) '_' num2str(dat.debris_counts(1))];
        end
        title(tstr, 'FontSize', 14);
        axis off

        boonlib('bsizewin', gcf, [1400 700])
        set(gcf, 'PaperPositionMode', 'Auto')
        cd(fig_dir)
            print(['smooth_' fig_name '_mean.png'], '-dpng')
        cd(base_dir)
    end
end



