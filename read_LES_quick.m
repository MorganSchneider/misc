% This script reads in the 4-dimensional LES data (x,y,z,t) from the binary
% files and plots horizontal slices at the desired height. Axisymmetric
% velocity statistics are also computed when axi_stats_flag=true.

% Options to experiment with: 
% 1) Change the height (variable: lev) and time (variable fnum and tme) to
% see how the simulated tornado evolves in height and time
% 2) Use the program to save figs (save_figs=true)
% 3) Compute the mean 3D winds over the entire simulation (set
% compute_mean=true and fnum_st=1 and fnum_end=16);
% 4) Compute the axisymmetric wind field (change axi_stats_flag to true)

clear
close all

% Name of simulation folder within the primary data directory
%drive_name = 'seagate1';
% OPTIONS: suctvort  suctvort_large  torgen  moore  onecell  breakdown
sim_base = 'twocell';

% Rename to your matlab directory
base_dir = '/Users/schneider/Documents';
save_dir = [base_dir '/sims/les/' sim_base];
% Location of the LES data 
data_dir = [base_dir '/tables/les/' sim_base '/'];
% Save output from axisymmetric calculation

%%%%%%%%%% You can modify this structure to the appropriate directory
%%%%%%%%%% or simplify it. I use it this way because I switch drives often 
%%%%%%%%%% have multiple subdirectories of experiments

if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end


vtrans = 0;
v0 = 100; % This is the characteristic velocity, which scales the non-dimensional model results.  
% In general, we leave it at 100 m/s
g = 9.8;
num_times = 10; % Number of times stored in the file (does not change)
name_type = 1; % Different OSs or environments can produce different file conventions, so you may have to change this to 1, 2, or 3
view_fact = 0.9; % fraction of the domain in view (smaller values zoom in on the tornado)

axi_stats_flag = true; % Compute and save axisymmetric statistics
compute_mean = false; % Computes time-averaged values (use if fnum_st and fnum_end are different)
convert_mat = false; % Saves a mat file of the output
save_figs = false;

% File selection
fnum_st = 1; % Choose starting file number
fnum_end = 1; % Choose ending file number

% Plotting
% size(ustore) = [176 176 80 10]
% 176x176 horizontal grid, 80 grid levels, 10 time indices
lev = 10; % Choose height index for plotting (e.g., 5 is the 5th grid level) 
tme = 1; % time index to plot (can't exceed num_times)

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the LES grid %
%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the script reads in the binary file containing the LES model
% grid information
cd(data_dir)
fid = fopen('fort.10_2', 'r', 'ieee-le');
fread(fid, 1, 'int');
Nx = fread(fid, 1, 'int'); % Number of points in the x-dim
Ny = fread(fid, 1, 'int'); % Number of points in the y-dim
Nz = fread(fid, 1, 'int'); % Number of points in the z-dim
fread(fid, 2, 'real*4');
Xc = fread(fid, Nx*Ny*Nz, 'real*4'); % x coordinates
fread(fid, 2, 'real*4');
Yc = fread(fid, Nx*Ny*Nz, 'real*4'); % y coordinates
fread(fid, 2, 'real*4');
Zc = fread(fid, Nx*Ny*Nz, 'real*4'); % z coordinates

% Dimensionalize variables
Xc = squeeze(Xc(1:Nx));
Yc = Xc;
Zc = squeeze(Zc(1+(Nx-1)*(Ny-1) : (Nx-1)*(Ny-1) : (Nz)*(Nx-1)*(Ny-1)+1));
Xc = Xc * v0^2 / 9.80;
Yc = Yc * v0^2 / 9.80;
Zc = Zc * v0^2 / 9.80;

xc(1:max(size(Xc)), 1, 1) = Xc;
yc(1, 1:max(size(Yc)), 1) = Yc;
zc(1, 1, 1:max(size(Zc))) = Zc;

% Form a three-dimensional grid
Xm = repmat(xc, [1 Ny Nz]);
Ym = repmat(yc, [Nx 1 Nz]);
Zm = repmat(zc, [Nx Ny 1]);

fclose(fid);

% These variables define the minimum and maximum values of the indices
% saved. For 1,Nx, all of the data are saved but some times only we save a
% subdomain to save space.
ix1 = 1;
ix2 = Nx;
iy1 = 1;
iy2 = Nx;
iz1 = 1;
iz2 = Nz;

% Reconstruct for the subdomain
Xmf = Xm(ix1:ix2, iy1:iy2, iz1:iz2); 
Ymf = Ym(ix1:ix2, iy1:iy2, iz1:iz2);
Zmf = Zm(ix1:ix2, iy1:iy2, iz1:iz2);

% Loops through the LES data

u_LES = zeros(size(Xmf,1), size(Xmf,2), size(Xmf,3), num_times*fnum_end);
v_LES = u_LES;
w_LES = u_LES;
tke_LES = u_LES;
p_LES = u_LES;
t_LES = zeros(1,num_times*fnum_end);

x_LES = Xmf;
y_LES = Ymf;
z_LES = Zmf;

r_LES = sqrt(Xmf.^2 + Ymf.^2);
phi_LES = atan(Xmf ./ Ymf); % Azimuth angles w.r.t TOR center (domain center)
% adjust angle values for each quadrant
phi_LES(Ymf < 0) = phi_LES(Ymf < 0) + pi; % Lower 2 quadrants
phi_LES(Xmf < 0 & Ymf >= 0) = phi_LES(Xmf < 0 & Ymf >= 0) + 2*pi; % Upper left quadrant


for fnum = fnum_st:fnum_end
    %LES file name
    if name_type == 1
        LES_name = ['LES_mean_1_6_fnum' num2str(fnum) '.dat'];
    elseif name_type == 2
        LES_name = ['LES_mean_fnum' num2str(fnum) '.dat'];
    elseif name_type == 3
        LES_name = ['LES_mean_1 _6 _fnum' num2str(fnum) '.dat'];
    end

    % Spatial dimensions are scaled by a factor of v0^2/9.8, so increasing the
    % non-dimensional velocity scale increases the domain size.  Time also
    % increases with increasing velocity scale by a factor of v0/g.

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Read in the LES file %
    %%%%%%%%%%%%%%%%%%%%%%%%

    % Opens the file
    fid = fopen(LES_name, 'r', 'ieee-le');
    tmp = fread(fid, 1, 'int');
    
    time = zeros(1, num_times);
    ustore = zeros(size(Xmf,1), size(Xmf,2), size(Xmf,3), num_times);
    vstore = ustore;
    wstore = ustore;
    pstore = ustore;
    tkestore = ustore;
    urstore = ustore;
    vtstore = ustore;
    Lstore = ustore;

    % Loops through the LES data
    for kdx = 1:num_times
        % LES model time
        time(kdx) = fread(fid, 1, 'real*4') * v0 / g;
        tmp = fread(fid, 2, 'int');

        % U component of the wind (cartesian)
        u = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0 + vtrans;
        tmp = fread(fid, 2, 'int');

        % V component of the wind (cartesian)
        v = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0;
        tmp = fread(fid, 2 ,'int');

        % W component of the wind (cartesian)
        w = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0;
        tmp = fread(fid, 2, 'int');
        
        % Perturbation pressure 
        p = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0^2;
        tmp = fread(fid, 2, 'int');

        % Turbulent kinetic energy (resolved component)
        tke = single(fread(fid, max(size(Xmf(:))), 'real*4')) * v0^2;
        tmp = fread(fid, 2, 'int');
        

        % Reshape the matrices for a three-dimensional grid
        u = reshape(u, [size(Xmf,1), size(Xmf,2), size(Xmf,3)]);
        v = reshape(v, [size(Xmf,1), size(Xmf,2), size(Xmf,3)]);
        w = reshape(w, [size(Xmf,1), size(Xmf,2), size(Xmf,3)]);
        p = reshape(p, [size(Xmf,1), size(Xmf,2), size(Xmf,3)]);
        tke = reshape(tke, [size(Xmf,1), size(Xmf,2), size(Xmf,3)]);
        
        % Ur component of the wind (polar)
        ur = u.*sin(phi_LES) + v.*cos(phi_LES);
        
        % Vt component of the wind (polar)
        vt = v.*sin(phi_LES) - u.*cos(phi_LES);
        
        % Angular momentum
        L = vt .* r_LES;

        % Create a 4-D matrix (x,y,z,t) to store the data 
        ustore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = u; 
        vstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = v;
        wstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = w;
        pstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = p;
        tkestore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = tke;
        urstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = ur;
        vtstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = vt;
        Lstore(1:size(Xmf,1), 1:size(Xmf,2), 1:size(Xmf,3), kdx) = L;
        clear u v w p tke ur vt L;
    end
    fclose(fid);
    cd(base_dir)

    % Plot u,v,w,tke at the selected height level and time on lines 32 ? 33
    figure(1)
    feval('boonlib', 'bsizewin', 1, [700 600])
    
    subplot(2,2,1)
    tmp = squeeze(double(ustore(:,:,lev,tme)));
    hs = pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp);
    colormap(blib('rbmap'))
    shading flat
    axis square
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    caxis([-nanmax(abs(tmp(:))) nanmax(abs(tmp(:)))])
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    %title(['u (m s^{-1}) at z = ' num2str(roundn(Zm(1,1,lev), 1)) ' m'], 'FontSize', 14)
    title('(a) U wind', 'FontSize', 14)
    
    subplot(2,2,2)
    tmp = squeeze(double(vstore(:,:,lev,tme)));
    hs(2) = pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp);
    colormap(blib('rbmap'))
    shading flat
    axis square
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([-nanmax(abs(tmp(:))) nanmax(abs(tmp(:)))])
    %title(['v (m s^{-1}) at z = ' num2str(roundn(Zm(1,1,lev), 1)) ' m'], 'FontSize', 14)
    title('(b) V wind', 'FontSize', 14)
    
    subplot(2,2,[3,4])
    tmp = squeeze(double(wstore(:,:,lev,tme)));
    hs(3) = pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp);
    colormap(blib('rbmap'))
    shading flat
    axis square
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([-nanmax(abs(tmp(:))) nanmax(abs(tmp(:)))])
    %title(['w (m s^{-1}) at z = ' num2str(roundn(Zm(1,1,lev), 1)) ' m'], 'FontSize', 14)
    title('(c) W wind', 'FontSize', 14)
    
%     ha = subplot(2,2,4);
%     tmp=squeeze(double(tkestore(:,:,lev,tme)));
%     hs(4) = pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp);
%     colormap(ha, blib('rbmap'))
%     shading flat
%     colorbar
%     xlabel('x (m)')
%     ylabel('y (m)')
%     axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
%     %title(['TKE (m^2 s^{-2}) at z = ' num2str(roundn(Zm(1,1,lev), 1)) ' m'], 'FontSize', 14)
%     title('(d) TKE', 'FontSize', 14)
    
    %set(ha, 'DataAspect', [1 1 1])
    set(hs, 'EdgeColor', 'none')
    
    set(gcf, 'PaperPositionMode', 'auto', 'renderer', 'zbuffer'); % This makes the plot look nicer
    if save_figs
        cd(save_dir)
            print(['LES_z' num2str(lev) '_f' num2str(fnum) '_t' num2str(tme)], '-dpng');
        cd(base_dir)
    end
    
    
    figure(2)
    
    subplot(1,3,1)
    tmp = squeeze(mean(double(urstore(:,:,:,tme)),3));
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    colormap(blib('rbmap'))
    shading flat
    colorbar
    xlabel('x (m)')
    ylabel('y (m)')
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([-nanmax(abs(tmp(:))) nanmax(abs(tmp(:)))])
    title('(a) Height-averaged u_r', 'FontSize', 14)
    
    subplot(1,3,2)
    tmp = squeeze(mean(double(vtstore(:,:,:,tme)),3));
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    colormap(blib('rbmap'))
    shading flat
    colorbar
    xlabel('x (m)')
    ylabel('y (m)')
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([-nanmax(abs(tmp(:))) nanmax(abs(tmp(:)))])
    title('(b) Height-averaged v_t', 'FontSize', 14)
    
    subplot(1,3,3)
    tmp = squeeze(mean(double(Lstore(:,:,:,tme)),3));
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    xlabel('x (m)')
    ylabel('y (m)')
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([0 nanmax(abs(tmp(:)))])
    title('(c) Height-averaged L', 'FontSize', 14)

    set(gcf, 'PaperPositionMode', 'auto'); % This makes the plot look nicer
    if save_figs
        cd(save_dir)
            print(['LESpolar_z' num2str(lev) '_f' num2str(fnum) '_t' num2str(tme)], '-dpng');
        cd(base_dir)
    end
    
    figure(3)
    tmp = squeeze(mean(double(Lstore(:,:,:,tme)),3));
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    axis square
    c = colorbar;
    c.Label.String = 'm^2 s^{-1}';
    c.Label.FontSize = 13;
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    axis(view_fact * [-nanmax(Xmf(:)) nanmax(Xmf(:)) -nanmax(Xmf(:)) nanmax(Xmf(:))])
    caxis([0 nanmax(abs(tmp(:)))])
    title('Vertically averaged {\it L}', 'FontSize', 14)
    
    set(gcf, 'PaperPositionMode', 'auto'); % This makes the plot look nicer
    if save_figs
        cd(save_dir)
            print('LES_ang_mom_mean', '-dpng');
        cd(base_dir)
    end
    
    xmin = 0;
    ymin = 0;
    xmin = ones(size(Zm,3), 1) * xmin;
    ymin = ones(size(Zm,3), 1) * ymin;
    xtmp2 = 1:max(size(Xm));
    ytmp2 = 1:max(size(Ym));
    umn = nanmean(ustore, 4);
    vmn = nanmean(vstore, 4);
    wmn = nanmean(wstore, 4);
    tkemn = nanmean(tkestore, 4);
    pmn = nanmean(pstore, 4);
    urmn = nanmean(urstore, 4);
    vtmn = nanmean(vtstore, 4);
    Lmn = nanmean(Lstore, 4);
    avg_file = 1;

    cd(base_dir)
    if axi_stats_flag
        [ur_mean, vr_mean, wr_mean, pr_mean, ang_mom_mean, tke_mean, stats, rtmp2, ztmp2] = ...
            axy_stats_v3(Xm, Ym, Zm, xtmp2, ytmp2, xmin, ymin, umn, vmn, wmn, pmn, tkemn, avg_file, v0);
        ur_save(:,:,fnum) = ur_mean;
        vr_save(:,:,fnum) = vr_mean;
        wr_save(:,:,fnum) = wr_mean;
        pr_save(:,:,fnum) = pr_mean;
        tke_save(:,:,fnum) = tke_mean;
        L_save(:,:,fnum) = permute(ang_mom_mean, [2 1]);
    end
    
    figure()
    
    subplot(1,3,1)
    pcolor(rtmp2, ztmp2, permute(ur_mean, [2 1]))
    shading flat
    colorbar
    title('Mean u_r')
    xlabel('x (m)')
    ylabel('z (m)')
    
    subplot(1,3,2)
    pcolor(rtmp2, ztmp2, permute(vr_mean, [2 1]))
    shading flat
    colorbar
    title('Mean v_t')
    xlabel('x (m)')
    ylabel('z (m)')
    
    subplot(1,3,3)
    pcolor(rtmp2, ztmp2, permute(wr_mean, [2 1]))
    shading flat
    colorbar
    title('Mean w')
    xlabel('x (m)')
    ylabel('z (m)')
    
    
    
    if compute_mean
        if ~exist('umn_save', 'var')
            umn_save = zeros(size(umn));
            vmn_save = umn_save;
            wmn_save = umn_save;
            pmn_save = umn_save;
            tkemn_save = umn_save;
            Lmn_save = umn_save;
        end
       umn_save = umn_save + umn / (fnum_end - fnum_st + 1);
       vmn_save = vmn_save + vmn / (fnum_end - fnum_st + 1);
       wmn_save = wmn_save + wmn / (fnum_end - fnum_st + 1);
       pmn_save = pmn_save + pmn / (fnum_end - fnum_st + 1);
       tkemn_save = tkemn_save + tkemn / (fnum_end - fnum_st + 1);
       Lmn_save = Lmn_save + Lmn / (fnum_end - fnum_st + 1);
    end
    
    u_LES(:, :, :, fnum*10-9: fnum*10) = ustore;
    v_LES(:, :, :, fnum*10-9: fnum*10) = vstore;
    w_LES(:, :, :, fnum*10-9: fnum*10) = wstore;
    tke_LES(:, :, :, fnum*10-9: fnum*10) = tkestore;
    p_LES(:, :, :, fnum*10-9: fnum*10) = pstore;
    t_LES(fnum*10-9: fnum*10) = time;
    
    if convert_mat
        cd(save_dir)
        if fnum == fnum_st
            save('grid.mat', 'Xmf', 'Ymf', 'Zmf');
        end
        save(['LES_' num2str(fnum) '.mat'], 'ustore', 'vstore', 'wstore', 'tkestore', 'pstore', 'urstore', 'vtstore', 'Lstore', 'time');
        cd(base_dir)
    end
    %clear ustore vstore pstore wstore tkestore urstore vtstore Lstore time;
end

%save([save_dir '/LES_all.mat'], 'x_LES', 'y_LES', 'z_LES',...
%    't_LES', 'u_LES', 'v_LES', 'w_LES', 'tke_LES', 'p_LES', '-v7.3', '-nocompression')

if compute_mean
    cd(save_dir)
    save('mean.mat', 'umn_save', 'vmn_save', 'wmn_save', 'pmn_save', 'tkemn_save', 'Lmn_save');
    
    lev = 3;
    
    figure(151)
    cmax = nanmax(nanmax(abs(double(umn_save(:,:,lev)))));
    tmp = squeeze(double(umn_save(:,:,lev)));
    tmp(abs(tmp) < 0.001) = NaN;
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    caxis([-cmax cmax])
    title(['Mean U (m/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)
    if save_figs
        print([save_dir '/' 'umean_t0_z' num2str(roundn(Zm(1,1,lev), -1))], '-dpng')
    end
    
    figure(152)
    cmax = nanmax(nanmax(abs(double(vmn_save(:,:,lev)))));
    tmp = squeeze(double(vmn_save(:,:,lev)));
    tmp(abs(tmp) < 0.001) = NaN;
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    caxis([-cmax cmax])
    title(['Mean V (m/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)
    if save_figs
        print([save_dir '/' 'vmean_t0_z' num2str(roundn(Zm(1,1,lev), -1))], '-dpng')
    end
    
    figure(153)
    cmax = nanmax(nanmax(abs(double(wmn_save(:,:,lev)))));
    tmp = squeeze(double(wmn_save(:,:,lev)));
    tmp(abs(tmp) < 0.001) = NaN;
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    caxis([-cmax cmax])
    title(['Mean W (m/s) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)
    if save_figs
        print([save_dir '/' 'wmean_t0_z' num2str(roundn(Zm(1,1,lev), -1))], '-dpng')
    end
    
    figure(154)
    cmax = nanmax(nanmax(abs(double(tkemn_save(:,:,lev)))));
    tmp = squeeze(double(tkemn_save(:,:,lev)));
    tmp(abs(tmp) < 0.001) = NaN;
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    caxis([0 cmax])
    title(['Mean TKE (m^2/s^2) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)
    if save_figs
        print([save_dir '/' 'tkemean_t0_z' num2str(roundn(Zm(1,1,lev), -1))], '-dpng')
    end
    
    figure(155)
    cmax = nanmax(nanmax(abs(double(pmn_save(:,:,lev)))));
    tmp = squeeze(double(pmn_save(:,:,lev)));
    tmp(abs(tmp) < 0.001) = NaN;
    pcolor(squeeze(double(Xmf(:,:,lev))), squeeze(double(Ymf(:,:,lev))), tmp)
    shading flat
    colorbar
    caxis([-cmax cmax])
    title(['Mean P'' (Pa) at z = ' num2str(roundn(Zm(1,1,lev), -1)) ' m'], 'FontSize', 14)
    if save_figs
        print([save_dir '/' 'pmean_t0_z' num2str(roundn(Zm(1,1,lev), -1))], '-dpng')
    end
    
end

if axi_stats_flag
    cd(save_dir)
        save('axy.mat', 'ur_save', 'vr_save', 'wr_save', 'pr_save', 'L_save', 'tke_save', 'rtmp2', 'ztmp2');
    cd(base_dir)
end


