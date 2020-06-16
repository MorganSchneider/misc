%This code plots *.mat files that have been converted from WSR-88D data,
%and *.nc files from PX-1000 or PX-10K data.

clear;
close all; 

case_dir = uigetdir('~/Documents/nexrad'); %choose the folder containing case study netCDF files
dir_name = [case_dir '/mat/']; %gets the name of the directory where *.mat files are located

base_dir = '~/Documents/'; %This is your main directory where this script is saved

data_type = input('Enter data type: 1) .mat, 2) PX-1K/10K .nc ');
ouprime_flag = false;

rhohv_map; %sets up the color bar for rhohv

big_monitor = true; %big monitor is a flag I use if I have a separate large monitor, use false for a laptop
ext_monitor = true;

%plot_choice determines the variables that will be plotted on the screen
if big_monitor
    plot_choice = [1 2 3 4 6 5]; 
    winsize = [700 600]; %Winsize is the [width,height] size of all the figures in pixels;
else
    plot_choice = [1 2 3 4]; 
    winsize = [500 400]; %Winsize is the [width,height] size of all the figures in pixels;
end
%1: Reflectivity
%2: Radial Velocity
%3: ZDR
%4: rhohv
%5: TDS thresholded image
%6: Phidp
%7: KDP (may not be available)
%8: Specific attenuation (may not be available)
%9: Spectrum width

data_cursor_mode = true; 
att_corr_flag = false; %turns on attenuation correction and phidp filtering
KOUN_adjust_flag = false;

load carbone42.mat;

Z_thres = 43;
rhohv_thres = 0.9; 

%These are limits for the color bars; vr_lims sets the limits of the color bars for velocity
vr_lims = [-25 25];
zdr_lims = [-3 5];
vr_lims2 = [-60 60];

%Flags for plotting damage and mapping over the figure
damage_plot = false;
add_map = false; 

%These are flags used in bmapover
bmap_flags = [1 1 1 1 1 1 1 1 0 0 1 0];

%Set figure positions and sizes
if big_monitor
    adjust_pos(1) = 1500; %Shifts the figure 1500 pixels to the right
    adjust_pos(2) = -500; %Shifts the figure winsize(1) pixels to the right
else
    adjust_pos = [0 0];
end

%These numbers control the size and position of the figure windows
pos1 = [1 + adjust_pos(1), adjust_pos(2), winsize];
pos2 = [winsize(1) + 1 + adjust_pos(1), adjust_pos(2), winsize];
pos3 = [1 + adjust_pos(1), 1 + winsize(1), winsize];
pos4 = [1 + adjust_pos(1) + winsize(1), 1 + winsize(1), winsize];
pos5 = [1 + adjust_pos(1) + 2 * winsize(1), adjust_pos(2), winsize];
pos6 = [1 + adjust_pos(1), 1 + winsize(1), 700, 90];
pos7 = [winsize(1) + adjust_pos(1), adjust_pos(2), 260, 100];
pos8 = pos7;
pos9 = [1 + adjust_pos(1) + 2 * winsize(1), 1 + winsize(1), winsize];
pos10 = pos6;
pos10(3) = pos10(3) - 300;

%Sets different parameters for plotting
axis_label_size = 16;
options = 8;
is_running = true;
title_size = 16;
zoom_factor_x = 3;
zoom_factor_y = 3;
p1_check = false;
p2_check = false;
p3_check = false;
p4_check = false;
p5_check = false;
p6_check = false;
new_color_flag = [1 1 1 1 1 1];
ae = 6378.1370;

%ZDR calibration
zdr_calib = 0; %This may need to be changed if ZDR is poorly calibrated

%Set up different colormaps
cmap3 = carbone42;
cmap4 = flipud([0.4 0.4 0.4; 0.5 0.5 0.5; 0.6 0.6 0.6; 0.7 0.7 0.7; 0.8 0.8 0.8]);
cmap3 = [cmap3(1, 1:size(cmap3,2)); cmap4; cmap3(2:size(cmap3,1), 1:size(cmap3,2))]; 
cmap5 = feval('boonlib', 'rgmap', 32);

%This loops through the current directory and lists all of the files.  This
%loop obtains all of the file times and elevation scan numbers
cd(dir_name)
    if data_type == 1
        %List the mat files in the directory
        files = dir('*.mat'); 
        %This is the starting window coordinates
        lims = [-125 125 -125 125]; 
        
        nFiles = max(size(files));
        scan_number = zeros(1, nFiles);
        file_time_save = cell(1, nFiles);
        ele_save = zeros(1, nFiles);
        for idx = 1:nFiles
            file_tmp = files(idx).name;
            file_name_length = max(size(file_tmp));
            file_time = char(file_tmp(file_name_length-12 : file_name_length-7));  
            scan_number(idx) = idx;
            file_time_save{idx} = file_time;
            ele_save(idx) = str2double(file_tmp(file_name_length-5 : file_name_length-4)); 
        end
    else
        %List the file names in each directory by file type
        files = dir('*Z.nc');
        files_D = dir('*D.nc');
        files_V = dir('*V.nc');
        files_P = dir('*P.nc');
        files_W = dir('*W.nc');
        files_R = dir('*R.nc');
        if max(size(files)) == max(size(files_P)) && max(size(files)) == max(size(files_W)) && max(size(files)) == max(size(files_D)) && max(size(files)) == max(size(files_V))
            dir_ok = true;
        else
            dir_ok = false;
            keyboard;
        end
        %This is the starting window coordinates
        lims = [-50 50 -50 50]; 
        
        nFiles = max(size(files));
        scan_number = zeros(1, nFiles);
        file_time_save = cell(1, nFiles);
        ele_save = zeros(1, nFiles);
        for idx = 1:nFiles
            file_tmp = files(idx).name;
            time_tmp = ncreadatt(file_tmp, '/', 'Time');
            time_tmp = datenum([1970 1 1 0 0 double(time_tmp)]);
            scan_number(idx) = idx;
            file_time_save{idx} = datestr(time_tmp, 'HHMMSS');
            ele_save(idx) = nanmean(ncread(file_tmp, 'Elevation'));
        end
    end
cd(base_dir)

num_loops = max(size(plot_choice));
p1 = gobjects(1, num_loops);
t1 = gobjects(1, num_loops);

while (is_running)
    %To load a certain tilt, just find(scan_number==desired scan_number)
    if options == 8
        fprintf('Choose a time: \n');
        for idx = 1:max(size(file_time_save))
            fprintf([num2str(idx) ': ' num2str(round(ele_save(idx) * 10) / 10) ' deg ' file_time_save{idx}(1:2) ':' file_time_save{idx}(3:4) ':' file_time_save{idx}(5:6) '\n']);
        end
        fprintf([num2str(idx + 1) ': Previous Settings' '\n']);
        scan_num = input('Choose scan number: ');

        if scan_num == idx + 1
            load([dir_name 'last_settings_edge.mat']);
            [row] = find(scan_number == scan_num);
            time_choice = file_time_save{scan_num}(1:6);
        else
            [row] = find(scan_number == scan_num);
            time_choice = file_time_save{scan_num}(1:6);
        end
    end
    
    %Read in the data
    if data_type == 1
        cd(dir_name)
        load(files(row).name);
        cd(base_dir)
    else
        read_PX;
    end

    %Set colorbar limits
    if data_type == 1
        if sum(ismember(plot_choice, 9)) > 0
            data.SW = data.SW_h;
            sw_lims = [0 15];
        end
        if KOUN_adjust_flag
           zdr_lims = [-5 8]; 
        end
        clims3 = [-15 75];
    end
    if ~isfield(data, 'detr') || ~isfield(data, 'gw')
        data.detr = 250;
        clims3 = [-15 65];
    end

%     if radar_choice < 3
        %Apply any data corrections needed, replace missing data as NaN
        if data_type == 1
            [row2] = find(abs(imag(data.rho_hv)) > 0);
            data.rho_hv(row2) = NaN;
            if att_corr_flag
                [row2] = find(abs(imag(data.ZDR_corr)) > 0);
                data.ZDR_corr(row2) = NaN;
                clear row2;
                data.ZDR_corr = data.ZDR_corr + zdr_calib;
                data.ZDR_corr(data.ZDR_corr < -900) = NaN;
                data.Z_corr(data.Z_corr < -900) = NaN;
            end
            data.Z_H(data.Z_H < -900) = NaN;
            data.Z_DR(data.Z_DR < -900) = NaN;
            data.phi_dp(data.phi_dp < -900) = NaN;

            if isfield(data, 'phi_dp_filtered')
                data.phi_dp_filtered(data.phi_dp_filtered < -900) = NaN;
            end
            data.rho_hv(data.rho_hv < -900) = NaN;
    %     end

            if isfield(data, 'va')
                vr_lims = [-data.va data.va];
            end
        end
    
    %If a radar name hasn't been provided in the file, it will ask for it
    %here
    if ~isfield(data, 'Station')
        if radar_choice == 1
            data.Station = 'OU-PRIME';
        else
            radar_input = input('Choose radar: (1) OU-PRIME, (2) KOUN, (3) other: ');
            if radar_input == 1
                data.Station = 'OU-PRIME';
            elseif radar_input == 2
                data.Station = 'KOUN';
            else
                data.Station = input('Input radar name: ');
            end
        end
    end
    if strcmp(data.Station, 'NOP4')
        radar = 'KCRI';
    else
        radar = data.Station;
    end
%     keyboard;
    
    %Rename variables
    az = data.az; 
    if isfield(data, 'el')
        ele = median(data.el);
    else
        ele = median(data.ele);
    end

    %Get the number of range gates
    if size(data.vr_h, 2) == max(size(az))
        num_gates = size(data.vr_h, 1);
    else
        num_gates = size(data.vr_h, 2);
    end
    
    %Convert range and azimuth to x and y for plotting
    if data_type == 1
        if ouprime_flag
            r = (data.detr + data.detr * ((1:num_gates) - 1)) / 1000;
        else
            r = (data.detr * 8 + data.detr * ((1:num_gates) - 1)) / 1000;
        end
    else
        if size(data.Z_H, 1) == max(size(data.az))
            az_dims = size(data.Z_H, 1);
        else
            gat_dims = size(data.Z_H, 1);
        end
        if size(data.Z_H, 2) == max(size(data.az))
            az_dims = size(data.Z_H, 2);
        else
            gat_dims = size(data.Z_H, 2);
        end
        r = (data.gw(1) : data.gw(1) : gat_dims * data.gw(1)) / 1000;
    end
    az_rad = az * pi/180;
    az_rad = repmat(az_rad, 1, max(size(r)));
    r_km = repmat(r, size(data.az, 1), 1);
    if data_cursor_mode
        xx = r_km .* sin(az_rad) * cos(ele*pi/180);
        yy = r_km .* cos(az_rad) * cos(ele*pi/180);
        if isfield(data, 'az2')
            az_rad = data.az2 * pi / 180;
            az_rad = repmat(az_rad, 1, max(size(r)));
            r_km = repmat(r, size(data.az2, 1), 1);
            xx2 = r_km .* sin(az_rad);
            yy2 = r_km .* cos(az_rad);
        elseif isfield(data, 'az_vr')
            az_rad = data.az_vr * pi / 180;
            az_rad = repmat(az_rad, 1, max(size(r)));
            r_km = repmat(r, size(data.az_vr, 1), 1);
            xx2 = r_km .* sin(az_rad) * cos(ele*pi/180);
            yy2 = r_km .* cos(az_rad) * cos(ele*pi/180);
            xx2 = xx2';
            yy2 = yy2';
        end
    else
        xx = r_km .* sin(az_rad) * cos(ele*pi/180);
        yy = r_km .* cos(az_rad) * cos(ele*pi/180);
%         if ~ouprime_flag
%              xx = xx';
%              yy = yy';
%         end
        if isfield(data, 'az2')
            az_rad = data.az2 * pi / 180;
            az_rad = repmat(az_rad, 1, max(size(r)));
            r_km = repmat(r, size(data.az2, 1), 1);
            xx2 = r_km .* sin(az_rad) * cos(ele*pi/180);
            yy2 = r_km .* cos(az_rad) * cos(ele*pi/180);
        elseif isfield(data, 'az_vr')
            az_rad = data.az_vr * pi / 180;
            az_rad = repmat(az_rad, 1, max(size(r)));
            r_km = repmat(r, size(data.az2, 1), 1);
            xx2 = r_km .* sin(az_rad) * cos(ele*pi/180);
            yy2 = r_km .* cos(az_rad) * cos(ele*pi/180);
        end
    end
    %Transpose data if needed
    if size(xx, 1) ~= size(data.vr_h, 1)
        xx = xx';
        yy = yy';
    end
    

    %Loops through for each plot
    
    for pdx = 1:num_loops
        %Check to see if the figure exists already
        if pdx == 1
            if p1_check
                fig_exists = true;
            else
                fig_exists = false;
            end
        elseif pdx == 2
            if p2_check
                fig_exists = true;
            else
                fig_exists = false;
            end
        elseif pdx == 3
            if p3_check
                fig_exists = true;
            else
                fig_exists = false;
            end
        elseif pdx == 4
            if p4_check
                fig_exists = true;
            else
                fig_exists = false;
            end
        elseif pdx == 5
            if p5_check
                fig_exists = true;
            else
                fig_exists = false;
            end
        end
        %Create a new figure if it doesn't exist
        if ~fig_exists
            figure(pdx);

            %Set the figure size
            feval('boonlib', 'bsizewin', pdx, winsize)
            
            %Radar reflectivity factor
            if plot_choice(pdx) == 1
                if att_corr_flag
                    data_tmp = double(data.Z_corr);
                else
                    data_tmp = double(data.Z_H(1:size(data.vr_h, 1), 1:size(data.vr_h, 2)));
                end
                p1(pdx) = pcolor(xx, yy, data_tmp);
                colormap(cmap3)
                caxis(clims3)
                title([num2str(round(ele*10)/10) '^{\circ} Reflectivity (dBZ) at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
            
            %Doppler velocity
            elseif plot_choice(pdx) == 2
                if isfield(data, 'vr_dealiased')
                    data_tmp = double(data.vr_dealiased);
                elseif isfield(data, 'fold_int_map')
                    data_tmp = double(data.vr_h) + double(data.fold_int_map) * double(data.va);
                else
                    data_tmp = double(data.vr_h);
                end
                if exist('xx2', 'var')
                    p1(pdx) = pcolor(xx2, yy2, data_tmp);
                else
                    p1(pdx) = pcolor(xx, yy, data_tmp);
                end
                colormap(cmap5)
                if nanmax(nanmax(abs(data_tmp))) > 25
                    caxis(vr_lims2)
                else
                    caxis(vr_lims)
                end
                title([num2str(round(ele*10)/10) '^{\circ} Radial Velocity at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
           
            %Correlation coefficient (rhohv)
            elseif plot_choice(pdx) == 3
                data_tmp = double(data.rho_hv);
                p1(pdx) = pcolor(xx, yy, data_tmp);
                colormap(cmap2)
                caxis([0 1])
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ $\rho_{  HV}$ at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
            
            %Differential reflectivity
            elseif plot_choice(pdx) == 4
                if att_corr_flag
                    data_tmp = double(data.ZDR_corr);
                else
                    data_tmp = double(data.Z_DR);
                end
                p1(pdx) = pcolor(xx, yy, data_tmp);
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ $Z_{DR}$ at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')            
                colormap(cmap3)
                caxis(zdr_lims)
            
            %Thresholded TDS plot
            elseif plot_choice(pdx) == 5
                if att_corr_flag
                    yn1 = data.Z_corr > Z_thres & data.rho_hv < rhohv_thres;
                    yn2 = data.Z_corr > Z_thres & data.rho_hv < 0.5;
                    yn3 = data.Z_corr > Z_thres & data.rho_hv < 0.3;
                    yn4 = data.Z_corr > Z_thres & data.rho_hv < 0.1;
                else
                    yn1 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < rhohv_thres;
                    yn2 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.5;
                    yn3 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.3;
                    yn4 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.1;
                end
                data_tmp = double(yn1) + double(yn2) + double(yn3) + double(yn4);
                p1(pdx) = pcolor(xx, yy, data_tmp);
                caxis([0 4])
                colormap(jet(5))
                shading flat
                if add_map
                    boonlib('spidergrid')
                end
                hcb = colorbar('YTickLabel',...
                    {'>0.7'; '0.5-0.7'; '0.3-0.5'; '0.1-0.3'; '<0.1'});
                set(hcb, 'YTick', 0.4:0.8:3.6)
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ Thresholded image at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')            
            
            %Propagation Differential Phase (phidp)
            elseif plot_choice(pdx) == 6
                if att_corr_flag
                    data_tmp = double(data.phi_dp_filtered);
                else
                    if isfield(data, 'Phi_DP')
                        data_tmp = double(data.Phi_DP);
                    else
                        data_tmp = double(data.phi_dp);
                    end
                end
                p1(pdx) = pcolor(xx, yy, data_tmp);
                title([num2str(round(ele*10)/10) '^{\circ} \Phi_{DP} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                caxis([-180 180])
                
            %Specific differential phase (KDP)
            elseif plot_choice(pdx) == 7
                data_tmp = double(data.KDP);
                p1(pdx) = pcolor(xx2, yy2, data_tmp);
                title([num2str(round(ele*10)/10) '^{\circ} KDP at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                caxis([-2 10])
            
            %Attenuation    
            elseif plot_choice(pdx) == 8
                data_tmp = double(data.A_H);
                p1(pdx) = pcolor(xx, yy, data_tmp);
                title([num2str(round(ele*10)/10) '^{\circ} A_{H} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                caxis([-2 15])
               
            %Spectrum width
            elseif plot_choice(pdx) == 9
                data_tmp = double(data.SW);
                if exist('xx2', 'var')
                    p1(pdx) = pcolor(xx2, yy2, data_tmp);
                else
                    p1(pdx) = pcolor(xx, yy, data_tmp);
                end
                title([num2str(round(ele*10)/10) '^{\circ} SW at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
%                 caxis([-2 15])
                caxis(sw_lims)
                
            elseif plot_choice(pdx) == 10
                data_tmp = pcolor(xx, yy, double(sd_phidp(1:size(data.vr_h,1), 1:size(data.vr_h,2))) );
                title([num2str(round(ele*10)/10) '^{\circ} \sigma_{\Phi_{DP}} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                
            elseif plot_choice(pdx) == 11
                data_tmp = pcolor(xx, yy, double(sd_ZDR(1:size(data.vr_h,1), 1:size(data.vr_h,2))) );
                title([num2str(round(ele*10)/10) '^{\circ} \sigma_{Z_{DR}} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
            end
            
            %Add figure colorbars, axes, adjust figure limits and
            %properties
            shading flat
            if add_map
                boonlib('spidergrid')
            end
            if plot_choice(pdx) ~= 5
                hcb = colorbar;
            end
            xlabel('Zonal Distance (km)', 'FontSize', axis_label_size)
            ylabel('Meridonal Distance (km)', 'FontSize', axis_label_size)
            set(gca, 'DataAspectRatio', [1 1 1])
            if pdx == 1
                set(gcf,'Position',pos1)
                p1_check = true;
                h1 = gca;
                if ~exist('lims', 'var')
                    lims(1:2) = get(h1, 'xlim');
                    lims(3:4) = get(h1, 'ylim');
                end
                axis(lims)
            elseif pdx == 2
                set(gcf, 'Position', pos2)
                p2_check = true;
                h2 = gca;
                axis(lims)
            elseif pdx == 3
                set(gcf, 'Position', pos3)
                p3_check = true;
                h3 = gca;
                axis(lims)
            elseif pdx == 4
                set(gcf, 'Position', pos4)
                p4_check = true;
                h4 =gca;
                axis(lims)
            elseif pdx == 5
                set(gcf, 'Position', pos5)
                p5_check = true;
                h5 = gca;
                axis(lims)
            elseif pdx == 6
                set(gcf, 'Position', pos9)
                p6_check = true;
                h6 = gca;
                axis(lims)
            end

            %Add a background map
            if add_map
                if KOUN_adjust_flag
                    if strcmp(radar, 'OU-PRIME')
                        bmapover(gca, bmap_flags, 'OU-PRIME')
                    else
                        bmapover(gca, bmap_flags, radar)
                    end
                else
                    bmapover(gca, bmap_flags, {data.lon, data.lat, 'Radar', 'OK'})
                end
            end
        else
            figure(pdx)
            
            %Change x and y axes if needed
            if plot_choice(pdx) == 7 || (plot_choice(pdx) >= 9 && exist('xx2', 'var')) || (plot_choice(pdx) == 2 && exist('xx2', 'var'))
                set(p1(pdx), 'XData', cat(2, xx2, xx2(:,1)));
                set(p1(pdx), 'YData', cat(2, yy2, yy2(:,1)));
            else
                set(p1(pdx), 'XData', cat(2, xx, xx(:,1)));
                set(p1(pdx), 'YData', cat(2, yy, yy(:,1)));
            end
            xtmp = cat(2, xx, xx(:,1));
            ytmp = cat(2, yy, yy(:,1));
            if exist('xx2', 'var')
                xtmp2 = cat(2, xx2, xx2(:,1));
                ytmp2 = cat(2, yy2, yy2(:,1));
            end
            set(p1(pdx), 'ZData', double(zeros(size(xtmp,1), size(xtmp,2))) )
            %These change the background data but leave the figure axes and
            %other unchanged properties the same
            if plot_choice(pdx) == 1
                if att_corr_flag
                    set(p1(pdx), 'CData', double(data.Z_corr))
                else
                    set(p1(pdx), 'CData', double(data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2))) )
                end
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ Z (dBZ) at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')                
                if new_color_flag(pdx)
                    colormap(cmap3)
                    caxis(clims3)
                    colorbar
                    new_color_flag(pdx) = false;
                end
                
            elseif plot_choice(pdx) == 2
                if isfield(data, 'vr_dealiased')
                    set(p1(pdx), 'CData', double(data.vr_dealiased)); 
                elseif isfield(data, 'fold_int_map')
                    set(p1(pdx), 'CData', double(data.vr_h) + double(data.fold_int_map) * double(data.va)); 
                else
                    set(p1(pdx), 'CData', double(data.vr_h));
                end
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ $v_{r}$ at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')                
                if new_color_flag(pdx)
                    colormap(cmap5)
                    if isfield(data, 'vr_dealiased') || isfield(data, 'fold_int_map')
                        caxis(vr_lims2)
                    else
                        caxis(vr_lims)
                    end
                    colorbar
                    new_color_flag(pdx) = false;
                end
                
            elseif plot_choice(pdx) == 3
                set(p1(pdx), 'CData', double(data.rho_hv))
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ $\rho_{  HV}$ at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')                            
                if new_color_flag(pdx)
                    colormap(cmap2)
                    caxis([0 1])
                    colorbar
                    new_color_flag(pdx) = false;
                end
                
            elseif plot_choice(pdx) == 4
                if att_corr_flag
                    set(p1(pdx), 'CData', double(data.ZDR_corr))
                else
                    set(p1(pdx), 'CData', double(data.Z_DR))
                end
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '$^{\circ}$ $Z_{DR}$ at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')                
                if new_color_flag(pdx)
                    colormap(cmap3)
                    caxis(zdr_lims)
                    colorbar
                    new_color_flag(pdx) = false;
                end
                
            elseif plot_choice(pdx) == 5
                if att_corr_flag
                    yn1 = data.Z_corr > Z_thres & data.rho_hv < rhohv_thres;
                    yn2 = data.Z_corr > Z_thres & data.rho_hv < 0.5;
                    yn3 = data.Z_corr > Z_thres & data.rho_hv < 0.3;
                    yn4 = data.Z_corr > Z_thres & data.rho_hv < 0.1;
                else
                    yn1 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < rhohv_thres;
                    yn2 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.5;
                    yn3 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.3;
                    yn4 = data.Z_H(1:size(data.vr_h,1), 1:size(data.vr_h,2)) > Z_thres & data.rho_hv < 0.1;
                end
                data_tmp = double(yn1) + double(yn2) + double(yn3) + double(yn4); 
                set(p1(pdx), 'CData', data_tmp)
                title([num2str(round(ele*10)/10) '^{\circ} Thresholded image at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                if new_color_flag(pdx)
                    caxis([0 4])
                    colormap(jet(5))
                    hcb = colorbar('YTickLabel',...
                        {'>0.82'; '0.5-0.82'; '0.3-0.5'; '0.1-0.3'; '<0.1'});
                    set(hcb, 'YTick', 0.4:0.8:3.6)
                    new_color_flag(pdx) = false;
                end
                
            elseif plot_choice(pdx) == 6
                if att_corr_flag
                    set(p1(pdx), 'CData', double(data.phi_dp_filtered))
                else
                    set(p1(pdx), 'CData', double(data.phi_dp))
                end
                t1(pdx) = title('Hi');
                set(t1(pdx), 'String', [num2str(round(ele*10)/10) '^{\circ} \Phi_{DP} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')                
                colormap(cmap3)
                caxis([-180 180])
                
            elseif plot_choice(pdx) == 7
                set(p1(pdx), 'CData', double(data.KDP))
                title([num2str(round(ele*10)/10) '^{\circ} KDP at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                caxis([-2 10])
                
            elseif plot_choice(pdx) == 8
                set(p1(pdx), 'CData', double(data.A_H))
                title([num2str(round(ele*10)/10) '^{\circ} A_{H} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                caxis([-2 15])
                
            elseif plot_choice(pdx) == 9
                set(p1(pdx), 'CData', double(data.SW))
                title([num2str(round(ele*10)/10) '^{\circ} SW at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
%                 caxis([-2 15])
                caxis(sw_lims)
                
            elseif plot_choice(pdx) == 10
                set(p1(pdx), 'CData', double(sd_phidp(1:size(data.vr_h,1), 1:size(data.vr_h,2))) )
                title([num2str(round(ele*10)/10) '^{\circ} \sigma_{\Phi_{DP}} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
                
            elseif plot_choice(pdx) == 11
                set(p1(pdx), 'CData', double(sd_ZDR(1:size(data.vr_h,1), 1:size(data.vr_h,2))) )
                title([num2str(round(ele*10)/10) '^{\circ} \sigma_{\Z_{DR}} at ' time_choice ' UTC'],...
                    'FontSize', title_size, 'Interpreter', 'tex', 'HorizontalAlignment', 'Center')
                colormap(cmap3)
            end
        end
    end

    options = 8;
    %Save latest settings
    if exist('lims', 'var')
        save('last_settings_edge.mat', 'lims', 'time_choice', 'ele_save', 'scan_num');
    end  
    height_km = sqrt((sum(lims(1:2)) / 2) ^ 2 + (sum(lims(3:4)) / 2) ^ 2) * sin(ele * pi / 180) + 0.020;
    
    %This script provides several menu-based options to change times, elevations, zoom, save figures, etc. 
    getmenu_edge;
    %Save settings for easy reload
    if exist('lims', 'var')
        save('last_settings_edge.mat', 'lims', 'time_choice', 'ele_save', 'scan_num');
    end  
    
    if options < 6
        clear xx yy data yn1 yn2 yn3 yn4 r_km vel data_tmp;
    end
end

            
        


