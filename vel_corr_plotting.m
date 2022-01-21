%% this is temporary because I don't want to rerun run_velocity_correction manually again

% x = 1;
% debris_types = [1 3];
% debris_concentrations = [10000 100000];
% elevation_angles = [0.5 2.0 5.0];
% sim_date = '210614';
% LES = 'twocell';
% base_dir = ['~/Documents/sims/' LES];
% plot_flag = 0;
% var_save_flag = 1;
% plot_save_flag = 0;
% 
% for dx = debris_types
%     data_dir = [base_dir '/' sim_date '/debris' num2str(dx)];
%     for ex = elevation_angles
%         for cx = debris_concentrations
%             sim_name = ['sim-PPI' num2str(ex,'%.1f') '-DCU-d' num2str(dx) 'n' num2str(cx) '.mat'];
%             filename = [data_dir '/' sim_name];
%             
%             run_velocity_correction;
%         end
%     end
% end


%% Load data
clear

% SIMRADAR FIGURES:
% 1 - Sample PPIs
% 2 - Sample spectra
% 3 - Example DPSDs + DCA
% 4 - Example corrected spectra
% 5 - Bias PPIs + histograms
% 6 - Scatter plots with bias/correction & rhohv
% KOUN FIGURES:
% 1 - PPIs (mark gates for sample spectra)
% 2 - Sample spectra
% 3 - 
% 4 - 
% 5 - 
simr_plot = [1 1 1 1 0 1];
koun_plot = [0 0 0 0 0];
plot_save = 0;

LES = 'twocell';
sim_date = '210614';
dts = [1 3]; %debris types
cns = [10000 100000]; %debris concentration
els = [0.5 2.0 5.0]; %elevation angles

% Load all SimRadar DPSDs into a structure
% dd=debris type, conc=concentration, el=elevation
dd(length(dts)).conc(length(cns)).el(length(els)) = struct();
for dt = dts
    for cn = 1:length(cns)
        for el = 1:length(els)
            dpsd_dir = ['~/Documents/sims/' LES '/' sim_date '/debris' num2str(dt) '/'];
            dpsd_fname = ['dpsd-PPI' num2str(els(el),'%.1f') '-DCU-d' num2str(dt) 'n' num2str(cns(cn)) '.mat'];
            tmp = load([dpsd_dir dpsd_fname]);
            dd(dt).conc(cn).el(el) = tmp;
        end
    end
end

for dt = dts
    for cn = 1:length(cns)
        for el = 1:length(els)
            dd(dt).conc(cn).el(el).vcorr.dca = dd(dt).conc(cn).el(el).vr_new.dca - dd(dt).conc(cn).el(el).vr_old;
            dd(dt).conc(cn).el(el).vcorr.var = dd(dt).conc(cn).el(el).vr_new.var - dd(dt).conc(cn).el(el).vr_old;
            dd(dt).conc(cn).el(el).vcorr.agg = dd(dt).conc(cn).el(el).vr_new.agg - dd(dt).conc(cn).el(el).vr_old;
        end
    end
end
vvx = dd(1).conc(1).el(1).vvx; %velocity axis for plotting spectra
va = round(dd(1).conc(1).el(1).params.va);

% Load KOUN data into a structure
obs_el = 0.5;
obs_dir = '~/Documents/code/obsdata/';
obs_fname = ['dpsd-KOUN_data_' num2str(obs_el,'%.1f') '.mat'];
%obs_fname = 'dpsd-KOUN_data.mat';
koun = load([obs_dir obs_fname]);

img_dir = '~/Documents/articles/2021/';




%% SimRadar plots

%%% SECTION 2 %%%


% Plot SimRadar PPIs
if simr_plot(1)
    dt = 3;
    cn = 100000;
    el = 2.0;
    
    dd_dir = ['~/Documents/sims/' LES '/' sim_date '/debris' num2str(dt) '/'];
    dd_fname = ['sim-PPI' num2str(el,'%.1f') '-DCU-d' num2str(dt) 'n' num2str(cn) '.mat'];
    ddsim = load([dd_dir dd_fname]);
    nd_dir = ['~/Documents/sims/' LES '/' sim_date '/nodebris/'];
    nd_fname = ['sim-PPI' num2str(el,'%.1f') '-U-nodebris.mat'];
    ndsim = load([nd_dir nd_fname]);
    
    rho_clims = [0.4 1];
    zdr_clims = [-5 5];
    vr_clims = [-va-1 va+1];
    
    
    figure(1) %this is the less fancy version
    clf
    
    ha = subplot(2,4,1);
    pcolor(ndsim.xx, ndsim.yy, squeeze(ndsim.zh));
    caxis([0 80])
    colormap(ha, blib('zmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'dBZ';
    c.Title.FontSize = 12;
    %c.Label.VerticalAlignment = 'top';
    c.Ticks = 0:20:80;
    title('(a) Z_H', 'FontSize', 14)
    %xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 14)
    
    ha(2) = subplot(2,4,2);
    hs(2) = pcolor(ndsim.xx, ndsim.yy, squeeze(ndsim.vr));
    caxis(vr_clims);
    colormap(ha(2), blib('rgmap2'))
    shading flat
    c = colorbar;
%     c.Label.String = 'm s^{-1}';
%     c.Label.FontSize = 12;
%     c.Label.VerticalAlignment = 'middle';
    c.Ticks = -100:25:100;
    title('(b) v_r', 'FontSize', 14)
    %xlabel('x (m)', 'FontSize', 12)
    %ylabel('y (m)', 'FontSize', 12)
    
    ha(3) = subplot(2,4,3);
    hs(3) = pcolor(ndsim.xx, ndsim.yy, squeeze(ndsim.zdr));
    caxis(zdr_clims)
    colormap(ha(3), blib('nwsdmap'))
    c = colorbar;
%     c.Label.String = 'dB';
%     c.Label.FontSize = 12;
%     c.Label.VerticalAlignment = 'top';
    c.Ticks = -5:2:5;
    shading flat
    title('(c) Z_{DR}', 'FontSize', 14)
    %xlabel('x (m)', 'FontSize', 12)
    %ylabel('y (m)', 'FontSize', 12)
    
    ha(4) = subplot(2,4,4);
    hs(4) = pcolor(ndsim.xx, ndsim.yy, squeeze(real(ndsim.rhohv)));
    caxis(rho_clims)
    colormap(ha(4), blib('nwsrmap'))
    colorbar
    shading flat
    title('(d) \rho_{HV}', 'FontSize', 14)
    %xlabel('x (m)', 'FontSize', 12)
    %ylabel('y (m)', 'FontSize', 12)
    
    ha(5) = subplot(2,4,5);
    hs(5) = pcolor(ddsim.xx, ddsim.yy, squeeze(ddsim.zh));
    caxis([0 80])
    colormap(ha(5), blib('zmap'))
    c = colorbar;
%     c.Label.String = 'dBZ';
%     c.Label.FontSize = 12;
%     c.Label.VerticalAlignment = 'top';
    c.Ticks = 0:20:80;
    shading flat
    title('(e) Z_H', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    ha(6) = subplot(2,4,6);
    hs(6) = pcolor(ddsim.xx, ddsim.yy, squeeze(ddsim.vr));
    caxis(vr_clims)
    colormap(ha(6), blib('rgmap2'))
    c = colorbar;
%     c.Label.String = 'm s^{-1}';
%     c.Label.FontSize = 12;
%     c.Label.VerticalAlignment = 'middle';
    c.Ticks = -100:25:100;
    shading flat
    title('(f) v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    %ylabel('y (m)', 'FontSize', 12)
    
    ha(7) = subplot(2,4,7);
    hs(7) = pcolor(ddsim.xx, ddsim.yy, squeeze(ddsim.zdr));
    caxis(zdr_clims)
    colormap(ha(7), blib('nwsdmap'))
    c = colorbar;
%     c.Label.String = 'dB';
%     c.Label.FontSize = 12;
%     c.Label.VerticalAlignment = 'top';
    c.Ticks = -5:2:5;
    shading flat
    title('(g) Z_{DR}', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    %ylabel('y (m)', 'FontSize', 12)
    
    ha(8) = subplot(2,4,8);
    hs(8) = pcolor(ddsim.xx, ddsim.yy, squeeze(real(ddsim.rhohv)));
    caxis(rho_clims)
    colormap(ha(8), blib('nwsrmap'))
    colorbar
    shading flat
    title('(h) \rho_{HV}', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    %ylabel('y (m)', 'FontSize', 12)
    
    %set(ha, 'DataAspect', [1 1 1])
    %set(gcf, 'Units', 'inches', 'Position', [2 10 12 15])
    set(gcf, 'Units', 'inches', 'Position', [2 8 22 6])
    
%     annotation('textbox', [0.175 0.96 0.2 0.03], 'String', 'Rain only',...
%         'FontSize',18,'FontWeight','bold','EdgeColor','k','LineWidth',1,...
%         'HorizontalAlignment','center','VerticalAlignment','middle')
%     annotation('textbox', [0.615 0.96 0.2 0.03], 'String', 'Rain + debris',...
%         'FontSize',18,'FontWeight','bold','EdgeColor','k','LineWidth',1,...
%         'HorizontalAlignment','center','VerticalAlignment','middle')

    
    if plot_save
        print([img_dir 'simr_PPI'], '-dpng')
    end
    %clear ndsim ddsim
    
    
    %-------- This is the more fancy version because I hate myself --------
    
    %%% DO NOT TOUCH THESE LINES THEY ARE EXACTLY WHAT THEY NEED TO BE
    xwhitespace = '                                                                          ';
    ywhitespace = '                                                   ';
    %%% THEY ARE FOR THE X AND Y AXIS LABELS
    %%% I WILL NOT BE TAKING CONSTRUCTIVE CRITICISM OF MY CODING PRACTICES
    
    f = figure(20); clf
    set(gcf, 'Units', 'inches', 'Position', [2 8 16 7])
    p1 = uipanel('Position', [0.03 0.03 0.45 0.97], 'Units', 'normalized', 'BorderType', 'none',...
        'Title', 'Rain Simulation', 'TitlePosition', 'centertop', 'FontSize', 16);
    p2 = uipanel('Position', [0.52 0.03 0.45 0.97], 'Units', 'normalized', 'BorderType', 'none',...
        'Title', 'Realistic Simulation', 'TitlePosition', 'centertop', 'FontSize', 16);
    
    % subtightplot is a user-defined function:
    % https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
    subsubplot = @(m,n,p,pr) subtightplot(m,n,p, [0.06 0.04], [0.05 0.05], [0.07 0.02], 'Parent', pr);
    
    h = subsubplot(2,2,1,p1);
    pcolor(ndsim.xx/1000, ndsim.yy/1000, squeeze(ndsim.zh))
    caxis([0 80])
    colormap(h, blib('zmap'))
    c = colorbar;
    c.Label.String = 'Z (dBZ)';
    c.Label.FontSize = 12;
    c.Ticks = 0:20:80;
    c.Location = 'northoutside';
    xticklabels([])
    xlabel(' ', 'FontSize', 12)
    ylabel(' ', 'FontSize', 12)
    shading flat
    axis tight
    
    h(2) = subsubplot(2,2,2,p1);
    pcolor(ndsim.xx/1000, ndsim.yy/1000, squeeze(ndsim.vr))
    caxis(vr_clims)
    colormap(h(2), blib('rgmap2'))
    c = colorbar;
    c.Label.String = 'v_r (m s^{-1})';
    c.Label.FontSize = 12;
    c.Ticks = -100:25:100;
    c.Location = 'northoutside';
    xticklabels([])
    yticklabels([])
    xlabel(' ', 'FontSize', 12)
    ylabel(' ', 'FontSize', 12)
    shading flat
    axis tight
    h(2).Position = [h(2).Position(1:2), h(1).Position(3:4)];
    
    h(3) = subsubplot(2,2,3,p1);
    pcolor(ndsim.xx/1000, ndsim.yy/1000, squeeze(ndsim.zdr))
    caxis(zdr_clims)
    colormap(h(3), blib('nwsdmap'))
    c = colorbar;
    c.Ticks = -5:2:5;
    c.Location = 'southoutside';
    c.Label.String  = 'Z_{DR} (dB)';
    c.Label.FontSize = 12;
    xlabel([xwhitespace 'Zonal distance (km)'], 'FontSize', 12)
    ylabel([ywhitespace 'Meridional distance (km)'], 'FontSize', 12)
    shading flat
    h(3).Position = [h(3).Position(1:2), h(1).Position(3:4)];
    axis tight
    
    h(4) = subsubplot(2,2,4,p1);
    pcolor(ndsim.xx/1000, ndsim.yy/1000, squeeze(real(ndsim.rhohv)))
    caxis(rho_clims)
    colormap(h(4), blib('nwsrmap'))
    c = colorbar;
    c.Label.String = '\rho_{HV}';
    c.Label.FontSize = 12;
    c.Ticks = 0.4:0.1:1;
    c.Location = 'southoutside';
    yticklabels([])
    xlabel(' ', 'FontSize', 12)
    ylabel(' ', 'FontSize', 12)
    shading flat
    h(4).Position = [h(4).Position(1:2), h(1).Position(3:4)];
    axis tight
    
    h(5) = subsubplot(2,2,1,p2);
    pcolor(ddsim.xx/1000, ddsim.yy/1000, squeeze(ddsim.zh))
    caxis([0 80])
    colormap(h(5), blib('zmap'))
    c = colorbar;
    c.Label.String = 'Z (dBZ)';
    c.Label.FontSize = 12;
    c.Ticks = 0:20:80;
    c.Location = 'northoutside';
    xticklabels([])
    xlabel(' ', 'FontSize', 12)
    ylabel(' ', 'FontSize', 12)
    shading flat
    axis tight
    h(5).Position = [h(5).Position(1:2), h(1).Position(3:4)];
    
    h(6) = subsubplot(2,2,2,p2);
    pcolor(ddsim.xx/1000, ddsim.yy/1000, squeeze(ddsim.vr))
    caxis(vr_clims)
    colormap(h(6), blib('rgmap2'))
    c = colorbar;
    c.Label.String = 'v_r (m s^{-1})';
    c.Label.FontSize = 12;
    c.Ticks = -100:25:100;
    c.Location = 'northoutside';
    xticklabels([])
    yticklabels([])
    xlabel(' ', 'FontSize', 12)
    ylabel(' ', 'FontSize', 12)
    shading flat
    axis tight
    h(6).Position = [h(6).Position(1:2), h(5).Position(3:4)];
    
    h(7) = subsubplot(2,2,3,p2);
    pcolor(ddsim.xx/1000, ddsim.yy/1000, squeeze(ddsim.zdr))
    caxis(zdr_clims)
    colormap(h(7), blib('nwsdmap'))
    c = colorbar;
    c.Label.String = 'Z_{DR} (dB)';
    c.Label.FontSize = 12;
    c.Ticks = -5:2:5;
    c.Location = 'southoutside';
    xlabel([xwhitespace 'Zonal distance (km)'], 'FontSize', 12)
    ylabel([ywhitespace 'Meridional distance (km)'], 'FontSize', 12)
    shading flat
    axis tight
    h(7).Position = [h(7).Position(1:2), h(5).Position(3:4)];
    
    h(8) = subsubplot(2,2,4,p2);
    pcolor(ddsim.xx/1000, ddsim.yy/1000, squeeze(real(ddsim.rhohv)))
    caxis(rho_clims)
    colormap(h(8), blib('nwsrmap'))
    c = colorbar;
    c.Label.String = '\rho_{HV}';
    c.Label.FontSize = 12;
    c.Ticks = 0.4:0.1:1;
    c.Location = 'southoutside';
    yticklabels([])
    xlabel(' ', 'FontSize', 12)
    xlabel(' ', 'FontSize', 12)
    shading flat
    axis tight
    h(8).Position = [h(8).Position(1:2), h(5).Position(3:4)];
    
    set(h, 'DataAspect', [1 1 1])
    
    if plot_save
        print([img_dir 'simr_PPI_v2'], '-dpng')
    end
    clear ndsim ddsim
end


%%% SECTION 4 %%%

%parameters for plotting spectra
dt = 3; %debris type
cn = 1; %debris concentration INDEX (from [10000, 100000])
el = 2; %elevation INDEX (from [0.5, 2.0, 5.0])
ri = 15; %range index
azi = 19; %azimuth index

% Sample spectra for rain + debris
if simr_plot(2)
    caxr = [0.5 0.5 0.8];
    
    % Realistic
    psd_rl = dd(dt).conc(cn).el(el).PSD.old;
    szdr_rl = dd(dt).conc(cn).el(el).szdr;
    sphv_rl = dd(dt).conc(cn).el(el).sphv;
    M = size(psd_rl, 1);
    d = nuttallwin(M); % data window
    
    %clear iqh_nd iqh_dd
    
    % Debris only
    sim_dir = ['~/Documents/sims/' LES '/' sim_date '/'];
    if ~exist('iqh_dd','var')
        dd_dir = [sim_dir 'debris' num2str(dt) '/'];
        dd_fname = ['sim-PPI' num2str(els(el),'%.1f') '-TC-d' num2str(dt) 'n' num2str(cns(cn)) '.mat'];
        load([dd_dir dd_fname],'iqh');
        load([dd_dir dd_fname],'iqv');
        iqh_dd = permute(iqh, [3 1 2]);
        iqv_dd = permute(iqv, [3 1 2]);
        psd_dd = fftshift(abs(fft(iqh_dd.*d, M, 1)).^2, 1);
        [szdr_dd, sphv_dd, ~, ~] = dpsd_calc(iqh_dd, iqv_dd, d, M, 20, 1);
    end
    
    % Rain only
    if ~exist('iqh_nd','var')
        nd_dir = [sim_dir 'nodebris/'];
        nd_fname = ['sim-PPI' num2str(els(el),'%.1f') '-U-nodebris.mat'];
        load([nd_dir nd_fname],'iqh');
        load([nd_dir nd_fname],'iqv');
        iqh_nd = permute(iqh, [3 1 2]);
        iqv_nd = permute(iqv, [3 1 2]);
        psd_nd = fftshift(abs(fft(iqh_nd.*d, M, 1)).^2, 1);
        [szdr_nd, sphv_nd, ~, ~] = dpsd_calc(iqh_nd, iqv_nd, d, M, 20, 1);
    end
    
    
    figure(2)
    clf
    
    ax = subplot(3,2,1);
    yyaxis left
    plot(vvx, szdr_nd(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r', 'FontSize', 14)
    ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va va])
    ylim([-15 15])
    yyaxis right
    hs = semilogy(vvx, abs(psd_nd(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(a) sZ_{DR} (rain)', 'FontSize', 16)
    legend('sZ_{DR}', 'PSD', 'Location', 'northwest')
    
    ax = subplot(3,2,2);
    yyaxis left
    plot(vvx, sphv_nd(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r', 'FontSize', 14)
    ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va va])
    ylim([0 1])
    yticks(0:0.2:1)
    yyaxis right
    hs = semilogy(vvx, abs(psd_nd(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(b) s\rho_{HV} (rain)', 'FontSize', 16)
    legend('s\rho_{HV}', 'PSD', 'Location', 'northwest')
    
    ax = subplot(3,2,3);
    yyaxis left
    plot(vvx, szdr_dd(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r', 'FontSize', 14)
    ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va va])
    ylim([-15 15])
    yyaxis right
    hs = semilogy(vvx, abs(psd_dd(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(c) sZ_{DR} (debris)', 'FontSize', 16)
    legend('sZ_{DR}', 'PSD', 'Location', 'northwest')
    
    ax = subplot(3,2,4);
    yyaxis left
    plot(vvx, sphv_dd(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r', 'FontSize', 14)
    ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va va])
    ylim([0 1])
    yticks(0:0.2:1)
    yyaxis right
    hs = semilogy(vvx, abs(psd_dd(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(d) s\rho_{HV} (debris)', 'FontSize', 16)
    legend('s\rho_{HV}', 'PSD', 'Location', 'northwest')
    
    ax = subplot(3,2,5);
    yyaxis left
    plot(vvx, szdr_rl(:,ri,azi), '-k', 'LineWidth', 1.7)
    xlabel('v_r (m s^{-1})', 'FontSize', 16)
    ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va va])
    ylim([-15 15])
    yyaxis right
    hs = semilogy(vvx, abs(psd_rl(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(e) sZ_{DR} (rain & debris)', 'FontSize', 16)
    legend('sZ_{DR}', 'PSD', 'Location', 'northwest')
    
    ax = subplot(3,2,6);
    yyaxis left
    plot(vvx, sphv_rl(:,ri,azi), '-k', 'LineWidth', 1.7)
    xlabel('v_r (m s^{-1})', 'FontSize', 16)
    ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va va])
    ylim([0 1])
    yticks(0:0.2:1)
    yyaxis right
    hs = semilogy(vvx, abs(psd_rl(:,ri,azi)), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    ylabel('S_H', 'FontSize', 14)
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = caxr;
    title('(f) s\rho_{HV} (rain & debris)', 'FontSize', 16)
    legend('s\rho_{HV}', 'PSD', 'Location', 'northwest')
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 18 12])
    
    if plot_save
        print([img_dir 'simr_DPSD'], '-dpng')
    end
end


    
% Example DPSDs & DCA classification
if simr_plot(3)
    caxr = [0.5 0.5 0.8];
    
    figure(3)
    clf
    
    ax = subplot(3,1,1);
    %yyaxis left
    plot(vvx, dd(dt).conc(cn).el(el).szdr(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r (m s^{-1})', 'FontSize', 16)
    %ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va va])
    ylim([-15 15])
    title('(a) sZ_{DR}','FontSize', 16)
%     yyaxis right
%     hs = semilogy(vvx, abs(dd(dt).conc(cn).el(el).PSD.old(:,11,jj)), 'LineWidth', 2);
%     hs.Color = caxr;
%     hs.LineStyle = ':';
%     ylabel('S_H', 'FontSize', 14)
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = caxr;
%     legend('sZ_{DR}', 'PSD', 'Location', 'northeast')
%     
    ax = subplot(3,1,2);
    %yyaxis left
    plot(vvx, dd(dt).conc(cn).el(el).sphv(:,ri,azi), '-k', 'LineWidth', 1.7)
    %xlabel('v_r (m s^{-1})', 'FontSize', 16)
    %ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va va])
    ylim([0 1])
    yticks(0:0.2:1)
    title('(b) s\rho_{HV}','FontSize', 16)
%     yyaxis right
%     hs = semilogy(vvx, abs(dd(dt).conc(cn).el(el).PSD.old(:,11,jj)), 'LineWidth', 2);
%     hs.Color = caxr;
%     hs.LineStyle = ':';
%     ylabel('S_H', 'FontSize', 14)
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = caxr;
%     legend('s\rho_{HV}', 'PSD', 'Location', 'northeast')
    
    ax = subplot(3,1,3);
    %yyaxis left
    plot(vvx, dd(dt).conc(cn).el(el).obj_class(:,ri,azi), '-k', 'LineWidth', 1.7)
    xlabel('v_r (m s^{-1})', 'FontSize', 16)
    %ylabel('DCA', 'FontSize', 14)
    xlim([-va va])
    ylim([-0.2 1.2])
    yticks(0:0.2:1)
    title('(c) DCA classification','FontSize', 16)
    text(-95, 1.1, 'Debris', 'FontSize', 16)
    text(-95, 0, 'Rain', 'FontSize', 16)
%     yyaxis right
%     hs = semilogy(vvx, abs(dd(dt).conc(cn).el(el).PSD.old(:,11,jj)), 'LineWidth', 2);
%     hs.Color = caxr;
%     hs.LineStyle = ':';
%     ylabel('S_H', 'FontSize', 14)
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = caxr;
%     legend('DCA class', 'PSD', 'Location', 'northeast')
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 10 10])
    
    if plot_save
        print([img_dir 'simr_DCA'], '-dpng')
    end
end


%%% SECTION 5 %%%


% Example corrected spectra for 3 methods
% Look at figs 11/12/13, 17
if simr_plot(4)
    caxr = [0.6 0.6 0.6];
    
    %Filtered spectra with uncorrected, corrected, and truth velocities
    figure(4)
    clf
    
    subplot(3,1,1) %DCA correction
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_truth(ri,azi), [1e0 1e15], ':k', 'LineWidth', 2)
    hold on
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_old(ri,azi), [1e0 1e15], ':b', 'LineWidth', 2)
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_new.dca(ri,azi), [1e0 1e15], ':r', 'LineWidth', 2)
    hs = semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.old(:,ri,azi))), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.dca(:,ri,azi))), 'k', 'LineWidth', 2)
    hold off
    xlabel('v_r (m s^{-1})', 'FontSize', 14)
    ylabel('S_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e15])
    title('(a) DCA-based correction', 'FontSize', 14)
    legend(['v_{truth} = ' num2str(dd(dt).conc(cn).el(el).vr_truth(ri,azi),'%.1f')],...
        ['v_{old} = ' num2str(dd(dt).conc(cn).el(el).vr_old(ri,azi),'%.1f')],...
        ['v_{new} = ' num2str(dd(dt).conc(cn).el(el).vr_new.dca(ri,azi),'%.1f')],...
        'Original PSD', 'Filtered PSD', 'Location', 'northwest')
    
    subplot(3,1,2) %Variance correction
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_truth(ri,azi), [1e0 1e15], ':k', 'LineWidth', 2)
    hold on
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_old(ri,azi), [1e0 1e15], ':b', 'LineWidth', 2)
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_new.var(ri,azi), [1e0 1e15], ':r', 'LineWidth', 2)
    hs = semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.old(:,ri,azi))), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.var(:,ri,azi))), 'k', 'LineWidth', 1.5)
    hold off
    xlabel('v_r (m s^{-1})', 'FontSize', 14)
    ylabel('S_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e15])
    title('(b) Variance-based correction', 'FontSize', 14)
    legend(['v_{truth} = ' num2str(dd(dt).conc(cn).el(el).vr_truth(ri,azi),'%.1f')],...
        ['v_{old} = ' num2str(dd(dt).conc(cn).el(el).vr_old(ri,azi),'%.1f')],...
        ['v_{new} = ' num2str(dd(dt).conc(cn).el(el).vr_new.var(ri,azi),'%.1f')],...
        'Original PSD', 'Filtered PSD', 'Location', 'northwest')
    
    subplot(3,1,3) %Aggregation correction
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_truth(ri,azi), [1e0 1e15], ':k', 'LineWidth', 2)
    hold on
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_new.agg(ri,azi), [1e0 1e15], ':r', 'LineWidth', 2)
    semilogy(ones(1,2)*dd(dt).conc(cn).el(el).vr_old(ri,azi), [1e0 1e15], ':b', 'LineWidth', 2)
    hs = semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.old(:,ri,azi))), 'LineWidth', 2);
    hs.Color = caxr;
    hs.LineStyle = ':';
    semilogy(vvx(:), squeeze(abs(dd(dt).conc(cn).el(el).PSD.agg(:,ri,azi))), 'k', 'LineWidth', 1.5)
    hold off
    xlabel('v_r (m s^{-1})', 'FontSize', 14)
    ylabel('S_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e15])
    title('(c) Aggregation-based correction', 'FontSize', 14)
    legend(['v_{truth} = ' num2str(dd(dt).conc(cn).el(el).vr_truth(ri,azi),'%.1f')],...
        ['v_{old} = ' num2str(dd(dt).conc(cn).el(el).vr_old(ri,azi),'%.1f')],...
        ['v_{new} = ' num2str(dd(dt).conc(cn).el(el).vr_new.agg(ri,azi),'%.1f')],...
        'Original PSD', 'Filtered PSD', 'Location', 'northwest')
    
    set(gcf, 'Units', 'inches', 'Position', [1 8 12 12])
    
    if plot_save
        print([img_dir 'simr_newPSD'], '-dpng')
    end
end


% Correction PPIs and histograms
% left corrected velocity PPIs, middle bias PPIs, right bias histograms
% top row original velocity + bias
% Look at fig 8
if simr_plot(5)
    el1 = 2; %2.0
    el2 = 3; %5.0
    
    xx1 = dd(dt).conc(1).el(el1).params.xx(:,:,1);
    yy1 = dd(dt).conc(1).el(el1).params.yy(:,:,1);
    xx2 = dd(dt).conc(1).el(el2).params.xx(:,:,1);
    yy2 = dd(dt).conc(1).el(el2).params.yy(:,:,1);
    
    histc = [0.7 0.7 0.7];
    
    figure(5)
    clf
    
    % UNCORRECTED
    h = subplot(4,3,1);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).vr_old(:,:,1))
    caxis([-1 1] * va)
    colormap(h, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    h(2) = subplot(4,3,2);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).dv.uncorr.data(:,:,1))
    caxis([-1 1] * va)
    colormap(h(2), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    title('Pre-correction (10,000 wood boards, 2.0{\circ})', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(4,3,3)
    histogram(dd(dt).conc(1).el(el1).dv.uncorr.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.mean),'%.1f')], 'FontSize', 12)
    text(-100, 290, ['Q_1 = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.prc25),'%.1f')], 'FontSize', 12)
    text(-100, 250, ['Q_3 = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.prc75),'%.1f')], 'FontSize', 12)
    
    
    % DEBRIS 3, CONCENTRATION 10,000, ELEVATION 2.0
    h(3) = subplot(4,3,4);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).vr_new.agg(:,:,1))
    caxis([-1 1] * va)
    colormap(h(3), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    h(4) = subplot(4,3,5);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).dv.agg.data(:,:,1))
    caxis([-1 1] * va)
    colormap(h(4), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    title('10,000 wood boards, 2.0{\circ}', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(4,3,6)
    histogram(dd(dt).conc(1).el(el1).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-100, 290, ['Q_1 = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-100, 250, ['Q_3 = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    
    % DEBRIS 3, CONCENTRATION 100,000, ELEVATION 2.0
    h(5) = subplot(4,3,7);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).vr_new.agg(:,:,1))
    caxis([-1 1] * va)
    colormap(h(5), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    h(6) = subplot(4,3,8);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).dv.agg.data(:,:,1))
    caxis([-1 1] * va)
    colormap(h(6), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    title('100,000 wood boards, 2.0{\circ}', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(4,3,9)
    histogram(dd(dt).conc(2).el(el1).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-100, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-100, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    
    % DEBRIS 3, CONCENTRATION 100,000, ELEVATION 5.0
    h(7) = subplot(4,3,10);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).vr_new.agg(:,:,1))
    caxis([-1 1] * va)
    colormap(h(7), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    h(8) = subplot(4,3,11);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).dv.agg.data(:,:,1))
    caxis([-1 1] * va)
    colormap(h(8), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    title('100,000 wood boards, 5.0{\circ}', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(4,3,12)
    histogram(dd(dt).conc(2).el(el2).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-100, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-100, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 16 13])
    %set(h, 'DataAspect', [1 1 1])
    
    if plot_save
        print([img_dir 'simr_biasPPI'], '-dpng')
    end
    
    
    
    
    blim = 100;
    
    figure(21); clf % 10,000, 2.0 DEG
    
    ha = subplot(2,3,1);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).vr_old(:,:,1))
    caxis([-va-1 va+1])
    colormap(ha, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    ha(2) = subplot(2,3,2);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).dv.uncorr.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(ha(2), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(2,3,3)
    histogram(dd(dt).conc(1).el(el1).dv.uncorr.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(1).el(el1).dv.uncorr.prc75),'%.1f')], 'FontSize', 12)
    
    
    ha(3) = subplot(2,3,4);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).vr_new.agg(:,:,1))
    caxis([-va-1 va+1])
    colormap(ha(3), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Corrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    ha(4) = subplot(2,3,5);
    pcolor(xx1, yy1, dd(dt).conc(1).el(el1).dv.agg.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(ha(4), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Corrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(2,3,6)
    histogram(dd(dt).conc(1).el(el1).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(1).el(el1).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(1).el(el1).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 20 8]) %also good: [2 10 14 6]
    sgtitle('10,000 wood boards, 2{\circ} elevation', 'FontSize', 16, 'FontWeight', 'bold')
    %set(ha, 'DataAspect', [1 1 1])
    
    
    
    figure(22); clf % 100,000, 2.0 DEG
    
    hb = subplot(2,3,1);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).vr_old(:,:,1))
    caxis([-va-1 va+1])
    colormap(hb, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    hb(2) = subplot(2,3,2);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).dv.uncorr.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(hb(2), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(2,3,3)
    histogram(dd(dt).conc(2).el(el1).dv.uncorr.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el1).dv.uncorr.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el1).dv.uncorr.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el1).dv.uncorr.prc75),'%.1f')], 'FontSize', 12)
    
    
    hb(3) = subplot(2,3,4);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).vr_new.agg(:,:,1))
    caxis([-va-1 va+1])
    colormap(hb(3), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Corrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    hb(4) = subplot(2,3,5);
    pcolor(xx1, yy1, dd(dt).conc(2).el(el1).dv.agg.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(hb(4), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 12;
    c.Ticks = -100:25:100;
    title('Corrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    subplot(2,3,6)
    histogram(dd(dt).conc(2).el(el1).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el1).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 12)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el1).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 20 8]) %also good: [2 10 14 6]
    sgtitle('100,000 wood boards, 2{\circ} elevation', 'FontSize', 16, 'FontWeight', 'bold')
    
    
    
    figure(23); clf % 100,000, 5.0 DEG
    
    hc = subplot(2,3,1);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).vr_old(:,:,1))
    caxis([-va-1 va+1])
    colormap(hc, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 14;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    hc(2) = subplot(2,3,2);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).dv.uncorr.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(hc(2), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 14;
    c.Ticks = -100:25:100;
    title('Uncorrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    subplot(2,3,3)
    histogram(dd(dt).conc(2).el(el2).dv.uncorr.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 14)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el2).dv.uncorr.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el2).dv.uncorr.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el2).dv.uncorr.prc75),'%.1f')], 'FontSize', 12)
    
    
    hc(3) = subplot(2,3,4);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).vr_new.agg(:,:,1))
    caxis([-va-1 va+1])
    colormap(hc(3), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 14;
    c.Ticks = -100:25:100;
    title('Corrected v_r', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    hc(4) = subplot(2,3,5);
    pcolor(xx2, yy2, dd(dt).conc(2).el(el2).dv.agg.data(:,:,1))
    caxis([-va-1 va+1])
    colormap(hc(4), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Title.String = 'm s^{-1}';
    c.Title.FontSize = 14;
    c.Ticks = -100:25:100;
    title('Corrected v_r bias', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    subplot(2,3,6)
    histogram(dd(dt).conc(2).el(el2).dv.agg.data, -va:5:va, 'FaceColor', histc)
    hold on
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.mean, [0 350], 'k', 'LineWidth', 1.5)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.prc25, [0 350], '--k', 'LineWidth', 1)
    plot(ones(1,2)*dd(dt).conc(2).el(el2).dv.agg.prc75, [0 350], '--k', 'LineWidth', 1)
    hold off
    xlabel('v_r bias (m s^{-1})', 'FontSize', 14)
    xlim([-blim blim])
    ylim([0 350])
    text(-blim+5, 330, ['\mu = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.mean),'%.1f')], 'FontSize', 12)
    text(-blim+5, 290, ['Q_1 = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.prc25),'%.1f')], 'FontSize', 12)
    text(-blim+5, 250, ['Q_3 = ' num2str(double(dd(dt).conc(2).el(el2).dv.agg.prc75),'%.1f')], 'FontSize', 12)
    
    set(gcf, 'Units', 'inches', 'Position', [2 10 20 8]) %also good: [2 10 14 6]
    sgtitle('100,000 wood boards, 5{\circ} elevation', 'FontSize', 16, 'FontWeight', 'bold')
    
    
    
end


% Scatter plots with polarimetric vars (rhohv) and bias/correction magnitude
% top row corrected bias, bottom row correction magnitude
if simr_plot(6)
    
    phv = [reshape(dd(dt).conc(1).el(2).params.phv,[],1);
                    reshape(dd(dt).conc(2).el(2).params.phv,[],1);
                    reshape(dd(dt).conc(1).el(3).params.phv,[],1);
                    reshape(dd(dt).conc(2).el(3).params.phv,[],1)];
    vcorr = [reshape(dd(dt).conc(1).el(2).vcorr.agg,[],1);
                    reshape(dd(dt).conc(2).el(2).vcorr.agg,[],1);
                    reshape(dd(dt).conc(1).el(3).vcorr.agg,[],1);
                    reshape(dd(dt).conc(2).el(3).vcorr.agg,[],1)];
    vbias = [reshape(dd(dt).conc(1).el(2).dv.agg.data,[],1);
                    reshape(dd(dt).conc(2).el(2).dv.agg.data,[],1);
                    reshape(dd(dt).conc(1).el(3).dv.agg.data,[],1);
                    reshape(dd(dt).conc(2).el(3).dv.agg.data,[],1)];
    ucbias = [reshape(dd(dt).conc(1).el(2).dv.uncorr.data,[],1);
                    reshape(dd(dt).conc(2).el(2).dv.uncorr.data,[],1);
                    reshape(dd(dt).conc(1).el(3).dv.uncorr.data,[],1);
                    reshape(dd(dt).conc(2).el(3).dv.uncorr.data,[],1)];
    vtruth = [reshape(dd(dt).conc(1).el(2).vr_truth,[],1);
                    reshape(dd(dt).conc(2).el(2).vr_truth,[],1);
                    reshape(dd(dt).conc(1).el(3).vr_truth,[],1);
                    reshape(dd(dt).conc(2).el(3).vr_truth,[],1)];
    
    %stuff
    rsquared = @(y,yfit) 1 - sum((y-yfit).^2)/sum((y-mean(y)).^2);
    %rsquared = @(y,normr) 1 - (normr/norm(y-mean(y)))^2;
    
    cffs1 = polyfit(phv, abs(vbias), 1);
    rhofit = phv;
    vbfit1 = polyval(cffs1, rhofit);
    rsq1 = rsquared(abs(vbias), vbfit1(:));
    
    cffs2 = polyfit(phv, abs(vcorr), 1);
    vcfit1 = polyval(cffs2, rhofit);
    rsq2 = rsquared(abs(vcorr), vcfit1(:));
    
    cffs3 = polyfit(vtruth, vbias, 1);
    vtfit = vtruth;
    vbfit2 = polyval(cffs3, vtfit);
    rsq3 = rsquared(vbias, vbfit2(:));
    
    cffs4 = polyfit(vtruth, vcorr, 1);
    vcfit2 = polyval(cffs4, vtfit);
    rsq4 = rsquared(vcorr, vcfit2(:));
    
    cffs5 = polyfit(ucbias, vbias, 1);
    vbfit3 = polyval(cffs5, ucbias);
    rsq5 = rsquared(vbias, vbfit3(:));
    
    cffs6 = polyfit(ucbias, vcorr, 1);
    %ubfit = ucbias;
    vcfit3 = polyval(cffs6, ucbias);
    rsq6 = rsquared(vcorr, vcfit3(:));
    
    
    
    
    figure(6)
    clf
    
    subplot(3,2,1)
     scatter(phv, abs(vbias), '.')
%    scatter(phv, vbias, '.')
    hold on
    plot(rhofit, vbfit1, 'k', 'LineWidth', 2)
    hold off
    xlabel('rhohv')
    ylabel('v bias')
    ylim([-10 50])
    title(['r^2 = ' num2str(rsq1,'%.2f')])
    
    subplot(3,2,2)
     scatter(phv, abs(vcorr), '.')
 %   scatter(phv, vcorr, '.')
    hold on
    plot(rhofit, vcfit1, 'k', 'LineWidth', 2)
    hold off
    xlabel('rhohv')
    ylabel('v correction')
    ylim([-10 50])
    title(['r^2 = ' num2str(rsq2,'%.2f')])
    
    subplot(3,2,3)
    scatter(vtruth, vbias, '.')
    hold on
    plot(vtfit, vbfit2, 'k', 'LineWidth', 2)
    hold off
    xlabel('truth v_r')
    ylabel('v bias')
    ylim([-50 50])
    title(['r^2 = ' num2str(rsq3,'%.2f')])
    
    subplot(3,2,4)
    scatter(vtruth, vcorr, '.')
    hold on
    plot(vtfit, vcfit2, 'k', 'LineWidth', 2)
    hold off
    xlabel('truth v_r')
    ylabel('v correction')
    ylim([-50 50])
    title(['r^2 = ' num2str(rsq4,'%.2f')])
    
    subplot(3,2,5)
    scatter(ucbias, vbias, '.')
    hold on
    plot(ucbias, vbfit3, 'k', 'LineWidth', 2)
    hold off
    xlabel('original bias')
    ylabel('v bias')
    ylim([-50 50])
    title(['r^2 = ' num2str(rsq5,'%.2f')])
    
    subplot(3,2,6)
    scatter(ucbias, vcorr, '.')
    hold on
    plot(ucbias, vcfit3, 'k', 'LineWidth', 2)
    hold off
    xlabel('original bias')
    ylabel('v correction')
    xlim([-100 100])
    ylim([-50 50])
    title(['r^2 = ' num2str(rsq6,'%.2f')])
    
    if plot_save
        print([img_dir 'simr_biascorr'], '-dpng')
    end
    
end



%% KOUN plots


%%% SECTION 2 %%%

% Cropped PPIs w/ sample spectra gates marked
if koun_plot(1)
    figure(10)
    clf
    
    
    
    if plot_save
        print([img_dir 'koun_PPI'], '-dpng')
    end
end



%%% SECTION 4 %%%

% Example spectra from center + edge of TDS
if koun_plot(2)
    figure(11)
    clf
    
    
    
    if plot_save
        print([img_dir 'koun_DPSD'], '-dpng')
    end
end


% Example DPSDs & DCA classification
if koun_plot(3)
    figure(12)
    clf
    
    
    
    if plot_save
        print([img_dir 'koun_DCA'], '-dpng')
    end
end


% SECTION 6


% Scatter plots with polarimetric vars (rhohv) and correction magnitude
if koun_plot(4)
    figure(13)
    clf
    
    
    
    if plot_save
        print([img_dir 'koun_biascorr'], '-dpng')
    end
end


% Corrected velocity and correction magnitude PPIs
% include aliased and dealiased?
% include bulk measure of radial divergence pre-/post-correction?
if koun_plot(5)
    figure(14)
    clf
    
    
    
    if plot_save
        print([img_dir 'koun_biasPPI'], '-dpng')
    end
end





%%












