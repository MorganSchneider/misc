%This script was used in the final output for the JTech Frequency paper.
%It will take in the stats file from the radar simulation, and make Figures
%6 ? 9, and compute the data for the tables

clear;
close all;
v0=225;
% data_dir=['/Users/alexlyakhov/Documents/LES/stats_traj/traj_test/LES/squareplate/ravg/' num2str(v0)];
% data_dir=['/Users/alexlyakhov/Documents/LES/stats_traj/traj_test/LES/woodsphere/ravg/' num2str(v0)];
va=20;
vmax=100;
corr_ind=5;
addpath data_files

%Batch in sets
deb_scale_lg_set=[0.05 0.05 0.05 0.005 0.0005];
rain_scales_set=[0.001 0.01 0.1 0.1 0.1];
orient_case=1; %1: Completely random, 2: single plane

rcs_flags=[1 1 1 1 1 0 0 0 0 0 1 1 1 1 1];
dry_flags=[1 1 1 1 1 1 1 1 1 1 0 0 0 0 0];
deb_scale_lg=[deb_scale_lg_set deb_scale_lg_set deb_scale_lg_set];
rain_scales=[rain_scales_set rain_scales_set rain_scales_set];

base_corr_flag=0.01;
% deb_scale_large=0.05;
% rain_scale=0.1;
% dry_flag=1;
orient_case=1; %1: random, 2: single plane
% rcs_flag=1; %1 for POA, 0 for T-matrix

% if(dry_flag)
%     dry_title='dry';
% else
%     dry_title='wet';
% end
% 
% if(rcs_flag)
%     em_title='poa';
% else
%     em_title='tm';
% end

flow_type=2;
if(flow_type==1)
    flow_name='twocell';
elseif(flow_type==2)
    flow_name='vb2_equal';
    vfile=225;
elseif(flow_type==3)
    flow_name='twocell_fin';
end
data_dir=['/Users/dbodine/Documents/LES_radar/radar_out/' flow_name '/'];

stats_dir=data_dir;
if(~exist(stats_dir,'dir'))
    mkdir(stats_dir);
end

base_dir='/Users/dbodine/Documents/LES_radar/';

lb=0.107; %m
ls=0.031;

noise_db=1;
s_noise=10^(noise_db/10); b_noise=s_noise;
dwr_thres_db=1;
dws_thres=-10;
dwr_max=lb^4/ls^4*2000;
dwr_thres=10^(dwr_thres_db/10);
noise_floor_db=10^(-100/10);
sigma_deb=1; sigma_rain=1;
sigma1(1,1,1:15)=[ones(1,10)*sigma_deb ones(1,5)*sigma_rain];
% dwr_max=10^(dwr_max_db/10);

make_figs=true;
make_figs_analysis=false;
spectral_flag=false;
conf_flag=false;
new_title=true;
small_size_flag=0;

for case_idx=1:max(size(rcs_flags))
    
    rcs_flag=rcs_flags(case_idx);
    dry_flag=dry_flags(case_idx);
    deb_scale_large=deb_scale_lg(case_idx);
    rain_scale=rain_scales(case_idx);
    if(vfile==150)
        deb_scale_sm=[2.9497    0.4724    0.1584    0.0747    0.0052]*1e5;
    else
        deb_scale_sm=[2.9497    0.4724    0.1584    0.0747    0.0052]*1e5*(225^2/150^2)^3*rain_scale;
    end
    
    if(rcs_flag)
        em_title='poa';
    else
        em_title='tm';
    end

    if(dry_flag)
        dry_title='dry';
    else
        dry_title='wet';
    end
    fprintf(['Working on case: ' num2str(case_idx) '/' num2str(max(size(rcs_flags))) '\n'])
    
    if(new_title)
        if(small_size_flag)
            fname=['small_stats_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.mat'];
        else
            fname=['stats_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.mat'];
        end
    else
        fname=['stats_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(dry_flag) '.mat'];
    end

    cd(data_dir)
        load(fname)
    cd(base_dir)

    lims=[0 700 0 600]; label_x=50; label_y=550;

    % exp_name='square_plate';
    % exp_name='wood_sphere';
    if(min(size(strfind(data_dir,'vb2_equal')))>0)
        fig_dir=['/Users/schneider/Documents/imgs/'];
    else
        fig_dir=['/Users/schneider/Documents/imgs/'];
    end
    dr=rtmp2(2,1)-rtmp2(1,1);
    dz=ztmp2(1,2)-ztmp2(1,1);
    wi=0;
    w_int_s=nan(size(uds));
    w_int_c=nan(size(uds));
    w_int_x=nan(size(uds));
    w_int_ka=nan(size(uds));
    w_int_w=nan(size(uds));
    w_int_t=nan(size(uds));
    w_int_sc=nan(size(uds));
    w_int_cc=nan(size(uds));
    w_int_xc=nan(size(uds));
    w_int_kac=nan(size(uds));
    w_int_wc=nan(size(uds));
    w_int_tc=nan(size(uds));

    radii=[9.5209   19.042  28.563  38.084  47.605  57.126  66.647  76.168  85.688  95.209 0.25 0.5 0.75 1 2];

    if(v0==225)
        ulims=[-40 40];
        vlims=[0 70];
        wlims=[-40 40];
    end

    cmax=0;
    ind=2:size(rtmp2,1);
    u_adjust=cmax*(vdw(ind,:)-vrw(ind,:)).^2./rtmp2(ind,:)./(nanmax(nanmax((vdw(ind,:)-vrw(ind,:)).^2./rtmp2(ind,:))));
    % u_adjust=cmax*(vdw(2:size(rtmp2,1),:)).^2./rtmp2(2:size(rtmp2,1),:)./(nanmax(nanmax((vdw(2:size(rtmp2,1),:)).^2./rtmp2(2:size(rtmp2,1),:))));

    figure(15)
    pcolor(rtmp2(ind,:),ztmp2(ind,:),u_adjust)
    shading flat
    colorbar

    udsc(ind,:)=uds(ind,:)-u_adjust;
    udcc(ind,:)=udc(ind,:)-u_adjust;
    udxc(ind,:)=udx(ind,:)-u_adjust;
    udkac(ind,:)=udka(ind,:)-u_adjust;
    udwc(ind,:)=udw(ind,:)-u_adjust;

    udwt=udw-urw;
    for jdx=1:size(udwc,2)
        for idx=2:size(udwc,1)-2
            if(jdx==1)
                w_int_s(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udsc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udsc(idx-1,jdx)-0);
                w_int_c(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udcc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udcc(idx-1,jdx)-0);
                w_int_x(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udxc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udxc(idx-1,jdx)-0);
                w_int_ka(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udkac(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udkac(idx-1,jdx)-0);
                w_int_w(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udwc(idx-1,jdx)-0);
                w_int_t(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwt(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udwt(idx-1,jdx)-0);
            else
                w_int_s(idx,jdx)=w_int_s(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udsc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udsc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udsc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udsc(idx-1,jdx-1));
                w_int_c(idx,jdx)=w_int_c(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udcc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udcc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udcc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udcc(idx-1,jdx-1));
                w_int_x(idx,jdx)=w_int_x(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udxc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udxc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udxc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udxc(idx-1,jdx-1));
                w_int_ka(idx,jdx)=w_int_ka(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udkac(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udkac(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udkac(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udkac(idx-1,jdx-1));
                w_int_w(idx,jdx)=w_int_w(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udwc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udwc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udwc(idx-1,jdx-1));
                w_int_t(idx,jdx)=w_int_t(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwt(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udwt(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udwt(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udwt(idx-1,jdx-1));
            end
        end
    end

    u_adjust=squeeze(urt_avg_sm(2:size(rtmp2,1),:,corr_ind));
    udsc(ind,:)=uds(ind,:)-u_adjust;
    udcc(ind,:)=udc(ind,:)-u_adjust;
    udxc(ind,:)=udx(ind,:)-u_adjust;
    udkac(ind,:)=udka(ind,:)-u_adjust;
    udwc(ind,:)=udw(ind,:)-u_adjust;

    wi=0;
    udwt=udw-urw;
    for jdx=1:size(udwc,2)
        for idx=2:size(udwc,1)-2
            if(jdx==1)
                w_int_sc(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udsc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udsc(idx-1,jdx)-0);
                w_int_cc(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udcc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udcc(idx-1,jdx)-0);
                w_int_xc(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udxc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udxc(idx-1,jdx)-0);
                w_int_kac(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udkac(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udkac(idx-1,jdx)-0);
                w_int_wc(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwc(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udwc(idx-1,jdx)-0);
                w_int_tc(idx,jdx)=wi-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwt(idx+1,jdx)+0-...
                    rtmp2(idx-1,jdx)*udwt(idx-1,jdx)-0);
            else
                w_int_sc(idx,jdx)=w_int_sc(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udsc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udsc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udsc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udsc(idx-1,jdx-1));
                w_int_cc(idx,jdx)=w_int_cc(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udcc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udcc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udcc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udcc(idx-1,jdx-1));
                w_int_xc(idx,jdx)=w_int_xc(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udxc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udxc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udxc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udxc(idx-1,jdx-1));
                w_int_kac(idx,jdx)=w_int_kac(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udkac(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udkac(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udkac(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udkac(idx-1,jdx-1));
                w_int_wc(idx,jdx)=w_int_wc(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwc(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udwc(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udwc(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udwc(idx-1,jdx-1));
                w_int_tc(idx,jdx)=w_int_tc(idx,jdx-1)-dz/(4*dr)*1/rtmp2(idx,jdx)*(rtmp2(idx+1,jdx)*udwt(idx+1,jdx)+rtmp2(idx+1,jdx-1)*udwt(idx+1,jdx-1)-...
                    rtmp2(idx-1,jdx)*udwt(idx-1,jdx)-rtmp2(idx-1,jdx-1)*udwt(idx-1,jdx-1));
            end
        end
    end

    if(make_figs)
    if(~base_corr_flag)
        w_int_t=zeros(size(w_int_t));
    else
        w_int_t=wds-wrs-w_int_t;
    end
    clims1=[-100 100]; 
    cmap5=feval('boonlib','carbmap',41);
    cmap5=cmap5(2:size(cmap5,1),:);
    rtmp2_adjust=(rtmp2(2,1)-rtmp2(1,1))/2;
    rlab=650; zlab=550;
    % figure(1)
    % feval('boonlib','bsizewin',1,[1000 800])
    % subplot(3,2,1)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_s)
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('S-band w','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(a)','FontSize',12,'BackgroundColor',[1 1 1])
    % subplot(3,2,3)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_x)
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('X-band w','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(b)','FontSize',12,'BackgroundColor',[1 1 1])
    % subplot(3,2,5)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_w)
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('W-band w','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(c)','FontSize',12,'BackgroundColor',[1 1 1])
    % subplot(3,2,2)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_s-(wds-wrs))
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('S-band w error','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(d)','FontSize',12,'BackgroundColor',[1 1 1])
    % subplot(3,2,4)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_x-(wdx-wrx))
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('X-band w error','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(e)','FontSize',12,'BackgroundColor',[1 1 1])
    % subplot(3,2,6)
    % pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_w-(wdw-wrw))
    % shading flat
    % colorbar
    % colormap(cmap5)
    % title('W-band w error','FontSize',14)
    % xlabel('r(m)','FontSize',12)
    % ylabel('z(m)','FontSize',12)
    % axis(lims)
    % caxis(clims1)
    % text(rlab,zlab,'(f)','FontSize',12,'BackgroundColor',[1 1 1])
    figure(1)
    feval('boonlib','bsizewin',1,[1000 800])
    subplot(3,2,1)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_s)
    shading flat
    colorbar
    colormap(cmap5)
    title('S-band w (m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(a)','FontSize',12,'BackgroundColor',[1 1 1])
    subplot(3,2,3)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_x)
    shading flat
    colorbar
    colormap(cmap5)
    title('X-band w (m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(b)','FontSize',12,'BackgroundColor',[1 1 1])
    subplot(3,2,5)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_w)
    shading flat
    colorbar
    colormap(cmap5)
    title('W-band w (m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(c)','FontSize',12,'BackgroundColor',[1 1 1])
    subplot(3,2,2)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_sc)
    shading flat
    colorbar
    colormap(cmap5)
    title('S-band w (Corrected, m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(d)','FontSize',12,'BackgroundColor',[1 1 1])
    subplot(3,2,4)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,w_int_xc)
    shading flat
    colorbar
    colormap(cmap5)
    title('X-band w (Corrected, m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(e)','FontSize',12,'BackgroundColor',[1 1 1])
    subplot(3,2,6)
    pcolor(rtmp2-rtmp2_adjust,ztmp2,wds-wrs)
    shading flat
    colorbar
    colormap(cmap5)
    title('w (m s^{-1})','FontSize',14)
    xlabel('r(m)','FontSize',12)
    ylabel('z(m)','FontSize',12)
    axis(lims)
    caxis(clims1)
    text(rlab,zlab,'(f)','FontSize',12,'BackgroundColor',[1 1 1])

    rmax_w=300;
    [ind]=rtmp2(:,1)<rmax_w;
    rmse_swr=sqrt(nanmean(nanmean((w_int_s(ind,:)-(wds(ind,:)-wrs(ind,:))).^2)));
    rmse_xwr=sqrt(nanmean(nanmean((w_int_x(ind,:)-(wdx(ind,:)-wrx(ind,:))).^2)));
    rmse_wwr=sqrt(nanmean(nanmean((w_int_w(ind,:)-(wdw(ind,:)-wrw(ind,:))).^2)));
    rmse_swrc=sqrt(nanmean(nanmean((w_int_sc(ind,:)-(wds(ind,:)-wrs(ind,:))).^2)));
    rmse_xwrc=sqrt(nanmean(nanmean((w_int_xc(ind,:)-(wdx(ind,:)-wrx(ind,:))).^2)));
    rmse_wwrc=sqrt(nanmean(nanmean((w_int_wc(ind,:)-(wdw(ind,:)-wrw(ind,:))).^2)));

    dz=rtmp2(2,1)-rtmp2(1,1);
    clims=[nanmin(nanmin([wds-wrs-w_int_s wdc-wrc-w_int_c wdx-wrx-w_int_x wdka-wrka-w_int_ka wdw-wrw-w_int_w]))...
        nanmax(nanmax([wds-wrs-w_int_s wdc-wrc-w_int_c wdx-wrx-w_int_x wdka-wrka-w_int_ka wdw-wrw-w_int_w]))];
    figure(2)
    feval('boonlib','bsizewin',2,[1000 800])
    subplot(3,2,1)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_s)
    shading flat
    colorbar
    title('S-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,2)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_c)
    shading flat
    colorbar
    title('C-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,3)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_x)
    shading flat
    colorbar
    title('X-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,4)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_ka)
    shading flat
    colorbar
    title('Ka-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,5)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_w)
    shading flat
    colorbar
    title('W-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,6)
    pcolor(rtmp2-dz/2,ztmp2,wdw-wrw)
    shading flat
    colorbar
    title('Actual w')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)

    figure(211)
    feval('boonlib','bsizewin',211,[1000 800])
    subplot(3,2,1)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_sc)
    shading flat
    colorbar
    title('S-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,2)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_cc)
    shading flat
    colorbar
    title('C-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,3)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_xc)
    shading flat
    colorbar
    title('X-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,4)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_kac)
    shading flat
    colorbar
    title('Ka-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,5)
    pcolor(rtmp2-dz/2,ztmp2,wds-wrs-w_int_wc)
    shading flat
    colorbar
    title('W-band w error')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)
    subplot(3,2,6)
    pcolor(rtmp2-dz/2,ztmp2,wdw-wrw)
    shading flat
    colorbar
    title('Actual w')
    xlabel('r(m)')
    ylabel('z(m)')
    axis(lims)
    caxis(clims)


    if(make_figs_analysis)
    figure(3)
    feval('boonlib','bsizewin',3,[1000 800])
    subplot(4,3,1)
    pcolor(rtmp2,ztmp2,uds-udc)
    shading flat
    colorbar
    title('S-C DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,2)
    pcolor(rtmp2,ztmp2,uds-udx)
    shading flat
    colorbar
    title('S-X DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,3)
    pcolor(rtmp2,ztmp2,uds-udka)
    shading flat
    colorbar
    title('S-Ka DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,4)
    pcolor(rtmp2,ztmp2,uds-udw)
    shading flat
    colorbar
    title('S-W DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,5)
    pcolor(rtmp2,ztmp2,udc-udx)
    shading flat
    colorbar
    title('C-X DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,6)
    pcolor(rtmp2,ztmp2,udc-udka)
    shading flat
    colorbar
    title('C-Ka DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,7)
    pcolor(rtmp2,ztmp2,udc-udw)
    shading flat
    colorbar
    title('C-W DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,8)
    pcolor(rtmp2,ztmp2,udx-udka)
    shading flat
    colorbar
    title('X-Ka DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,9)
    pcolor(rtmp2,ztmp2,udx-udw)
    shading flat
    colorbar
    title('X-W DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,10)
    pcolor(rtmp2,ztmp2,udka-udw)
    shading flat
    colorbar
    title('Ka-W DDU (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')

    figure(4)
    feval('boonlib','bsizewin',4,[1000 800])
    subplot(4,3,1)
    pcolor(rtmp2,ztmp2,vds-vdc)
    shading flat
    colorbar
    title('S-C DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,2)
    pcolor(rtmp2,ztmp2,vds-vdx)
    shading flat
    colorbar
    title('S-X DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,3)
    pcolor(rtmp2,ztmp2,vds-vdka)
    shading flat
    colorbar
    title('S-Ka DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,4)
    pcolor(rtmp2,ztmp2,vds-vdw)
    shading flat
    colorbar
    title('S-W DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,5)
    pcolor(rtmp2,ztmp2,vdc-vdx)
    shading flat
    colorbar
    title('C-X DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,6)
    pcolor(rtmp2,ztmp2,vdc-vdka)
    shading flat
    colorbar
    title('C-Ka DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,7)
    pcolor(rtmp2,ztmp2,vdc-vdw)
    shading flat
    colorbar
    title('C-W DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,8)
    pcolor(rtmp2,ztmp2,vdx-vdka)
    shading flat
    colorbar
    title('X-Ka DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,9)
    pcolor(rtmp2,ztmp2,vdx-vdw)
    shading flat
    colorbar
    title('X-W DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,10)
    pcolor(rtmp2,ztmp2,vdka-vdw)
    shading flat
    colorbar
    title('Ka-W DDV (m s^{-1})')
    xlabel('r(m)')
    ylabel('z(m)')

    figure(5)
    pcolor(rtmp2,ztmp2,urx)
    shading flat
    colorbar
    title('U_{r} X-band')
    xlabel('r(m)')
    ylabel('z(m)')

    figure(6)
    pcolor(rtmp2,ztmp2,vrx)
    shading flat
    colorbar
    title('V_{r} X-band')
    xlabel('r(m)')
    ylabel('z(m)')


    end

    end
    figure(7)
    feval('boonlib','bsizewin',7,[800 600])
    subplot(2,2,1)
    pcolor(rtmp2,ztmp2,udw-urw)
    shading flat
    colorbar
    xlabel('r(m)')
    ylabel('z(m)')
    title('U (m s^{-1})','FontSize',14)
    caxis(ulims)
    subplot(2,2,2)
    pcolor(rtmp2,ztmp2,vdw-vrw)
    shading flat
    colorbar
    caxis(vlims)
    xlabel('r(m)')
    ylabel('z(m)')
    title('V (m s^{-1})','FontSize',14)
    subplot(2,2,3)
    pcolor(rtmp2,ztmp2,wds-wrs)
    shading flat
    colorbar
    caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    title('W (m s^{-1})','FontSize',14)

    manu_flag=true;
    figure(7)
    set(gcf,'PaperPositionMode','auto')
    if(manu_flag)    
        picname=[fig_dir '/mean_vel.pdf'];
        print(picname,'-dpdf');
    else
        picname=[fig_dir '/mean_vel.png'];
        print(picname,'-dpng'); 
    end

    for idx=1:max(size(radii))
       yn=uds_ind==idx; 
       uds_rad(yn)=radii(idx); 
       yn=vds_ind==idx; 
       vds_rad(yn)=radii(idx);
       yn=wds_ind==idx; 
       wds_rad(yn)=radii(idx);
       yn=udc_ind==idx; 
       udc_rad(yn)=radii(idx); 
       yn=vdc_ind==idx; 
       vdc_rad(yn)=radii(idx);
       yn=wdc_ind==idx; 
       wdc_rad(yn)=radii(idx);
       yn=udx_ind==idx; 
       udx_rad(yn)=radii(idx); 
       yn=vdx_ind==idx; 
       vdx_rad(yn)=radii(idx);
       yn=wdx_ind==idx; 
       wdx_rad(yn)=radii(idx);
       yn=udka_ind==idx; 
       udka_rad(yn)=radii(idx); 
       yn=vdka_ind==idx; 
       vdka_rad(yn)=radii(idx);
       yn=wdka_ind==idx; 
       wdka_rad(yn)=radii(idx);
       yn=udw_ind==idx; 
       udw_rad(yn)=radii(idx); 
       yn=vdw_ind==idx; 
       vdw_rad(yn)=radii(idx);
       yn=wdw_ind==idx; 
       wdw_rad(yn)=radii(idx);
    end
    uds_rad=reshape(uds_rad,size(uds_ind,1),size(uds_ind,2));
    udc_rad=reshape(udc_rad,size(uds_ind,1),size(uds_ind,2));
    udx_rad=reshape(udx_rad,size(uds_ind,1),size(uds_ind,2));
    udka_rad=reshape(udka_rad,size(uds_ind,1),size(uds_ind,2));
    udw_rad=reshape(udw_rad,size(udw_ind,1),size(udw_ind,2));
    vds_rad=reshape(vds_rad,size(uds_ind,1),size(uds_ind,2));
    vdc_rad=reshape(vdc_rad,size(uds_ind,1),size(uds_ind,2));
    vdx_rad=reshape(vdx_rad,size(uds_ind,1),size(uds_ind,2));
    vdka_rad=reshape(vdka_rad,size(uds_ind,1),size(uds_ind,2));
    vdw_rad=reshape(vdw_rad,size(udw_ind,1),size(udw_ind,2));
    wds_rad=reshape(wds_rad,size(uds_ind,1),size(uds_ind,2));
    wdc_rad=reshape(wdc_rad,size(uds_ind,1),size(uds_ind,2));
    wdx_rad=reshape(wdx_rad,size(uds_ind,1),size(uds_ind,2));
    wdka_rad=reshape(wdka_rad,size(uds_ind,1),size(uds_ind,2));
    wdw_rad=reshape(wdw_rad,size(udw_ind,1),size(udw_ind,2));

    uds_radf=uds_rad;
    udc_radf=udc_rad;
    udx_radf=udx_rad;
    udka_radf=udka_rad;
    udw_radf=udw_rad;
    %Apply a median filter to reduce noisiness
    for idx=2:size(uds_rad,1)-1
        for jdx=2:size(uds_rad,2)-1
            tmp_array=uds_rad;
            tmp=[tmp_array(idx-1:idx+1,jdx);tmp_array(idx,jdx-1);tmp_array(idx,jdx+1)];
            uds_radf(idx,jdx)=median(tmp);
            tmp_array=udc_rad;
            tmp=[tmp_array(idx-1:idx+1,jdx);tmp_array(idx,jdx-1);tmp_array(idx,jdx+1)];
            udc_radf(idx,jdx)=median(tmp);
            tmp_array=udx_rad;
            tmp=[tmp_array(idx-1:idx+1,jdx);tmp_array(idx,jdx-1);tmp_array(idx,jdx+1)];
            udx_radf(idx,jdx)=median(tmp);
            tmp_array=udka_rad;
            tmp=[tmp_array(idx-1:idx+1,jdx);tmp_array(idx,jdx-1);tmp_array(idx,jdx+1)];
            udka_radf(idx,jdx)=median(tmp);
            tmp_array=udw_rad;
            tmp=[tmp_array(idx-1:idx+1,jdx);tmp_array(idx,jdx-1);tmp_array(idx,jdx+1)];
            udw_radf(idx,jdx)=median(tmp);
        end
    end

    %Handle corner points
    uds_radf(1,1)=median([uds_rad(1,1); uds_rad(2,1); uds_rad(1,2)]);
    udc_radf(1,1)=median([udc_rad(1,1); udc_rad(2,1); udc_rad(1,2)]);
    udx_radf(1,1)=median([udx_rad(1,1); udx_rad(2,1); udx_rad(1,2)]);
    udka_radf(1,1)=median([udka_rad(1,1); udka_rad(2,1); udka_rad(1,2)]);
    udw_radf(1,1)=median([udw_rad(1,1); udw_rad(2,1); udw_rad(1,2)]);

    uds_rad=uds_radf;
    udc_rad=udc_radf;
    udx_rad=udx_radf;
    udka_rad=udka_radf;
    udw_rad=udw_radf;

    if(make_figs)
    cmap1=parula(nanmax(round(radii/10)*10));
    % cmap1(1:2,:)=repmat([1 1 1],[2 1]);
    figure(7)
    feval('boonlib','bsizewin',7,[800 600])
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,uds_rad)
    shading flat
    colorbar
    colormap(cmap1)
    caxis([0 nanmax(round(radii/10)*10)])

    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,2)
    pcolor(rtmp2,ztmp2,udc_rad)
    shading flat
    colorbar
    colormap(cmap1)
    caxis([0 nanmax(round(radii/10)*10)])
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,udx_rad)
    shading flat
    colorbar
    colormap(cmap1)
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,4)
    pcolor(rtmp2,ztmp2,udka_rad)
    shading flat
    colorbar
    colormap(cmap1)
    caxis([0 nanmax(round(radii/10)*10)])
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,udw_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    colormap(cmap1)
    caxis([0 nanmax(round(radii/10)*10)])

    figure(8)
    feval('boonlib','bsizewin',8,[800 600])
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,vds_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,2)
    pcolor(rtmp2,ztmp2,vdc_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,vdx_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,4)
    pcolor(rtmp2,ztmp2,vdka_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,vdw_rad)
    shading flat
    colorbar
    % caxis(wlims)
    xlabel('r(m)')
    ylabel('z(m)')

    end

    zhhs_smt=nansum(zhhs_sm,3);
    zhhc_smt=nansum(zhhc_sm,3);
    zhhx_smt=nansum(zhhx_sm,3);
    zhhka_smt=nansum(zhhka_sm,3);
    zhhw_smt=nansum(zhhw_sm,3);
    figure(15)
    feval('boonlib','bsizewin',15,[1000 800])
    subplot(4,3,1)
    pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_c_tot+zhhc_smt))
    shading flat
    colorbar
    title('S-C Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,2)
    pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_x_tot+zhhx_smt))
    shading flat
    colorbar
    title('S-X Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,3)
    pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_ka_tot+zhhka_smt))
    shading flat
    colorbar
    title('S-Ka Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,4)
    pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_w_tot+zhhw_smt))
    shading flat
    colorbar
    title('S-W Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,5)
    pcolor(rtmp2,ztmp2,10*log10(zhh_c_tot+zhhc_smt)-10*log10(zhh_x_tot+zhhx_smt))
    shading flat
    colorbar
    title('C-X Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,6)
    pcolor(rtmp2,ztmp2,10*log10(zhh_c_tot+zhhc_smt)-10*log10(zhh_ka_tot+zhhka_smt))
    shading flat
    colorbar
    title('C-Ka Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,7)
    pcolor(rtmp2,ztmp2,10*log10(zhh_c_tot+zhhc_smt)-10*log10(zhh_w_tot+zhhw_smt))
    shading flat
    colorbar
    title('C-W Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,8)
    pcolor(rtmp2,ztmp2,10*log10(zhh_x_tot+zhhx_smt)-10*log10(zhh_ka_tot+zhhka_smt))
    shading flat
    colorbar
    title('X-Ka Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,9)
    pcolor(rtmp2,ztmp2,10*log10(zhh_x_tot+zhhx_smt)-10*log10(zhh_w_tot+zhhw_smt))
    shading flat
    colorbar
    title('X-W Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')
    subplot(4,3,10)
    pcolor(rtmp2,ztmp2,10*log10(zhh_ka_tot+zhhka_smt)-10*log10(zhh_w_tot+zhhw_smt))
    shading flat
    colorbar
    title('Ka-W Z_{HH} (dB)')
    xlabel('r(m)')
    ylabel('z(m)')

    if(conf_flag)
        title_size=24;
        color_bar_size=18;
        tick_size=18;
        label_size=24;
    else
        title_size=14;
        color_bar_size=14;
        tick_size=12;
        label_size=14;
    end

    cmap2=feval('boonlib','carbmap',nanmax(ceil(radii/5)));
    cmap2=cmap2(2:max(size(cmap2)),:);
    urlims=ceil(nanmax(nanmax(abs([urs;urx;urw])))/5)*5;
    urlims=[-urlims urlims];
    vrlims=ceil(nanmax(nanmax(abs([vrs;vrx;vrw])))/5)*5;
    vrlims=[-vrlims vrlims];
    if(make_figs_analysis)
    % cmap2=parula(nanmax(round(radii/10)*10));

    figure(16)
    feval('boonlib','bsizewin',16,[1000 800])
    subplot(3,3,1)
    pcolor(rtmp2,ztmp2,urs)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('S-band u_{r} (m s^{-1})','FontSize',title_size)
    title('u_{r} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    colormap(cmap2)
    caxis(urlims)
    if(~conf_flag)
    text(label_x,label_y,'(a)','FontSize',14,'BackgroundColor',[1 1 1])
    else
        text(-500,350,'S','FontSize',36)
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,2)
    pcolor(rtmp2,ztmp2,vrs)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('S-band v_{r} (m s^{-1})','FontSize',title_size)
    title('v_{r} (m s^{-1})','FontSize',title_size)

    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    axis(lims)
    caxis(vrlims)
    colormap(cmap2)
    if(~conf_flag)
    text(label_x,label_y,'(b)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,3)
    pcolor(rtmp2,ztmp2,uds_rad)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('S-band Dominant Scatt. D (mm)','FontSize',title_size)
    title('Dominant Scatt. R (mm)','FontSize',title_size)

    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    axis(lims)
    colormap(cmap2)
    caxis([0 nanmax(round(radii/10)*10)])
    if(~conf_flag)
    text(label_x,label_y,'(c)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,4)
    pcolor(rtmp2,ztmp2,urx)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('X-band u_{r} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    colormap(cmap2)
    if(~conf_flag)
    text(label_x,label_y,'(d)','FontSize',14,'BackgroundColor',[1 1 1])
    else
        text(-500,350,'X','FontSize',36)
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,5)
    pcolor(rtmp2,ztmp2,vrx)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('X-band v_{r} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    axis(lims)
    caxis(vrlims)
    colormap(cmap2)
    if(~conf_flag)
    text(label_x,label_y,'(e)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,6)
    pcolor(rtmp2,ztmp2,udx_rad)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('X-band Dominant Scatt. D (mm)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    if(~conf_flag)
    text(label_x,label_y,'(f)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)
    axis(lims)
    colormap(cmap2)
    caxis([0 nanmax(round(radii/10)*10)])
    subplot(3,3,7)
    pcolor(rtmp2,ztmp2,urw)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('W-band u_{r} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    colormap(cmap2)
    if(~conf_flag)
    text(label_x,label_y,'(g)','FontSize',14,'BackgroundColor',[1 1 1])
    else
        text(-500,350,'W','FontSize',36)
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,8)
    pcolor(rtmp2,ztmp2,vrw)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('W-band v_{r} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    axis(lims)
    caxis(vrlims)
    colormap(cmap2)
    if(~conf_flag)
    text(label_x,label_y,'(h)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)
    subplot(3,3,9)
    pcolor(rtmp2,ztmp2,udw_rad)
    shading flat
    colorbar('FontSize',color_bar_size)
    % title('W-band Dominant Scatt. D (mm)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    % ylabel('z(m)')
    axis(lims)
    colormap(cmap2)
    caxis([0 nanmax(round(radii/10)*10)])
    if(~conf_flag)
    text(label_x,label_y,'(i)','FontSize',14,'BackgroundColor',[1 1 1])
    end
    set(gca,'FontSize',tick_size)

    save_fig=true;
    if(save_fig)
        cd(fig_dir)
            set(gcf,'PaperPositionMode','auto')
            print(['velr_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        cd(base_dir)
    end

    end

    zhh_lims=[0 60]; ddu_lims=[0 35]; ddv_lims=[-25 25];
    % figure(17)
    % feval('boonlib','bsizewin',17,[1000 800])
    % subplot(3,3,1)
    % pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_w_tot+zhhw_smt))
    % shading flat
    % hcb(1)=colorbar('FontSize',12);
    % title('S-W Z_{HH} Diff (dB)','FontSize',14)
    % % xlabel('r(m)')
    % ylabel('z(m)','FontSize',14)
    % axis(lims)
    % colormap(cmap2)
    % caxis(zhh_lims)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(a)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,2)
    % pcolor(rtmp2,ztmp2,uds-udw)
    % shading flat
    % colorbar
    % title('S-W DDU (m s^{-1})','FontSize',14)
    % % xlabel('r(m)')
    % % ylabel('z(m)')
    % axis(lims)
    % caxis(ddu_lims)
    % colormap(cmap2)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(b)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,3)
    % pcolor(rtmp2,ztmp2,vds-vdw)
    % shading flat
    % colorbar
    % title('S-W DDV (m s^{-1})','FontSize',14)
    % % xlabel('r(m)')
    % % ylabel('z(m)')
    % axis(lims)
    % colormap(cmap2)
    % caxis(ddv_lims)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(c)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,4)
    % pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_x_tot+zhhx_smt))
    % shading flat
    % colorbar
    % title('S-X Z_{HH} Diff (dB)','FontSize',14)
    % % xlabel('r(m)')
    % ylabel('z(m)','FontSize',14)
    % axis(lims)
    % caxis(zhh_lims)
    % colormap(cmap2)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(d)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,5)
    % pcolor(rtmp2,ztmp2,uds-udx)
    % shading flat
    % colorbar
    % title('S-X DDU (m s^{-1})','FontSize',14)
    % % xlabel('r(m)')
    % % ylabel('z(m)')
    % axis(lims)
    % caxis(ddu_lims)
    % colormap(cmap2)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(e)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,6)
    % pcolor(rtmp2,ztmp2,vds-vdx)
    % shading flat
    % colorbar
    % title('S-X DDV (m s^{-1})','FontSize',14)
    % % xlabel('r(m)')
    % % ylabel('z(m)')
    % text(label_x,label_y,'(f)','FontSize',14,'BackgroundColor',[1 1 1])
    % set(gca,'FontSize',14)
    % axis(lims)
    % colormap(cmap2)
    % caxis(ddv_lims)
    % subplot(3,3,7)
    % pcolor(rtmp2,ztmp2,10*log10(zhh_x_tot+zhhx_smt)-10*log10(zhh_w_tot+zhhw_smt))
    % shading flat
    % colorbar
    % title('X-W Z_{HH} Diff (dB)','FontSize',14)
    % xlabel('r(m)','FontSize',14)
    % ylabel('z(m)','FontSize',14)
    % axis(lims)
    % caxis(zhh_lims)
    % set(gca,'FontSize',14)
    % colormap(cmap2)
    % text(label_x,label_y,'(g)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,8)
    % pcolor(rtmp2,ztmp2,udx-udw)
    % shading flat
    % colorbar
    % title('X-W DDU (m s^{-1})','FontSize',14)
    % xlabel('r(m)','FontSize',14)
    % % ylabel('z(m)')
    % axis(lims)
    % caxis(ddu_lims)
    % colormap(cmap2)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(h)','FontSize',14,'BackgroundColor',[1 1 1])
    % subplot(3,3,9)
    % pcolor(rtmp2,ztmp2,vdx-vdw)
    % shading flat
    % colorbar
    % title('X-W DDV (m s^{-1})','FontSize',14)
    % xlabel('r(m)','FontSize',14)
    % % ylabel('z(m)')
    % axis(lims)
    % colormap(cmap2)
    % caxis(ddv_lims)
    % set(gca,'FontSize',14)
    % text(label_x,label_y,'(i)','FontSize',14,'BackgroundColor',[1 1 1])

    if(make_figs)
    figure(17)
    feval('boonlib','bsizewin',17,[1000 800])
    subplot(2,2,1)
    pcolor(rtmp2,ztmp2,uds-udw)
    shading flat
    colorbar('FontSize',16)
    title('S-W DDU (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    caxis(ddu_lims)
    colormap(cmap2)
    set(gca,'FontSize',16)
    subplot(2,2,2)
    pcolor(rtmp2,ztmp2,urs)
    shading flat
    colorbar('FontSize',16)
    title('U_{r} (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    colormap(cmap2)
    caxis(ddu_lims)
    set(gca,'FontSize',16)

    subplot(2,2,3)
    pcolor(rtmp2,ztmp2,udx-udw)
    shading flat
    colorbar('FontSize',16)
    title('X-W DDU (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    caxis(ddu_lims)
    colormap(cmap2)
    set(gca,'FontSize',16)
    subplot(2,2,4)
    pcolor(rtmp2,ztmp2,urx)
    shading flat
    colorbar('FontSize',16)
    title('U_{r} (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    colormap(cmap2)
    caxis(ddu_lims)
    set(gca,'FontSize',16)

    figure(18)
    feval('boonlib','bsizewin',18,[1000 800])
    subplot(2,2,1)
    pcolor(rtmp2,ztmp2,vds-vdw)
    shading flat
    colorbar('FontSize',16)
    title('S-W DDV (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    caxis(ddv_lims)
    colormap(cmap2)
    set(gca,'FontSize',16)
    subplot(2,2,2)
    pcolor(rtmp2,ztmp2,vrs)
    shading flat
    colorbar('FontSize',16)
    title('V_{r} (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    colormap(cmap2)
    caxis(ddv_lims)
    set(gca,'FontSize',16)

    subplot(2,2,3)
    pcolor(rtmp2,ztmp2,vdx-vdw)
    shading flat
    colorbar('FontSize',16)
    title('X-W DDV (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    caxis(ddv_lims)
    colormap(cmap2)
    set(gca,'FontSize',16)
    subplot(2,2,4)
    pcolor(rtmp2,ztmp2,vrx)
    shading flat
    colorbar('FontSize',16)
    title('V_{r} (m s^{-1})','FontSize',24)
    xlabel('r(m)','FontSize',18)
    ylabel('z(m)','FontSize',18)
    axis(lims)
    colormap(cmap2)
    caxis(ddv_lims)
    set(gca,'FontSize',16)

    save_fig_paper=false;
    if(save_fig_paper)
        cd(fig_dir)
            set(gcf,'PaperPositionMode','auto')
            print(['ddv_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        cd(base_dir)
    else
        figure(17)
        cd(fig_dir)
            set(gcf,'PaperPositionMode','auto')
            print(['ddu_conf_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        cd(base_dir)
        figure(18)
        cd(fig_dir)
            set(gcf,'PaperPositionMode','auto')
            print(['ddv_conf_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        cd(base_dir)
    end

    end

    rmse_su=sqrt(nanmean(nanmean(((uds-udw)-urs).^2)));
    rmse_sv=sqrt(nanmean(nanmean(((vds-vdw)-vrs).^2)));
    rmse_xu=sqrt(nanmean(nanmean(((udx-udw)-urx).^2)));
    rmse_xv=sqrt(nanmean(nanmean(((vdx-vdw)-vrx).^2)));
    rmse_cu=sqrt(nanmean(nanmean(((udc-udw)-urc).^2)));
    rmse_cv=sqrt(nanmean(nanmean(((vdc-vdw)-vrc).^2)));
    rmse_kau=sqrt(nanmean(nanmean(((udka-udw)-urka).^2)));
    rmse_kav=sqrt(nanmean(nanmean(((vdka-vdw)-vrka).^2)));

    udlims=[-40 40]; vdlims=[0 70]; urlims=[0 40]; vrlims=[-30 30]; zhdiff_lims=[0 60]; dom_lims=[0 100];
    title_size=14;
    color_bar_size=14;
    tick_size=12;
    label_size=14;
    label_x=50; label_y=550;

    if(make_figs)
    figure(200)
    feval('boonlib','bsizewin',200,[1000 1050])
    subplot(4,2,1)
    pcolor(rtmp2,ztmp2,uds)
    % title('u_{d} (m s^{-1})','FontSize',title_size)
    title('u_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(udlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size)
    set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(a)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,2)
    pcolor(rtmp2,ztmp2,vds)
    title('v_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vdlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(b)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,3)
    pcolor(rtmp2,ztmp2,urs)
    title('u_{dr} - U (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(c)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,4)
    pcolor(rtmp2,ztmp2,vrs)
    title('v_{dr} - V (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vrlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(d)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,5)
    pcolor(rtmp2,ztmp2,uds-udw)
    title('S-W DDU (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(e)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,6)
    pcolor(rtmp2,ztmp2,vds-vdw)
    title('S-W DDV (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vrlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(f)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,7)
    pcolor(rtmp2,ztmp2,10*log10(zhh_s_tot+zhhs_smt)-10*log10(zhh_w_tot+zhhw_smt))
    title('S-W Z_{e} (dB)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(zhdiff_lims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(g)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,8)
    pcolor(rtmp2,ztmp2,uds_rad)
    title('Dom. Scatt. Rad. (mm)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(dom_lims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(h)','FontSize',14,'BackgroundColor',[1 1 1])

    figure(201)
    feval('boonlib','bsizewin',201,[1000 1050])
    subplot(4,2,1)
    pcolor(rtmp2,ztmp2,udx)
    title('u_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(udlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(a)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,2)
    pcolor(rtmp2,ztmp2,vdx)
    title('v_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vdlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(b)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,3)
    pcolor(rtmp2,ztmp2,urx)
    title('u_{dr} - U (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(c)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,4)
    pcolor(rtmp2,ztmp2,vrx)
    title('v_{dr} - V (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vrlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(d)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,5)
    pcolor(rtmp2,ztmp2,udx-udw)
    title('X-W DDU (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(e)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,6)
    pcolor(rtmp2,ztmp2,vdx-vdw)
    title('X-W DDV (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vrlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(f)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,7)
    pcolor(rtmp2,ztmp2,10*log10(zhh_x_tot+zhhx_smt)-10*log10(zhh_w_tot+zhhw_smt))
    title('X-W Z_{e} (dB)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(zhdiff_lims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(g)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(4,2,8)
    pcolor(rtmp2,ztmp2,udx_rad)
    title('Dom. Scatt. Rad. (mm)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(dom_lims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(h)','FontSize',14,'BackgroundColor',[1 1 1])

    figure(202)
    feval('boonlib','bsizewin',202,[1000 850])
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,udw)
    title('u_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(udlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(a)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(3,2,2)
    pcolor(rtmp2,ztmp2,vdw)
    title('v_{dr} (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vdlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(b)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,urw)
    title('u_{dr} - U (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(urlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(c)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(3,2,4)
    pcolor(rtmp2,ztmp2,vrw)
    title('v_{dr} - V (m s^{-1})','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(vrlims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(d)','FontSize',14,'BackgroundColor',[1 1 1])
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,udw_rad)
    title('Dom. Scatt. Rad. (mm)','FontSize',title_size)
    xlabel('r(m)','FontSize',label_size)
    ylabel('z(m)','FontSize',label_size)
    axis(lims)
    caxis(dom_lims)
    shading flat; colormap(cmap2)
    colorbar('FontSize',color_bar_size); set(gca,'FontSize',tick_size)
    text(label_x,label_y,'(e)','FontSize',14,'BackgroundColor',[1 1 1])

    cd(fig_dir)
        figure(200)
        set(gcf,'PaperPositionMode','auto')
        print(['all_S_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        figure(201)
        set(gcf,'PaperPositionMode','auto')
        print(['all_X_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
        figure(202)
        set(gcf,'PaperPositionMode','auto')
        print(['all_W_' em_title '_' dry_title '_' num2str(deb_scale_large) '_' num2str(rain_scale) '_' num2str(orient_case) '.png'],'-dpng');
    cd(base_dir)

    end

    ind2=1:9;
    urs_mean=nanmean(nanmean(urs(:,ind2)));
    urc_mean=nanmean(nanmean(urc(:,ind2)));
    urx_mean=nanmean(nanmean(urx(:,ind2)));
    urka_mean=nanmean(nanmean(urka(:,ind2)));
    urw_mean=nanmean(nanmean(urw(:,ind2)));
    % vrs_mean=nanmean(nanmean(vrs(:,ind2)));
    % vrc_mean=nanmean(nanmean(vrc(:,ind2)));
    % vrx_mean=nanmean(nanmean(vrx(:,ind2)));
    % vrka_mean=nanmean(nanmean(vrka(:,ind2)));
    % vrw_mean=nanmean(nanmean(vrw(:,ind2)));

    % rmax=nanmax(nanmax(rtmp2))+10; zmax=nanmax(nanmax(ztmp2))+10;
    rmax=700; zmax=300;

    rind=1:sum(rtmp2(:,1)<rmax);
    zind=1:sum(ztmp2(1,:)<zmax);
    xwdu=udx(rind,zind)-udw(rind,zind); xwdu=reshape(xwdu,size(xwdu,1)*size(xwdu,2),1);
    swdu=uds(rind,zind)-udw(rind,zind); swdu=reshape(swdu,size(swdu,1)*size(swdu,2),1);
    xwdv=vdx(rind,zind)-vdw(rind,zind); xwdv=reshape(xwdv,size(xwdv,1)*size(xwdv,2),1);
    swdv=vds(rind,zind)-vdw(rind,zind); swdv=reshape(swdv,size(swdv,1)*size(swdv,2),1);

    zhd_sw=10*log10(zhh_s_tot(rind,zind)+zhhs_smt(rind,zind))-10*log10(zhh_w_tot(rind,zind)+zhhw_smt(rind,zind));
    % zhd_sw=repmat(zhd_sw,size(zhd_sw,1)*size(zhd_sw,2),1);
    zhd_xw=10*log10(zhh_x_tot(rind,zind)+zhhx_smt(rind,zind))-10*log10(zhh_w_tot(rind,zind)+zhhw_smt(rind,zind));
    % zhd_xw=repmat(zhd_xw,size(zhd_xw,1)*size(zhd_xw,2),1);

    xw_thres=20; sw_thres=-1000;
    urwt=reshape(urx(rind,zind),size(urx(rind,zind),1)*size(urx(rind,zind),2),1);
    yn=isfinite(xwdu)&isfinite(urwt)&reshape(zhd_xw,size(zhd_xw,1)*size(zhd_xw,2),1)>xw_thres;
    xwucoeff=corrcoef(xwdu(yn),urwt(yn));

    urst=reshape(urs(rind,zind),size(urs(rind,zind),1)*size(urs(rind,zind),2),1);
    yn=isfinite(swdu)&isfinite(urst)&reshape(zhd_sw,size(zhd_sw,1)*size(zhd_sw,2),1)>sw_thres;
    swucoeff=corrcoef(swdu(yn),urst(yn));

    vrwt=reshape(vrx(rind,zind),size(vrx(rind,zind),1)*size(vrx(rind,zind),2),1);
    yn=isfinite(xwdv)&isfinite(vrwt)&reshape(zhd_xw,size(zhd_xw,1)*size(zhd_xw,2),1)>xw_thres;
    xwvcoeff=corrcoef(xwdv(yn),vrwt(yn));

    vrst=reshape(vrs(rind,zind),size(vrs(rind,zind),1)*size(vrs(rind,zind),2),1);
    yn=isfinite(swdv)&isfinite(vrst)&reshape(zhd_sw,size(zhd_sw,1)*size(zhd_sw,2),1)>sw_thres;
    swvcoeff=corrcoef(swdv(yn),vrst(yn));

    % zhd_xw=10*log10(zhh_x_tot(rind,zind)+zhhx_smt(rind,zind))-10*log10(zhh_w_tot(rind,zind)+zhhw_smt(rind,zind));
    zhd_xw=repmat(zhd_xw,size(zhd_xw,1)*size(zhd_xw,2),1);
    udx_radt=repmat(udx_rad(rind,zind),size(udx_rad(rind,zind),1)*size(udx_rad(rind,zind),2),1);
    yn=isfinite(zhd_xw)&isfinite(udx_radt);
    zhd_xw_coeff=corrcoef(zhd_xw(yn),udx_radt(yn));

    % zhd_sw=10*log10(zhh_s_tot(rind,zind)+zhhs_smt(rind,zind))-10*log10(zhh_w_tot(rind,zind)+zhhw_smt(rind,zind));
    zhd_sw=repmat(zhd_sw,size(zhd_sw,1)*size(zhd_sw,2),1);
    uds_radt=repmat(uds_rad(rind,zind),size(uds_rad(rind,zind),1)*size(uds_rad(rind,zind),2),1);
    yn=isfinite(zhd_sw)&isfinite(uds_radt);
    zhd_sw_coeff=corrcoef(zhd_sw(yn),uds_radt(yn));

    clear uds_radt;
    zhd_sx=10*log10(zhh_s_tot(rind,zind)+zhhs_smt(rind,zind))-10*log10(zhh_x_tot(rind,zind)+zhhx_smt(rind,zind));
    zhd_sx=repmat(zhd_sx,size(zhd_sx,1)*size(zhd_sx,2),1);
    uds_radt=repmat(uds_rad(rind,zind),size(uds_rad(rind,zind),1)*size(uds_rad(rind,zind),2),1);
    yn=isfinite(zhd_sx)&isfinite(uds_radt);
    zhd_sx_coeff=corrcoef(zhd_sx(yn),uds_radt(yn));

    rind=1:20;
    zind=1:15;
    urs_max=nanmax(nanmax(urs(rind,zind)));
    urc_max=nanmax(nanmax(urc(rind,zind)));
    urx_max=nanmax(nanmax(urx(rind,zind)));
    urka_max=nanmax(nanmax(urka(rind,zind)));
    urw_max=nanmax(nanmax(urw(rind,zind)));
    vrs_min=nanmin(nanmin(vrs(rind,zind)));
    vrc_min=nanmin(nanmin(vrc(rind,zind)));
    vrx_min=nanmin(nanmin(vrx(rind,zind)));
    vrka_min=nanmin(nanmin(vrka(rind,zind)));
    vrw_min=nanmin(nanmin(vrw(rind,zind)));
    urs_tmp=urs(rind,zind); urs_tmp=urs_tmp(:);
    % urc_tmp=urc(rind,zind); urc_tmp=urc_tmp(:);
    urx_tmp=urx(rind,zind); urx_tmp=urx_tmp(:);
    % urka_tmp=urka(rind,zind); urka_tmp=urka_tmp(:);
    urw_tmp=urw(rind,zind); urw_tmp=urw_tmp(:);
    ddu_sw=uds(rind,zind)-udw(rind,zind); ddu_sw_tmp=ddu_sw(:);
    ddu_xw=udx(rind,zind)-udw(rind,zind); ddu_xw_tmp=ddu_xw(:);
    urs_90=prctile(urs_tmp,90);
    urs_50=prctile(urs_tmp,50);
    % urc_90=prctile(urc_tmp,90);
    % urc_50=prctile(urc_tmp,50);
    urx_90=prctile(urx_tmp,90);
    urx_50=prctile(urx_tmp,50);
    % urka_90=prctile(urka_tmp,90);
    % urka_50=prctile(urka_tmp,50);
    urw_90=prctile(urw_tmp,90);
    urw_50=prctile(urw_tmp,50);
    ddu_sw_90=prctile(ddu_sw_tmp,90);
    ddu_sw_50=prctile(ddu_sw_tmp,50);
    ddu_xw_90=prctile(ddu_xw_tmp,90);
    ddu_xw_50=prctile(ddu_xw_tmp,50);

    %New Theory test
    % lb=0.107; %m
    % ls=0.055;
    u_act=uds-urs;
    v_act=vds-vrs;
    load Ze_wood_rcs.mat;
    if(dry_flag)
        Ze_wood=10.^(Ze_dry./10);
    else
        Ze_wood=10.^(Ze_wet./10);
    end
    rain_ind=[15 18 20 21 25];
    deb_scale_sm=[2.9497    0.4724    0.1584    0.0747    0.0052]*1e5*(225^2/150^2)^3*rain_scale;

    if(lb==0.107)
        udb_tmp=uds;
        vdb_tmp=vds;
        zhb_sm=zhhs_smt;
        zhb_tot=zhh_s_tot;
        zhb_sm_all=zhhs_sm;
        zhb_tot_all=zhh_s_store;
        zhb_sm_ind=load('rain_S.mat');
        zhb_tot_ind=Ze_wood(:,5);
    elseif(lb==0.055)
        udb_tmp=udc;
        vdb_tmp=vdc;
        zhb_sm=zhhc_smt;
        zhb_tot=zhh_c_tot;
        zhb_sm_all=zhhc_sm;
        zhb_tot_all=zhh_c_store;
        zhb_sm_ind=load('rain_C.mat');
        zhb_tot_ind=Ze_wood(:,4);
    elseif(lb==0.031)
        udb_tmp=udx;
        vdb_tmp=vdx;
        zhb_sm=zhhx_smt;
        zhb_tot=zhh_x_tot;
        zhb_sm_all=zhhx_sm;
        zhb_tot_all=zhh_x_store;
        zhb_sm_ind=load('rain_X.mat');
        zhb_tot_ind=Ze_wood(:,3);
    end
    if(ls==0.055)
        uds_tmp=udc;
        vds_tmp=vdc;
        zhs_sm=zhhc_smt;
        zhs_tot=zhh_c_tot;
        zhs_sm_all=zhhc_sm;
        zhs_tot_all=zhh_c_store;
        zhs_sm_ind=load('rain_C.mat');
        zhs_tot_ind=Ze_wood(:,4);
    elseif(ls==0.031)
        uds_tmp=udx;
        vds_tmp=vdx;
        zhs_sm=zhhx_smt;
        zhs_tot=zhh_x_tot;
        zhs_sm_all=zhhx_sm;
        zhs_tot_all=zhh_x_store;
        zhs_sm_ind=load('rain_X.mat');
        zhs_tot_ind=Ze_wood(:,3);
    elseif(ls==0.008)
        uds_tmp=udka;
        vds_tmp=vdka;
        zhs_sm=zhhka_smt;
        zhs_tot=zhh_ka_tot;
        zhs_sm_all=zhhka_sm;
        zhs_tot_all=zhh_ka_store;
        zhs_sm_ind=load('rain_Ka.mat');
        zhs_tot_ind=Ze_wood(:,2);
    elseif(ls==0.003)
        uds_tmp=udw;
        vds_tmp=vdw;
        zhs_sm=zhhw_smt;
        zhs_tot=zhh_w_tot;
        zhs_sm_all=zhhw_sm;
        zhs_tot_all=zhh_w_store;
        zhs_sm_ind=load('rain_W.mat');
        zhs_tot_ind=Ze_wood(:,1);
    end
    c1=ls^4/lb^4; km=0.93;
    dwr=(zhb_tot+zhb_sm)./(zhs_tot+zhs_sm);
    dwr_noise=(zhb_tot+zhb_sm+b_noise*randn(size(zhb_tot)))./(zhs_tot+zhs_sm+s_noise*randn(size(zhb_tot)));

    zlb=zhb_tot+zhb_sm;
    % dwr_noise(dwr_noise>1/c1*0.95)=1/c1*0.95;
    zlb2_divide_zlb1=(dwr_noise-1)./(1-dwr_noise*c1);

    zlb1=zlb./(1+zlb2_divide_zlb1);
    zlb2=zlb-zlb1;

    zls1=zlb1;
    zls2=c1*zlb2;

    nls1=pi^5/ls^4*km*zls1;
    nls2=pi^5/ls^4*km*zls2;
    nlb1=pi^5/lb^4*km*zlb1;
    nlb2=pi^5/lb^4*km*zlb2;
    nls=nls1+nls2;
    nlb=nlb1+nlb2;

    % vls=(vr*nls1+vm*nls2)./(nls1+nls2);
    % vlb=(vr*nlb1+vm*nlb2)./(nlb1+nlb2);
    %Subscript 1: Rayleigh, subscript 2: Mie
    uls=uds_tmp; ulb=udb_tmp;
    vls=vds_tmp; vlb=vdb_tmp;

    vm_app=(vls.*nls-vlb.*nlb*1/c1)./(nls2-nlb2*1/c1); %Retrieved velocity for Mie/optical scatterers
    vr_app=(vlb.*nlb-vm_app.*nlb2)./nlb1; %Retrived velocity for Rayleigh scatterers
    um_app=(uls.*nls-ulb.*nlb*1/c1)./(nls2-nlb2*1/c1); %Retrieved velocity for Mie/optical scatterers
    ur_app=(ulb.*nlb-um_app.*nlb2)./nlb1; %Retrived velocity for Rayleigh scatterers

    % ur_app=(uls.*nls./nls1-ulb.*nlb./nls1)./(1-nlb1./nls1);
    % vr_app=(vls.*nls./nls1-vlb.*nlb./nls1)./(1-nlb1./nls1);
    % um_app=ulb.*nlb./nlb2-nlb1./nlb2.*ur_app;
    % vm_app=vlb.*nlb./nlb2-nlb1./nlb2.*vr_app;

    yn=dwr_noise<dwr_thres|dwr_noise>dwr_max;
    um_app(yn)=NaN;
    yn2=nlb1==0|dwr_noise<dwr_thres|dwr_noise>dwr_max;
    ur_app(yn2)=NaN;
    vm_app(yn)=NaN;
    vr_app(yn2)=NaN;

    figure(1003)
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,um_app)
    shading flat
    colorbar
    caxis([-35 35])
    title('Retrieved Mie U_{d}')
    subplot(3,2,2)
    pcolor(rtmp2,ztmp2,udb_tmp)
    shading flat
    colorbar
    title(['U_{d} at lambda=' num2str(lb*1000) ' mm'])
    caxis([-35 35])
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,ur_app)
    shading flat
    colorbar
    caxis([-35 35])
    title('Retrieved Rayleigh U_{d}')
    subplot(3,2,4)
    pcolor(rtmp2,ztmp2,uds_tmp)
    shading flat
    colorbar
    caxis([-35 35])
    title(['U_{d} at lambda=' num2str(ls*1000) ' mm'])
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,u_act)
    shading flat
    colorbar
    caxis([-35 35])
    title('Actual U')

    figure(1004)
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,vm_app)
    shading flat
    colorbar
    caxis([0 70])
    title(['Ret. Mie V_{d} at lambda=' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,2)
    pcolor(rtmp2,ztmp2,vdb_tmp)
    shading flat
    colorbar
    caxis([0 70])
    title(['V_{d} at lambda=' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,vr_app)
    shading flat
    colorbar
    caxis([0 70])
    title(['Ret. Rayleigh V_{d} at lambda=' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,4)
    pcolor(rtmp2,ztmp2,vds_tmp)
    shading flat
    colorbar
    caxis([0 70])
    title(['V_{d} at lambda=' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,v_act)
    shading flat
    colorbar
    caxis([0 70])
    title(['V'])
    xlabel('R(m)')
    ylabel('Z(m)')

    figure(1005)
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,10*log10(zhb_sm))
    shading flat
    caxis([0 70])
    colorbar
    title(['Small Scatterer Ze at ' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,2)
    zlb1(zlb1<0)=NaN;
    zlb1_log=10*log10(zlb1);
    zlb1_log(~isfinite(zlb1_log))=NaN;
    pcolor(rtmp2,ztmp2,zlb1_log)
    shading flat
    title(['Ret. Rayleigh Scatterer Ze at' num2str(lb*1000) ' mm'])
    caxis([0 70])
    colorbar
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,10*log10(zhb_tot))
    shading flat
    caxis([0 70])
    colorbar
    title(['Large Scatterer Ze at ' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,4)
    zlb2(zlb2<0)=NaN;
    zlb2_log=10*log10(zlb2);
    zlb2_log(~isfinite(zlb2_log))=NaN;
    pcolor(rtmp2,ztmp2,zlb2_log)
    shading flat
    caxis([0 70])
    colorbar
    title(['Ret. Mie. Ze at ' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,10*log10(zhb_tot+zhb_sm))
    shading flat
    caxis([0 70])
    colorbar
    title(['Z_{e} at ' num2str(lb*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')

    figure(1006)
    subplot(3,2,1)
    pcolor(rtmp2,ztmp2,10*log10(zhs_sm))
    shading flat
    caxis([0 70])
    colorbar
    title(['Small Scatter Ze at ' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,2)
    zls1(zls1<0)=NaN;
    zls1_log=10*log10(zls1);
    zls1_log(~isfinite(zls1_log))=NaN;
    pcolor(rtmp2,ztmp2,zls1_log)
    shading flat
    title(['Ret. Small Scatterer Ze at ' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    caxis([0 70])
    colorbar
    subplot(3,2,3)
    pcolor(rtmp2,ztmp2,10*log10(zhs_tot))
    shading flat
    caxis([0 70])
    colorbar
    title(['Large Scatterer Ze at ' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,4)
    zls2(zls2<0)=NaN;
    zls2_log=10*log10(zls2);
    zls2_log(~isfinite(zls2_log))=NaN;
    pcolor(rtmp2,ztmp2,zls2_log)
    shading flat
    caxis([0 70])
    colorbar
    title(['Ret. Mie Scatter Ze at ' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')
    subplot(3,2,5)
    pcolor(rtmp2,ztmp2,10*log10(zhs_tot+zhs_sm))
    shading flat
    caxis([0 70])
    colorbar
    title(['Z_{e} at ' num2str(ls*1000) ' mm'])
    xlabel('R(m)')
    ylabel('Z(m)')

    %Velocity estimator
    rayu.rmse_lb=sqrt(nanmean(nanmean((nanmean(udt_avg_sm(rind,zind,:),3)-ur_app(rind,zind)).^2)));
    mieu.rmse_lb=sqrt(nanmean(nanmean((nanmean(udt_avg(rind,zind),3)-um_app(rind,zind)).^2)));
    rayu.rmse_dop=sqrt(nanmean(nanmean((nanmean(udt_avg_sm(rind,zind),3)-uds_tmp(rind,zind)).^2)));
    mieu.rmse_dop=sqrt(nanmean(nanmean((nanmean(udt_avg(rind,zind,:),3)-udb_tmp(rind,zind)).^2)));
    rayv.rmse_lb=sqrt(nanmean(nanmean((nanmean(vdt_avg_sm(rind,zind,:),3)-vr_app(rind,zind)).^2)));
    miev.rmse_lb=sqrt(nanmean(nanmean((nanmean(vdt_avg(rind,zind),3)-vm_app(rind,zind)).^2)));
    rayv.rmse_dop=sqrt(nanmean(nanmean((nanmean(vdt_avg_sm(rind,zind,:),3)-vds_tmp(rind,zind)).^2)));
    miev.rmse_dop=sqrt(nanmean(nanmean((nanmean(vdt_avg(rind,zind),3)-vdb_tmp(rind,zind)).^2)));

    zhb_sm_orig=zhb_sm;
    zhs_sm_orig=zhs_sm;
    zhs_tot_orig=zhs_tot;
    zhb_tot_orig=zhb_tot;

    zhb_sm(zhb_sm<=0)=NaN;
    zhs_sm(zhs_sm<=0)=NaN;
    zhs_tot(zhs_tot<=0)=NaN;
    zhb_tot(zhb_tot<=0)=NaN;
    %Reflectivity estimator
    %Retrieval error in Rayleigh scatter at long wavelength
    rayz.rmse_lb=sqrt(nanmean(nanmean((10*log10(zhb_sm(rind,zind))-zlb1_log(rind,zind)).^2))); 
    rayz.bias_lb=nanmean(nanmean(10*log10(zhb_sm(rind,zind))-zlb1_log(rind,zind)));
    %Retrieval error in Rayleigh scatter at small wavelength
    rayz.rmse_ls=sqrt(nanmean(nanmean((10*log10(zhs_sm(rind,zind))-zls1_log(rind,zind)).^2))); 
    rayz.bias_ls=nanmean(nanmean(10*log10(zhs_sm(rind,zind))-zls1_log(rind,zind)));
    %Retrieval error in Mie scatter at long wavelength
    miez.rmse_lb=sqrt(nanmean(nanmean((10*log10(zhb_tot(rind,zind))-zlb2_log(rind,zind)).^2))); 
    miez.bias_lb=nanmean(nanmean(10*log10(zhb_tot(rind,zind))-zlb2_log(rind,zind)));
    %Retrieval error in Mie scatter at small wavelength
    miez.rmse_ls=sqrt(nanmean(nanmean((10*log10(zhs_tot(rind,zind))-zls2_log(rind,zind)).^2)));
    miez.bias_ls=nanmean(nanmean(10*log10(zhs_tot(rind,zind))-zls2_log(rind,zind)));
    %Retrieval error assuming short wavelength Z is Rayleigh
    rayz.rmse_dop=sqrt(nanmean(nanmean((10*log10(zhs_sm(rind,zind))-10*log10(zhs_tot(rind,zind)+zhs_sm(rind,zind))).^2)));
    rayz.bias_dop=nanmean(nanmean(10*log10(zhs_tot(rind,zind))-10*log10(zhs_tot(rind,zind)+zhs_sm(rind,zind))));
    %Retrieval error assuming long wavelength Z is Mie
    miez.rmse_dop=sqrt(nanmean(nanmean((10*log10(zhb_tot(rind,zind))-10*log10(zhb_tot(rind,zind)+zhb_sm(rind,zind))).^2)));
    miez.bias_dop=nanmean(nanmean(10*log10(zhb_tot(rind,zind))-10*log10(zhb_tot(rind,zind)+zhb_sm(rind,zind))));

    figure(1010)
    subplot(2,2,1)
    pcolor(rtmp2(rind,zind),ztmp2(rind,zind),10*log10(zhb_tot(rind,zind))-10*log10(zhb_tot(rind,zind)+zhb_sm(rind,zind)))
    shading flat
    colorbar
    subplot(2,2,2)
    pcolor(rtmp2(rind,zind),ztmp2(rind,zind),10*log10(zhb_tot(rind,zind)))
    shading flat
    colorbar; caxis([0 70])
    subplot(2,2,3)
    pcolor(rtmp2(rind,zind),ztmp2(rind,zind),10*log10(zhb_sm(rind,zind)))
    shading flat
    colorbar; caxis([0 70])
    subplot(2,2,4)
    pcolor(rtmp2(rind,zind),ztmp2(rind,zind),10*log10(zhb_tot(rind,zind)+zhb_sm(rind,zind)))
    shading flat
    colorbar; caxis([0 70])

    figure(1011)
    subplot(2,2,1)
    pcolor(rtmp2(rind,zind),ztmp2(rind,zind),10*log10(zhs_sm(rind,zind))-10*log10(zhs_tot(rind,zind)+zhs_sm(rind,zind)))
    shading flat
    colorbar

    dwr_mean=10*log10(nanmean(nanmean((zhb_tot(rind,zind)+zhb_sm(rind,zind))./(zhs_tot(rind,zind)+zhs_sm(rind,zind)))));

    if(spectral_flag)
    %Spectra (p.275, eq. 8.77 DZ); Sekelsky (2001)
    r=10000;
    Pt=10000;
    npts=2*va*10;
    kw=0.93;
    mm6_to_m6=1e-18;

    % keyboard;
    sm_Z_lb=zhb_sm_ind.zh_r(rain_ind);
    sm_Z_ls=zhs_sm_ind.zh_r(rain_ind);
    eval('cd spec/woodsphere3_radar/')
    for idx=1:size(zhb_tot_all,3)
        load(['spec_' num2str(idx) '.mat'])
        if(~exist('area','var'))
            area=pi*(rcoords(2:max(size(rcoords))).^2-rcoords(1:max(size(rcoords))-1).^2);
            area=repmat(area',1,max(size(zcoords))-1,size(uspec,3));
            height=zcoords(2:max(size(zcoords)))-zcoords(1:max(size(zcoords))-1);
            height=repmat(height',max(size(rcoords))-1,1,size(uspec,3));
            vol_factor=1./(area.*height);
        end
        Su_deb_lb(:,:,:,idx)=pi^5/lb^4*kw*zhb_tot_ind(idx)*mm6_to_m6*uspec.*vol_factor;
        Su_deb_ls(:,:,:,idx)=pi^5/ls^4*kw*zhs_tot_ind(idx)*mm6_to_m6*uspec.*vol_factor;
        Sv_deb_lb(:,:,:,idx)=pi^5/lb^4*kw*zhb_tot_ind(idx)*mm6_to_m6*vspec.*vol_factor;
        Sv_deb_ls(:,:,:,idx)=pi^5/ls^4*kw*zhs_tot_ind(idx)*mm6_to_m6*vspec.*vol_factor;
    end
    eval('cd ..')
    scale_fact_sm=100;
    eval('cd smalldeb3_radar')
    for idx=size(zhs_tot_all,3)+1:size(zhs_tot_all,3)+size(zhs_sm_all,3)
        load(['spec_' num2str(idx-size(zhs_tot_all,3)) '.mat'])
        Su_deb_lb(:,:,:,idx)=deb_scale_sm(idx-size(zhs_tot_all,3))*scale_fact_sm*pi^5/lb^4*kw*sm_Z_lb(idx-size(zhs_tot_all,3))*mm6_to_m6*uspec.*vol_factor;
        Su_deb_ls(:,:,:,idx)=deb_scale_sm(idx-size(zhs_tot_all,3))*scale_fact_sm*pi^5/ls^4*kw*sm_Z_ls(idx-size(zhs_tot_all,3))*mm6_to_m6*uspec.*vol_factor;
        Sv_deb_lb(:,:,:,idx)=deb_scale_sm(idx-size(zhs_tot_all,3))*scale_fact_sm*pi^5/lb^4*kw*sm_Z_lb(idx-size(zhs_tot_all,3))*mm6_to_m6*vspec.*vol_factor;
        Sv_deb_ls(:,:,:,idx)=deb_scale_sm(idx-size(zhs_tot_all,3))*scale_fact_sm*pi^5/ls^4*kw*sm_Z_ls(idx-size(zhs_tot_all,3))*mm6_to_m6*vspec.*vol_factor;
    end

    ridx=4; zidx=3;
    figure(1500)
    plot(vctrs,squeeze(10*log10(nansum(Su_deb_ls(ridx,zidx,:,1:15),4))),'-k',...
    vctrs,squeeze(10*log10(nansum(Su_deb_lb(ridx,zidx,:,1:15),4))),'-b')

    figure(1501)
    plot(vctrs,squeeze(10*log10(nansum(Sv_deb_ls(ridx,zidx,:,1:15),4))),'-k',...
    vctrs,squeeze(10*log10(nansum(Sv_deb_lb(ridx,zidx,:,1:15),4))),'-b')

    c0=ls^4/lb^4;
    dwr_su=nansum(Su_deb_lb(ridx,zidx,:,1:15),4)./nansum(Su_deb_ls(ridx,zidx,:,1:15),4)*c0;
    figure(1502)
    plot(vctrs,squeeze(10*log10(dwr_su)))

    eval('cd ../..')

    v=linspace(-va,va,npts+1);
    % vr=vctrs;
    vr_rep(1,1,1,1:max(size(vr)))=vr;

    S_deb_lb(:,:,1:size(zhb_tot_all,3))=pi^5/lb^4*kw*zhb_tot_all*mm6_to_m6;
    S_deb_ls(:,:,1:size(zhs_tot_all,3))=pi^5/ls^4*kw*zhs_tot_all*mm6_to_m6;
    S_deb_lb(:,:,size(zhb_tot_all,3)+1:size(zhb_tot_all,3)+size(zhb_sm_all,3))=pi^5/lb^4*kw*zhb_sm_all*mm6_to_m6;
    S_deb_ls(:,:,size(zhs_tot_all,3)+1:size(zhs_tot_all,3)+size(zhs_sm_all,3))=pi^5/ls^4*kw*zhs_sm_all*mm6_to_m6;

    Sn_deb_lb=S_deb_lb./repmat(nansum(S_deb_lb,3),[1 1 size(S_deb_lb,3)]); 
    Sn_deb_ls=S_deb_ls./repmat(nansum(S_deb_ls,3),[1 1 size(S_deb_ls,3)]);

    ud=cat(3,udt_avg,udt_avg_sm);
    vd=cat(3,vdt_avg,vdt_avg_sm);

    vr_rep=repmat(vr_rep,[size(Sn_deb_lb,1) size(Sn_deb_lb,2) size(Sn_deb_lb,3) 1]);
    sigma1_rep=repmat(sigma1,[size(Sn_deb_lb,1) size(Sn_deb_lb,2) 1 size(vr_rep,4)]);
    %U spectra (model as Gaussians)
    %p1=a1*1/(sigma1*2*sqrt(pi))*exp(-(v-mu1).^2/(2*sigma1^2));
    Snu_all_lb=repmat(Sn_deb_lb,[1 1 1 max(size(vr))])./(sigma1_rep*sqrt(2*pi)).*exp(-(vr_rep-repmat(ud,[1 1 1 max(size(vr))])).^2./(2*sigma1_rep.^2))*(vr(2)-vr(1));
    Snu_all_ls=repmat(Sn_deb_ls,[1 1 1 max(size(vr))])./(sigma1_rep*sqrt(2*pi)).*exp(-(vr_rep-repmat(ud,[1 1 1 max(size(vr))])).^2./(2*sigma1_rep.^2))*(vr(2)-vr(1));
    Snu_lb=squeeze(nansum(Snu_all_lb,3));
    Snu_ls=squeeze(nansum(Snu_all_ls,3));

    Snv_all_lb=repmat(Sn_deb_lb,[1 1 1 max(size(vr))])./(sigma1_rep*sqrt(2*pi)).*exp(-(vr_rep-repmat(vd,[1 1 1 max(size(vr))])).^2./(2*sigma1_rep.^2)); %*(vr(2)-vr(1));
    Snv_all_ls=repmat(Sn_deb_ls,[1 1 1 max(size(vr))])./(sigma1_rep*sqrt(2*pi)).*exp(-(vr_rep-repmat(vd,[1 1 1 max(size(vr))])).^2./(2*sigma1_rep.^2)); %*(vr(2)-vr(1));
    Snv_lb=squeeze(nansum(Snv_all_lb,3));
    Snv_ls=squeeze(nansum(Snv_all_ls,3));
    clear Snu_all_lb Snu_all_ls Snv_all_lb Snv_all_ls sigma1_rep;

    %Rayleigh/Mie retrievals
    % con_fact=lb^4/ls^4;
    con_fact=1;
    Snu_lb=Snu_lb*con_fact+noise_floor_db;
    Snv_lb=Snv_lb*con_fact+noise_floor_db;
    Snu_ls=Snu_ls+noise_floor_db;
    Snv_ls=Snv_ls+noise_floor_db;

    dwsu=Snu_lb./Snu_ls;
    c2=(ls^4/lb^4-dwsu)./(dwsu-1);
    Snu_ls_ray=Snu_ls./(1+c2);
    Snu_ls_mie=Snu_ls-Snu_ls_ray;
    Snu_lb_mie=Snu_ls_mie;
    Snu_lb_ray=Snu_ls_ray*ls^4/lb^4;
    dwsv=Snv_lb./Snv_ls;
    c2=(ls^4/lb^4-dwsv)./(dwsv-1);
    Snv_ls_ray=Snv_ls./(1+c2);
    Snv_ls_mie=Snv_ls-Snv_ls_ray;
    Snv_lb_mie=Snv_ls_mie;
    Snv_lb_ray=Snv_ls_ray*ls^4/lb^4;

    %Compute Aliased version
    Snu_lb_al=zeros(size(Snu_lb,1),size(Snu_lb,2),npts);
    Snu_ls_al=Snu_lb_al; Snv_lb_al=Snu_lb_al; Snv_ls_al=Snu_lb_al;
    Snu_lb_al_ray=zeros(size(Snu_lb,1),size(Snu_lb,2),npts);
    Snu_ls_al_ray=Snu_lb_al_ray; Snv_lb_al_ray=Snu_lb_al_ray; Snv_ls_al_ray=Snu_lb_al_ray;
    Snu_lb_al_mie=zeros(size(Snu_lb,1),size(Snu_lb,2),npts);
    Snu_ls_al_mie=Snu_lb_al_mie; Snv_lb_al_mie=Snu_lb_al_mie; Snv_ls_al_mie=Snu_lb_al_mie;
    for idx=1:floor(vmax/va)
        Snu_lb_al(:,:,1:npts)=Snu_lb_al(:,:,1:npts)+Snu_lb(:,:,1+(idx-1)*npts:idx*npts);
        Snu_ls_al(:,:,1:npts)=Snu_ls_al(:,:,1:npts)+Snu_ls(:,:,1+(idx-1)*npts:idx*npts);
        Snv_lb_al(:,:,1:npts)=Snv_lb_al(:,:,1:npts)+Snv_lb(:,:,1+(idx-1)*npts:idx*npts);
        Snv_ls_al(:,:,1:npts)=Snv_ls_al(:,:,1:npts)+Snv_ls(:,:,1+(idx-1)*npts:idx*npts);
        Snu_lb_al_ray(:,:,1:npts)=Snu_lb_al_ray(:,:,1:npts)+Snu_lb_ray(:,:,1+(idx-1)*npts:idx*npts);
        Snu_ls_al_ray(:,:,1:npts)=Snu_ls_al_ray(:,:,1:npts)+Snu_ls_ray(:,:,1+(idx-1)*npts:idx*npts);
        Snv_lb_al_ray(:,:,1:npts)=Snv_lb_al_ray(:,:,1:npts)+Snv_lb_ray(:,:,1+(idx-1)*npts:idx*npts);
        Snv_ls_al_ray(:,:,1:npts)=Snv_ls_al_ray(:,:,1:npts)+Snv_ls_ray(:,:,1+(idx-1)*npts:idx*npts);
        Snu_lb_al_mie(:,:,1:npts)=Snu_lb_al_mie(:,:,1:npts)+Snu_lb_mie(:,:,1+(idx-1)*npts:idx*npts);
        Snu_ls_al_mie(:,:,1:npts)=Snu_ls_al_mie(:,:,1:npts)+Snu_ls_mie(:,:,1+(idx-1)*npts:idx*npts);
        Snv_lb_al_mie(:,:,1:npts)=Snv_lb_al_mie(:,:,1:npts)+Snv_lb_mie(:,:,1+(idx-1)*npts:idx*npts);
        Snv_ls_al_mie(:,:,1:npts)=Snv_ls_al_mie(:,:,1:npts)+Snv_ls_mie(:,:,1+(idx-1)*npts:idx*npts);
    end

    snr_thres=-1000;
    v_tmp(1,1,1:max(size(vr)))=vr;
    vrep=repmat(v_tmp,[size(Snu_lb,1) size(Snu_lb,2) 1]);
    clear v_tmp;
    v_tmp(1,1,1:max(size(v))-1)=v(1:max(size(v))-1);
    varep=repmat(v_tmp,[size(Snu_lb,1) size(Snu_lb,2) 1]);

    %Compute velocity
    dwsu_al=Snu_ls./Snu_lb;
    yn=isfinite(dwsu_al)|10*log10(dwsu_al)>dws_thres; %&10*log10(Snu_ls)>snr_thres;
    u_est=nansum(Snu_ls.*dwsu_al.*vrep.*yn,3)./nansum(Snu_ls.*dwsu_al.*yn,3);
    dwsv_al=Snv_ls./Snv_lb;
    yn=isfinite(dwsv_al)|10*log10(dwsv_al)>dws_thres; %&10*log10(Snv_ls)>snr_thres;
    v_est=nansum(Snv_ls.*dwsv_al.*vrep.*yn,3)./nansum(Snv_ls.*dwsv_al.*yn,3);
    % clear dwsu_al dwsv_al;

    %Calculate aliased version
    dwsu_alf=Snu_lb_al./Snu_ls_al;
    % yn=isfinite(dwsu_alf)|10*log10(dwsu_alf)>dws_thres; %&10*log10(Snu_ls)>snr_thres;
    % u_estf=nansum(Snu_ls_al.*dwsu_alf.*varep.*yn,3)./nansum(Snu_ls_al.*dwsu_alf.*yn,3);
    dwsv_alf=Snv_lb_al./Snv_ls_al;
    % yn=isfinite(dwsv_alf)|10*log10(dwsv_alf)>dws_thres; %&10*log10(Snv_ls)>snr_thres;
    % v_estf=nansum(Snv_ls_al.*dwsv_alf.*varep.*yn,3)./nansum(Snv_ls_al.*dwsv_alf.*yn,3);

    figure(101)
    subplot(2,2,1)
    pcolor(rtmp2,ztmp2,u_est)
    shading flat
    colorbar
    caxis([-50 50])
    subplot(2,2,2)
    pcolor(rtmp2,ztmp2,uds-urs)
    shading flat
    colorbar
    caxis([-50 50])
    subplot(2,2,3)
    pcolor(rtmp2,ztmp2,v_est)
    shading flat
    colorbar
    caxis([0 70])
    subplot(2,2,4)
    pcolor(rtmp2,ztmp2,vds-vrs)
    shading flat
    colorbar
    caxis([0 70])

    %Try weighting velocity uniformally for dwr > dwr_thres
    dwr_thres_rmse=0;
    yn=nanmax(dwsu_al,[],3)>dwr_thres_rmse;
    rain_ind=5;
    uray=nanmean(udt_avg_sm(:,:,rain_ind),3); vray=nanmean(vdt_avg_sm(:,:,rain_ind),3);
    rms_vel.uest_true=sqrt(nanmean(nanmean(((uds(yn)-urs(yn))-u_est(yn)).^2)));
    rms_vel.ux_true=sqrt(nanmean(nanmean(((uds(yn)-urs(yn))-uds_tmp(yn)).^2)));
    rms_vel.uest_ray=sqrt(nanmean(nanmean((uray(yn)-u_est(yn)).^2)));
    rms_vel.ux_ray=sqrt(nanmean(nanmean((uray(yn)-uds_tmp(yn)).^2)));
    yn=nanmax(dwsv_al,[],3)>dwr_thres_rmse;
    rms_vel.vest_true=sqrt(nanmean(nanmean(((vds(yn)-vrs(yn))-v_est(yn)).^2)));
    rms_vel.vx_true=sqrt(nanmean(nanmean(((vds(yn)-vrs(yn))-vds_tmp(yn)).^2)));
    rms_vel.vest_ray=sqrt(nanmean(nanmean((vray(yn)-v_est(yn)).^2)));
    rms_vel.vx_ray=sqrt(nanmean(nanmean((vray(yn)-vds_tmp(yn)).^2)));

    % rmse_su=sqrt(nanmean(nanmean(((uds-udw)-urs).^2)));
    ridx=4; zidx=3;
    %ridx=8; zidx=10;
    figure(1055)
    plot(vr,squeeze(10*log10(Snu_lb(ridx,zidx,:))),'-k',vr,squeeze(10*log10(Snu_ls(ridx,zidx,:))),'-b'...
        ,vr,squeeze(10*log10(Snu_lb(ridx,zidx,:)))-squeeze(10*log10(Snu_ls(ridx,zidx,:))),'--k')
    title('S(u)')
    legend(['\lambda = ' num2str(1000*lb) ' mm'],['\lambda = ' num2str(1000*ls) ' mm'],'DWR','Location','Best')
    % 
    figure(1056)
    plot(vr,squeeze(10*log10(Snv_lb(ridx,zidx,:))),'-k',vr,squeeze(10*log10(Snv_ls(ridx,zidx,:))),'-b'...
        ,vr,squeeze(10*log10(Snv_lb(ridx,zidx,:)))-squeeze(10*log10(Snv_ls(ridx,zidx,:))),'--k')
    title('S(v)')
    legend(['\lambda = ' num2str(1000*lb) ' mm'],['\lambda = ' num2str(1000*ls) ' mm'],'DWR','Location','Best')

    figure(1057)
    plot(v(yn(ridx,zidx,:)),squeeze(10*log10(dwsu_alf(ridx,zidx,:))),'-k',...
    v(yn(ridx,zidx,:)),squeeze(10*log10(Snu_lb_al(ridx,zidx,:))),'--k',...
    v(yn(ridx,zidx,:)),squeeze(10*log10(Snu_ls_al(ridx,zidx,:))),'-b')
    title('S(u)')
    legend('DWR',['\lambda = ' num2str(1000*lb) ' mm'],['\lambda = ' num2str(1000*ls) ' mm'],'Location','Best')
    figure(1058)
    plot(v(yn(ridx,zidx,:)),squeeze(10*log10(dwsv_alf(ridx,zidx,:))),'-k',...
    v(yn(ridx,zidx,:)),squeeze(10*log10(Snv_lb_al(ridx,zidx,:))),'--k',...
    v(yn(ridx,zidx,:)),squeeze(10*log10(Snv_ls_al(ridx,zidx,:))),'-b')
    title('S(v)')
    legend('DWR',['\lambda = ' num2str(1000*lb) ' mm'],['\lambda = ' num2str(1000*ls) ' mm'],'Location','Best')

    % 
    % v=v(1:360);
    % figure(1057)
    % plot(v,squeeze(10*log10(Snu_lb_al(ridx,zidx,:))),'-k',v,squeeze(10*log10(Snu_lb_al_ray(ridx,zidx,:))),'-b'...
    %     ,v,squeeze(10*log10(Snu_lb_al_mie(ridx,zidx,:))),'-r',v,squeeze(10*log10(Snu_lb_al(ridx,zidx,:)))-squeeze(10*log10(Snu_ls_al(ridx,zidx,:))),'--k')
    % title('Aliased S(u)')
    % legend(['\lambda = ' num2str(1000*lb) ' mm'],['\lambda_{ray} = ' num2str(1000*lb) ' mm'],['\lambda_{mie} = ' num2str(1000*lb) ' mm'],'DWR','Location','Best')
    % 
    % figure(1058)
    % plot(v,squeeze(10*log10(Snv_lb_al(ridx,zidx,:))),'-k',v,squeeze(10*log10(Snv_lb_al_ray(ridx,zidx,:))),'-b'...
    %     ,v,squeeze(10*log10(Snv_lb_al_mie(ridx,zidx,:))),'-r',v,squeeze(10*log10(Snv_lb_al(ridx,zidx,:)))-squeeze(10*log10(Snv_ls_al(ridx,zidx,:))),'--k')
    % title('Aliased S(v)')
    % legend(['\lambda = ' num2str(1000*lb) ' mm'],['\lambda_{ray} = ' num2str(1000*lb) ' mm'],['\lambda_{mie} = ' num2str(1000*lb) ' mm'],'DWR','Location','Best')
    % 
    % figure(1059)
    % plot(v,squeeze(10*log10(Snu_ls_al(ridx,zidx,:))),'-k',v,squeeze(10*log10(Snu_ls_al_ray(ridx,zidx,:))),'-b'...
    %     ,v,squeeze(10*log10(Snu_ls_al_mie(ridx,zidx,:))),'-r',v,squeeze(10*log10(Snu_lb_al(ridx,zidx,:)))-squeeze(10*log10(Snu_ls_al(ridx,zidx,:))),'--k')
    % title('Aliased S(u)')
    % legend(['\lambda = ' num2str(1000*ls) ' mm'],['\lambda_{ray} = ' num2str(1000*ls) ' mm'],['\lambda_{mie} = ' num2str(1000*ls) ' mm'],'DWR','Location','Best')
    % 
    % figure(1060)
    % plot(v,squeeze(10*log10(Snv_ls_al(ridx,zidx,:))),'-k',v,squeeze(10*log10(Snv_ls_al_ray(ridx,zidx,:))),'-b'...
    %     ,v,squeeze(10*log10(Snv_ls_al_mie(ridx,zidx,:))),'-r',v,squeeze(10*log10(Snv_lb_al(ridx,zidx,:)))-squeeze(10*log10(Snv_ls_al(ridx,zidx,:))),'--k')
    % title('Aliased S(v)')
    % legend(['\lambda = ' num2str(1000*ls) ' mm'],['\lambda_{ray} = ' num2str(1000*ls) ' mm'],['\lambda_{mie} = ' num2str(1000*ls) ' mm'],'DWR','Location','Best')

    % %Velocity retrieval based on DWR being positive or negative

    cd(stats_dir)
        save(['final_' fname],'swucoeff','swvcoeff','xwucoeff','xwvcoeff','zhd_sw_coeff','zhd_sx_coeff','zhd_xw_coeff',...
        'rmse_cu','rmse_cv','rmse_kau','rmse_kav','rmse_su','rmse_sv','rmse_xu','rmse_xv',...
        'urs_mean','urc_mean','urx_mean','urka_mean','urw_mean','urs_max','urc_max','urx_max','urka_max','urw_max',...
        'vrs_min','vrc_min','vrx_min','vrka_min','vrw_min','urs_90','urs_50','urx_90','urx_50',...
        'urw_90','urw_50','ddu_sw_90','ddu_sw_50','ddu_xw_90','ddu_xw_50','lb','ls','dwr_mean','rms_vel','rayz','miez');
    cd(base_dir)

    else
        cd(stats_dir)
        save(['final_' fname],'swucoeff','swvcoeff','xwucoeff','xwvcoeff','zhd_sw_coeff','zhd_sx_coeff','zhd_xw_coeff',...
        'rmse_cu','rmse_cv','rmse_kau','rmse_kav','rmse_su','rmse_sv','rmse_xu','rmse_xv',...
        'urs_mean','urc_mean','urx_mean','urka_mean','urw_mean','urs_max','urc_max','urx_max','urka_max','urw_max',...
        'vrs_min','vrc_min','vrx_min','vrka_min','vrw_min','urs_90','urs_50','urx_90','urx_50',...
        'urw_90','urw_50','ddu_sw_90','ddu_sw_50','ddu_xw_90','ddu_xw_50','lb','ls','dwr_mean','rayz','miez');
        cd(base_dir)
    end

    clear ddu_sw ddu_sw_50 duu_sw_90 ddu_sw_tmp ddu_xw ddu_xw_50 ddu_xw_90 ddu_xw_tmp dwr dwr_mean dwr_noise...
        ind lamt1 mieu miev miez nlb nlb1 nlb2 nls nls1 nls2 radii rain_ind rayu rayv rayz rind rtmp2 sigma1 swdu swdv swucoeff...
        swvcoeff tmp tmp_array u_act u_adjust udb_tmp udc udc_ind udc_rad udc_radf udcc udka udka_ind udka_rad udka_radf udkac uds...
        uds_rad uds_ind uds_radf uds_radt uds_tmp udsc udt_avg udt_avg_sm udw udw_ind udw_rad udw_radf udwc udwt udx udx_ind...
        udx_rad udx_radf udx_radt udxc ulb uls um_app ur_app urc urka urs urs_tmp urst urt_avg urt_avg_sm urw_tmp urwt urx urx_tmp...
        v_act vdb_tmp vdc_ind vdc vdc_rad vdka vdka_ind vdka_rad vds vds_ind vds_rad vds_tmp vdt_avg vdt_avg_sm vdw vdw_ind vdw_rad...
        vdx vdx_ind vdx_rad vlb vls vm_app vr_app vrc vrka vrs vrst vrt_avg vrt_avg_sm vrw vrwt vrx w_int_c w_int_cc w_int_ka w_int_kac...
        w_int_s w_int_sc w_int_t w_int_tc w_int_w w_int_wc w_int_x w_int_xc wdc wdc_ind wdc_rad wdka wdka_ind wdka_rad wds wds_ind...
        wds_rad wdt_avg wdt_avg_sm wdw wdw_ind wdw_rad wdx wdx_ind wdx_rad wrc wrka wrs wrt wrt_avg wrt_avg_sm wrw wrx xwdu xwdv...
        xwucoeff xwvcoeff yn2 yn Ze_dry Ze_wet Ze_wood zhb_sm zhb_sm_all zhb_sm_ind zhb_sm_orig zhb_tot zhb_tot_ind zhb_tot_orig zhd_sw...
        zhd_sw_coeff zhd_sx zhd_sx_coeff zhd_xw zhd_xw_coeff zhb_tot zhh_c_store zhh_c_tot zhh_ka_store zhh_ka_tot zhh_s_store zhh_s_tot...
        zhh_w_store zhh_w_tot zhh_x_store zhh_x_tot zhhc_sm zhhc_smt zhhka_sm zhhka_smt zhhs_sm zhhs_smt zhhw_sm zhhw_smt zhhx_sm zhhx_smt...
        zhs_sm zhs_sm_all zhs_sm_ind zhs_sm_orig zhs_tot zhs_tot_all zhs_tot_ind zhs_tot_orig zind zlb zlb1 zlb1_log...
        zlb2 zlb2_divide_zlb1 zlb2_log zls1 zls1_log zls2 zls2_log zmax ztmp2 zhb_tot_all urw dbrs_avg dbrs_cnt_small;


end




