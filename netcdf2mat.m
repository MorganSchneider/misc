%This program converts netcdf files from the NCDC converter to MATLAB
%format, which can later be converted to sweep.  The program handles both
%super res and regular data.  

%1. Get files from NCDC.  Convert them to netcdf
%2. Run this script for super res data flag==0 and regular data
%flag==1

clear;
close all;

% You'll need to edit the code to match your own folder names
% in line 37 (and potentially lines 38 and 39 if you organized your folders
% differently from me), and also line 65 if you decide not to put
% your .mat output files into a subfolder called 'mat'
%
% For reference, my folders are structured like this:
% > Users/schneider
%   > Documents
%       > nexrad
%           > (radar station name, e.g. KTLX, KGWX, etc)
%               > (date of case study, YYYYMMDD format)
%                   > mat
%                       > *this is where I save the .mat files generated
%                           by this script.
%                   > *this is where I put all the NCDC Data Viewer files
%                   (.ar2v) from the event. This folder also contains the 
%                   corresponding netCDF files (.nc) which you export from
%                   the Data Viewer.
%
% So for my case study, the full path to the folder containing my input (.nc) files is:
%   ~/Documents/nexrad/KGWX/20190414/
% and the full path to the folder where my output (.mat) files will go is:
%   ~/Documents/nexrad/KGWX/20190414/mat/
% '~' is shorthand for your home directory.  For me: '~' = 'Users/schneider'

case_dir = uigetdir('~/Documents/nexrad/'); %choose the folder containing your netCDF files
radar_name = case_dir(end-12:end-9); %pull the radar call sign from the path name
case_date = case_dir(end-7:end); %pull the case study date from the path name

flag = 1; %flag = 0 to process super res data, flag = 1 to process regular data
dual_pol_flag=1; %1 for dual-pol, 0 for only single-pol moments

%%

cd(case_dir)
    files=dir('*.nc');
    %file_tmp1='';
    for idx=1:max(size(files))
%     for idx=25:25
        file_tmp=files(idx).name;
        tmp_len=max(size(file_tmp));
        file_time=file_tmp(tmp_len-8:tmp_len-3);
        % **I decided these lines weren't necessary lol whoops**
%         case_date=file_tmp(tmp_len-17:tmp_len-10);
%         if(str2num(date_save)<20000000)
%             keyboard;
%         end
%         if(str2num(file_time)>300000)
%             keyboard;
%         end
%         if(strcmp(file_tmp(1:4),'NOP4'))
%             output_dir=[case_dir '/'];
%             if(~exist(output_dir,'dir'))
%                 mkdir(output_dir)
%             end
%         end
        output_dir=[case_dir '/mat/']; % where to save output files
        if(~exist(output_dir,'dir')) % check if output folder does not exist
            mkdir(output_dir) %if true (i.e. if output folder does not exist), create it
        end
        %If they are the same, it's the same tilt
%         ncid1=netcdf.open(files(idx).name,'NC_WRITE')
%        cd(case_dir)
        ncid1=netcdf.open(files(idx).name,'NC_NOWRITE');
        
        if(round(idx/10)*10==idx)
            fprintf(['Working on file: ' num2str(idx) ' of ' num2str(max(size(files))) '\n']);
        end
        [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(ncid1);
        vname=cell(1,nvars); numatts=cell(1,nvars); type=cell(1,nvars);
        for jdx=1:nvars
            [varname,xtype,dimids,natts]=netcdf.inqVar(ncid1,jdx-1);
            vname{jdx}=varname; numatts{jdx}=natts; type{jdx}=xtype;
        end
%         keyboard;
        gattname=cell(1,ngatts);
        for jdx=1:ngatts
            gattname{jdx} = netcdf.inqAttName(ncid1,netcdf.getConstant('NC_GLOBAL'),jdx-1);
        end
        [ind]=find(strcmp(gattname,'Station'));
        data.Station = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'StationLatitude'));
        data.lat = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'StationLongitude'));
        data.lon = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'StationElevationInMeters'));
        data.elevation = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'base_date'));
        data.base_date = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'time_coverage_start'));
        data.time_coverage_start = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});
        [ind]=find(strcmp(gattname,'time_coverage_end'));
        data.time_coverage_end = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),gattname{ind});

        if(flag==0)
            [ind]=find(strcmp(vname,'Reflectivity_HI'));
            Z_Hs=netcdf.getVar(ncid1,ind-1);
            Z_Hs=int16(Z_Hs);
            tic;

            Z_Hs(Z_Hs<0)=Z_Hs(Z_Hs<0)+256;
            t1=toc;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
%                 fprintf(['Att name: ' attname '\n'])
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'_FillValue'))
                    fill_Value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'units'))
                    units=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            Z_Hs=single(Z_Hs); 
            for kdx=1:max(size(missing_value))
                Z_Hs(Z_Hs==missing_value(kdx))=NaN;
            end
            Z_Hs=Z_Hs*scale_factor+add_offset;
    %         figure(1)
    %         pcolor(squeeze(double(Z_Hs(:,:,1))))
    %         shading flat
    %         colorbar

            [ind]=find(strcmp(vname,'RadialVelocity_HI'));
            vr_hs=netcdf.getVar(ncid1,ind-1);
            vr_hs=int16(vr_hs);
            vr_hs(vr_hs<0)=vr_hs(vr_hs<0)+256;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            vr_hs=single(vr_hs); 
            for kdx=1:max(size(missing_value))
                vr_hs(vr_hs<2)=NaN;
            end
            vr_hs=vr_hs*scale_factor+add_offset;

            [ind]=find(strcmp(vname,'SpectrumWidth_HI'));
            SW_hs=netcdf.getVar(ncid1,ind-1);
            SW_hs=int16(SW_hs);
            SW_hs(SW_hs<0)=SW_hs(SW_hs<0)+256;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            SW_hs=single(SW_hs);
            SW_hs(SW_hs<2)=NaN;
            SW_hs=SW_hs*scale_factor+add_offset;

            if(dual_pol_flag)
                [ind]=find(strcmp(vname,'DifferentialReflectivity_HI'));
                Z_DRs=netcdf.getVar(ncid1,ind-1);
                Z_DRs=int16(Z_DRs);
                Z_DRs(Z_DRs<0)=Z_DRs(Z_DRs<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
        %         Z_DRs(Z_DRs==missing_value)=-9999;
                Z_DRs=single(Z_DRs);
                Z_DRs(Z_DRs<2)=NaN;
                Z_DRs=Z_DRs*scale_factor+add_offset;

                [ind]=find(strcmp(vname,'CorrelationCoefficient_HI'));
                rho_hvs=netcdf.getVar(ncid1,ind-1);
                rho_hvs=int16(rho_hvs);
                rho_hvs(rho_hvs<0)=rho_hvs(rho_hvs<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
                rho_hvs=single(rho_hvs); 
                rho_hvs(rho_hvs<2)=NaN;
        %         for kdx=1:max(size(missing_value))
        %             rho_hvs(rho_hvs==int16(missing_value(kdx)))=NaN;
        %         end
                rho_hvs=rho_hvs*scale_factor+add_offset;

                [ind]=find(strcmp(vname,'DifferentialPhase_HI'));
                Phi_DPs=netcdf.getVar(ncid1,ind-1);
                Phi_DPs=int16(Phi_DPs);
                Phi_DPs(Phi_DPs<0)=Phi_DPs(Phi_DPs<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
                Phi_DPs=single(Phi_DPs);
        %         for kdx=1:max(size(missing_value))
        %             Phi_DPs(Phi_DPs==int16(missing_value(kdx)))=NaN;
        %         end
                Phi_DPs(Phi_DPs<2)=NaN;
                Phi_DPs=Phi_DPs*scale_factor+add_offset;
                [ind]=find(strcmp(vname,'azimuthD_HI'));
                if(min(size(ind))>0)
                    azs=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'elevationD_HI'));
                if(min(size(ind))>0)
                    els=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'timeD_HI'));
                if(min(size(ind))>0)
                    times=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'DistanceD_HI'));
                if(min(size(ind))>0)
                    first_gate=netcdf.getVar(ncid1,ind-1);
                    first_gates=first_gate(1);
                end
            
            else
                [ind]=find(strcmp(vname,'azimuthV_HI'));
                if(min(size(ind))>0)
                    azs=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'elevationV_HI'));
                if(min(size(ind))>0)
                    els=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'timeV_HI'));
                if(min(size(ind))>0)
                    times=netcdf.getVar(ncid1,ind-1);
                end

                [ind]=find(strcmp(vname,'distanceV_HI'));
                if(min(size(ind))>0)
                    first_gate=netcdf.getVar(ncid1,ind-1);
                    first_gates=first_gate(1);
                end
            end
            
            %Get azimuths for super res mode
            [ind]=find(strcmp(vname,'azimuthR_HI'));
            if(min(size(ind))>0)
                azsR=netcdf.getVar(ncid1,ind-1);
            end
            
            %USE THIS TO FIX SUPER RES AZIMUTHS
            [ind]=find(strcmp(vname,'azimuthV_HI'));
            if(min(size(ind))>0)
                azsV=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'elevationR_HI'));
            if(min(size(ind))>0)
                elsR=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'timeR_HI'));
            if(min(size(ind))>0)
                timesR=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'DistanceR_HI'));
            if(min(size(ind))>0)
                first_gate=netcdf.getVar(ncid1,ind-1);
                first_gatesR=first_gate(1);
            end
        end
            
        if(flag==1)
            [ind]=find(strcmp(vname,'distanceR'));
            r_zh=netcdf.getVar(ncid1,ind-1);
            [ind]=find(strcmp(vname,'distanceV'));
            r_vr=netcdf.getVar(ncid1,ind-1);
            [ind]=find(strcmp(vname,'Reflectivity'));
            Z_H=netcdf.getVar(ncid1,ind-1);
            Z_H=int16(Z_H);
            tic;

            Z_H(Z_H<0)=Z_H(Z_H<0)+256;
            t1=toc;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'_FillValue'))
                    fill_Value=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            Z_H=single(Z_H); 
            for kdx=1:max(size(missing_value))
                Z_H(Z_H==missing_value(kdx))=NaN;
            end
            Z_H=Z_H*scale_factor+add_offset;
    %         figure(1)
    %         pcolor(squeeze(double(Z_H(:,:,1))))
    %         shading flat
    %         colorbar

            [ind]=find(strcmp(vname,'RadialVelocity'));
            vr_h=netcdf.getVar(ncid1,ind-1);
            vr_h=int16(vr_h);
            vr_h(vr_h<0)=vr_h(vr_h<0)+256;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            vr_h=single(vr_h); 
            for kdx=1:max(size(missing_value))
                vr_h(vr_h<2)=NaN;
            end
            vr_h=vr_h*scale_factor+add_offset;

            [ind]=find(strcmp(vname,'SpectrumWidth'));
            SW_h=netcdf.getVar(ncid1,ind-1);
            SW_h=int16(SW_h);
            SW_h(SW_h<0)=SW_h(SW_h<0)+256;
            for jdx=1:numatts{ind}
                attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                if(strcmp(attname,'scale_factor'))
                    scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'missing_value'))
                    missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                elseif(strcmp(attname,'add_offset'))
                    add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                end
            end
            SW_h=single(SW_h);
            SW_h(SW_h<2)=NaN;
            SW_h=SW_h*scale_factor+add_offset;
            
            if(dual_pol_flag)
                [ind]=find(strcmp(vname,'DifferentialReflectivity'));
                Z_DR=netcdf.getVar(ncid1,ind-1);
                Z_DR=int16(Z_DR);
                Z_DR(Z_DR<0)=Z_DR(Z_DR<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
        %         Z_DR(Z_DR==missing_value)=-9999;
                Z_DR=single(Z_DR);
                Z_DR(Z_DR<2)=NaN;
                Z_DR=Z_DR*scale_factor+add_offset;

                [ind]=find(strcmp(vname,'CorrelationCoefficient'));
                rho_hv=netcdf.getVar(ncid1,ind-1);
                rho_hv=int16(rho_hv);
                rho_hv(rho_hv<0)=rho_hv(rho_hv<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
                rho_hv=single(rho_hv); 
                rho_hv(rho_hv<2)=NaN;
        %         for kdx=1:max(size(missing_value))
        %             rho_hv(rho_hv==int16(missing_value(kdx)))=NaN;
        %         end
                rho_hv=rho_hv*scale_factor+add_offset;
                [ind]=find(strcmp(vname,'DifferentialPhase'));
                Phi_DP=netcdf.getVar(ncid1,ind-1);
                Phi_DP=int16(Phi_DP);
                Phi_DP(Phi_DP<0)=Phi_DP(Phi_DP<0)+256;
                for jdx=1:numatts{ind}
                    attname=netcdf.inqAttName(ncid1,ind-1,jdx-1);
                    if(strcmp(attname,'scale_factor'))
                        scale_factor=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'missing_value'))
                        missing_value=netcdf.getAtt(ncid1,ind-1,attname);
                    elseif(strcmp(attname,'add_offset'))
                        add_offset=netcdf.getAtt(ncid1,ind-1,attname);
                    end
                end
                Phi_DP=single(Phi_DP);
        %         for kdx=1:max(size(missing_value))
        %             Phi_DP(Phi_DP==int16(missing_value(kdx)))=NaN;
        %         end
                Phi_DP(Phi_DP<2)=NaN;
                Phi_DP=Phi_DP*scale_factor+add_offset;
            end
            [ind]=find(strcmp(vname,'azimuthR'));
            if(min(size(ind))>0)
                az=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'elevationR'));
            if(min(size(ind))>0)
                el=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'timeR'));
            if(min(size(ind))>0)
                time=netcdf.getVar(ncid1,ind-1);
            end

            [ind]=find(strcmp(vname,'DistanceR'));
            if(min(size(ind))>0)
                first_gate=netcdf.getVar(ncid1,ind-1);
                first_gate=first_gate(1);
            end
        end
        
%         [ind]=find(strcmp(vname,'numGatesR'));
%         num_gates=netcdf.getVar(ncid1,ind-1);
%         figure(2)
%         pcolor(squeeze(double(vr_h(:,:,1))))
%         shading flat
%         colorbar
%         figure(3)
%         pcolor(squeeze(double(rho_hv(:,:,1))))
%         shading flat
%         colorbar
%         figure(4)
%         pcolor(squeeze(double(Phi_DP(:,:,1))))
%         shading flat
%         colorbar
%         figure(5)
%         pcolor(squeeze(double(Z_DR(:,:,1))))
%         shading flat
%         colorbar
        
        netcdf.close(ncid1);
        %cd(case_dir)
        if(flag==0)
            for jdx=1:size(azs,2)
                for kdx=1:size(azsR,2)
                    diff=azsR(1,kdx)-azs(1,jdx);
                    if(abs(diff)<0.001)
                    	Zind=kdx;
                    end
                end
                data.Z_H=squeeze(Z_Hs(:,:,Zind));
                data.vr_h=squeeze(vr_hs(:,:,jdx));
                data.SW_h=squeeze(SW_hs(:,:,jdx));
                if(dual_pol_flag)
                    data.Z_DR=squeeze(Z_DRs(:,:,jdx));
                    data.rho_hv=squeeze(rho_hvs(:,:,jdx));
                    data.phi_dp=squeeze(Phi_DPs(:,:,jdx));
                end
                data.az=squeeze(azs(:,jdx));
                data.el=squeeze(els(:,jdx));
                data.time=squeeze(times(:,jdx));
                data.az_ZH=squeeze(azsR(:,Zind));
                data.az_vr=squeeze(azsV(:,jdx));
                if(exist('first_gates','var'))
                    data.first_gate=first_gates;
                end
%                 if(~isfield(data,'Station'))
%                     data.Station = radar_name;
%                 end
                %if(mean(data.el)<10)
                if(jdx<10)
                    save([output_dir case_date '_' file_time '_0' num2str(jdx) '.mat'], 'data');
                    %save([output_dir case_date '_' file_time '_0' num2str(round(mean(data.el),1)) '.mat'],'data');
                else
                    save([output_dir case_date '_' file_time '_' num2str(jdx) '.mat'], 'data');
                    %save([output_dir case_date '_' file_time '_' num2str(round(mean(data.el),1)) '.mat'],'data');
                end
            end
            clear Z_Hs vr_hs SW_hs Z_DRs rho_hvs Phi_DPs azs azsV els times vname gattname numatts type data;
        end        
        
        if(flag==1)
            data.r_zh=squeeze(r_zh);
            data.r_vr=squeeze(r_vr);
            file_offset = max(size(dir([output_dir case_date '_' file_time '_*.mat'])));
            for jdx=1:size(Z_H,3)
                data.Z_H=squeeze(Z_H(:,:,jdx));
                data.vr_h=squeeze(vr_h(:,:,jdx));
                data.SW_h=squeeze(SW_h(:,:,jdx));
                if(dual_pol_flag)
                    data.Z_DR=squeeze(Z_DR(:,:,jdx));
                    data.rho_hv=squeeze(rho_hv(:,:,jdx));
                    data.phi_dp=squeeze(Phi_DP(:,:,jdx));
                end
                data.az=squeeze(az(:,jdx));
                data.el=squeeze(el(:,jdx));
                data.time=squeeze(time(:,jdx));
                if(exist('first_gate','var'))
                    data.first_gate=first_gate;
                end
%                 if(~isfield(data,'Station'))
%                     data.Station = radar_name;
%                 end
                ele_number=jdx+file_offset;
                if(ele_number<10)
                    save([output_dir case_date '_' file_time '_0' num2str(ele_number) '.mat'], 'data');
                    %save([output_dir case_date '_' file_time '_0' num2str(round(mean(data.el),1)) '.mat'],'data');
                else
                    save([output_dir case_date '_' file_time '_' num2str(ele_number) '.mat'], 'data');
                    %save([output_dir case_date '_' file_time '_' num2str(round(mean(data.el),1)) '.mat'],'data');
                end
                
            end
            clear Z_H vr_h SW_h Z_DR rho_hv Phi_DP az el r_vr r_zh time vname gattname numatts type data;
        end
        
        
        %cd('~')
    end
    
    cd('~')
            
            
        
    