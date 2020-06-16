    while options == 8
        %Calculate range
        range_km = sqrt((sum(lims(1:2))/2)^2 + (sum(lims(3:4))/2)^2);
%         height_km=sqrt((sum(lims(1:2))/2)^2+(sum(lims(3:4))/2)^2)*sin(ele*pi/180)+0.020;
        %Calculate beam height
        height_km = sqrt(range_km^2 + (ae)^2 + 2 * range_km * ae * sin(ele*pi/180)) - ae + 0.020;
        
        %Create menu GUI
        figure(100)
        set(100,'Position',pos6)
        h11 = uicontrol('Style', 'pushbutton', 'String', '+el','Position', [2 42 50 40], 'Callback', 'options=1; uiresume');
        h21 = uicontrol('Style', 'pushbutton', 'String', '-el','Position', [2 2 50 40], 'Callback', 'options=2; uiresume');
        h31 = uicontrol('Style', 'pushbutton', 'String', '+vol','Position', [52 42 50 40], 'Callback', 'options=3; uiresume');
        h41 = uicontrol('Style', 'pushbutton', 'String', '-vol','Position', [52 2 50 40], 'Callback', 'options=4; uiresume');
        h51 = uicontrol('Style', 'pushbutton', 'String', '+zoom','Position', [102 42 50 40], 'Callback', 'options=6; options2=2; uiresume');
        h61 = uicontrol('Style', 'pushbutton', 'String', '-zoom','Position', [102 2 50 40], 'Callback', 'options=6; options2=1; uiresume');
        h71 = uicontrol('Style', 'pushbutton', 'String', '+x','Position', [152 42 50 40], 'Callback', 'options=6; options2=5; uiresume');
        h81 = uicontrol('Style', 'pushbutton', 'String', '-x','Position', [152 2 50 40], 'Callback', 'options=6; options2=6; uiresume');
        h91 = uicontrol('Style', 'pushbutton', 'String', '+y','Position', [202 42 50 40], 'Callback', 'options=6; options2=3; uiresume');
        h101 = uicontrol('Style', 'pushbutton', 'String', '-y','Position', [202 2 50 40], 'Callback', 'options=6; options2=4; uiresume');
        h111 = uicontrol('Style', 'pushbutton', 'String', 'Match Windows','Position', [252 2 100 40], 'Callback', 'options=6; options2=8; uiresume');
        h71 = uicontrol('Style', 'pushbutton', 'String', 'Add Marker','Position', [252 42 100 40], 'Callback', 'options=5; uiresume');
%         h81 = uicontrol('Style', 'pushbutton', 'String', 'Change Window','Position', [202 2 100 50], 'Callback', 'options=6; uiresume');
        h91 = uicontrol('Style', 'pushbutton', 'String', 'Data Cursor','Position', [352 42 100 40], 'Callback', 'options=7; uiresume');
        h101 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [452 2 50 40], 'Callback', 'options=9; uiresume');
%         h91 = uicontrol(100,'Style', 'pushbutton', 'String', 'Chg Left Var','Position', [402 62 100 50], 'Callback', 'options=11; uiresume;');
%         h101 = uicontrol(100,'Style', 'pushbutton', 'String', 'Chg Right Var','Position', [402 2 100 50], 'Callback', 'options=12; uiresume;');
        h111 = uicontrol(100,'Style','pushbutton','String','Save Figs','Position',[352 2 100 40],'Callback','options=13; uiresume;');
%         h121 = uicontrol(100,'Style','pushbutton','String','TDS','Position',[502 62 100 50],'Callback','options=14; uiresume;');
%         h131 = uicontrol(100,'Style','pushbutton','String','Filter','Position',[602 62 100 50],'Callback','options=15; uiresume;');
        h141 = uicontrol(100,'Style','pushbutton','String','Dealias','Position',[452 42 50 40],'Callback','options=16; uiresume;');
%         h141 = uicontrol(100,'Style','pushbutton','String','VAD','Position',[452 42 50 40],'Callback','options=17; uiresume;');
        h151 = uicontrol(100,'Style','Text','String',['Height: ' num2str(round(height_km*1000)/1000) ' km'],'Position',[135 115 130 15]);
        h1001 = uicontrol(100,'Style','Text','String',['Range: ' num2str(round(sqrt(nanmean(lims(1:2))^2+nanmean(lims(3:4))^2)*10)/10) ' km'],'Position',[265 115 100 15]);
        figure(100)
        uiwait;
        close(100)
        
%         figure(101)
%         set(101,'Position',pos6+[0 0 50 0])
%         h11 = uicontrol('Style', 'pushbutton', 'String', 'Zoom Out','Position', [2 62 100 50], 'Callback', 'options2=1; uiresume');
%         h21 = uicontrol('Style', 'pushbutton', 'String', 'Zoom In','Position', [2 2 100 50], 'Callback', 'options2=2; uiresume');
%         h31 = uicontrol('Style', 'pushbutton', 'String', 'Y+','Position', [102 62 100 50], 'Callback', 'options2=3; uiresume');
%         h41 = uicontrol('Style', 'pushbutton', 'String', 'Y-','Position', [102 2 100 50], 'Callback', 'options2=4; uiresume');
%         h51 = uicontrol('Style', 'pushbutton', 'String', 'X+','Position', [202 62 100 50], 'Callback', 'options2=5; uiresume');
%         h61 = uicontrol('Style', 'pushbutton', 'String', 'X-','Position', [202 2 100 50], 'Callback', 'options2=6; uiresume');
%         h71 = uicontrol('Style', 'pushbutton', 'String', 'Chg Zoom Factor','Position', [302 62 100 50], 'Callback', 'options2=7; uiresume');
%         h81 = uicontrol('Style', 'pushbutton', 'String', 'Match Windows','Position', [302 2 100 50], 'Callback', 'options2=8; uiresume');
%         h91 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [402 2 50 50], 'Callback', 'options2=9; uiresume');
%         figure(101)
%         uiwait;
%         zoom_factor_x=(lims(2)-lims(1))/2*0.2;
%         zoom_factor_y=(lims(4)-lims(3))/2*0.2;
%         if(options2==9)
%             close(101);
%         end
      
        if(options==11)
            %Changes variable type for left panel
            if(plot_choice(1)==5)
                figure(1)
                colorbar off
            end
            new_color_flag=[1 0];
            figure(102)
            set(102,'Position',pos7)
    %         text(pos7(1),pos7(2)/2,'Left Panel Control')
            h12 = uicontrol(102,'Style', 'pushbutton', 'String', 'Z','Position',[2 2 60 40], 'Callback', 'plot_choice(1)=1; options=10; uiresume;');
            h22 = uicontrol(102,'Style', 'pushbutton', 'String', 'Vr','Position',[2 44 60 40], 'Callback', 'plot_choice(1)=2; options=10; uiresume;');
            h32 = uicontrol(102,'Style', 'pushbutton', 'String', 'rhohv','Position',[64 2 60 40], 'Callback', 'plot_choice(1)=3; options=10; uiresume;');
            h42 = uicontrol(102,'Style', 'pushbutton', 'String', 'ZDR','Position',[64 44 60 40], 'Callback', 'plot_choice(1)=4; options=10; uiresume;');
            h52 = uicontrol(102,'Style', 'pushbutton', 'String', 'Thres','Position',[126 44 60 40], 'Callback', 'plot_choice(1)=5; options=10; uiresume;');
            h62 = uicontrol(102,'Style', 'pushbutton', 'String', 'PhiDP','Position',[126 2 60 40],'Callback','plot_choice(1)=6; options=10; uiresume;');
            h72 = uicontrol(102,'Style', 'pushbutton', 'String', 'KDP','Position',[190 44 60 40],'Callback','plot_choice(1)=7; options=10; uiresume;');
            uiwait;
            close(102)
        elseif(options==12)
            %Changes variable type for right panel
            if(plot_choice(2)==5)
                figure(2)
                colorbar off
            end
            new_color_flag=[0 1];
            figure(103)
            set(103,'Position',pos8)
    %         text(pos8(1),pos8(2)/2,'Right Panel Control')
            h13 = uicontrol(103,'Style', 'pushbutton', 'String', 'Z','Position',[2 2 60 40], 'Callback', 'plot_choice(2)=1; options=10; uiresume;');
            h14 = uicontrol(103,'Style', 'pushbutton', 'String', 'Vr','Position',[2 44 60 40], 'Callback', 'plot_choice(2)=2; options=10; uiresume;');
            h15 = uicontrol(103,'Style', 'pushbutton', 'String', 'rhohv','Position',[64 2 60 40], 'Callback', 'plot_choice(2)=3; options=10; uiresume;');
            h16 = uicontrol(103,'Style', 'pushbutton', 'String', 'ZDR','Position',[64 44 60 40], 'Callback', 'plot_choice(2)=4; options=10; uiresume;');
            h17 = uicontrol(103,'Style', 'pushbutton', 'String', 'Thres','Position',[126 44 60 40], 'Callback', 'plot_choice(2)=5; options=10; uiresume;');
            h18 = uicontrol(103,'Style', 'pushbutton', 'String', 'PhiDP','Position',[126 2 60 40],'Callback','plot_choice(2)=6; options=10; uiresume;');
            h19 = uicontrol(103,'Style', 'pushbutton', 'String', 'KDP','Position',[190 44 60 40],'Callback','plot_choice(2)=7; options=10; uiresume;');
            uiwait;
            close(103)
        elseif(options==13)
            %Saves figures
            title_str={'Z','vr','rhohv','ZDR','thres','Phidp','','','','Phidp_var','Zdr_var'};
            if(isfield(data,'time'))
                dnum=datenum([2019 4 14 0 0 double(data.time(1))/1000]);
            end
            figure(1)
            picname=[datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10) '_' title_str{plot_choice(1)} '.png'];
            eval('cd imgs/')
            print(picname,'-dpng');
            eval('cd ..')
            figure(2)
            picname=[datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10) '_' title_str{plot_choice(2)} '.png'];
            eval('cd imgs/')
            print(picname,'-dpng');
            eval('cd ..')
            if(big_monitor)
                figure(3)
                picname=[datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10) '_' title_str{plot_choice(3)} '.png'];
                eval('cd imgs/')
            print(picname,'-dpng');
            eval('cd ..')
                figure(4)
                picname=[datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10) '_' title_str{plot_choice(4)} '.png'];
                eval('cd imgs/')
            print(picname,'-dpng');
            eval('cd ..')
                figure(5)
                picname=[datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10) '_' title_str{plot_choice(5)} '.png'];
                eval('cd imgs/')
                print(picname,'-dpng');
                eval('cd ..')
            end
        elseif(options==14)
            %This code is used to save the center of a TVS or mesocyclone
            if(ouprime_flag)
                if(~exist('storm_num','var') && strcmp(date_save,'20110524'))
                    storm_num=input(['Which storm? (1: El Reno, 2: Chickasha/Newcastle, 3: Goldsby): ']);
                    while(min(size(storm_num))<1)
                        storm_num=input(['Which storm? (1: El Reno, 2: Chickasha/Newcastle, 3: Goldsby): ']);
                    end
                    save_ctr_file=['20110524_' num2str(storm_num) '_' time_choice '_' num2str(ele_save(scan_num)) '.mat'];
                elseif(strcmp(date_save,'20110524'))
                    save_ctr_file=['20110524_' num2str(storm_num) '_' time_choice '_' num2str(ele_save(scan_num)) '.mat'];
                end
                if(~exist('storm_num','var') && strcmp(date_save,'20100510'))
                    storm_num=input(['Enter tornado number: ']);
                    while(min(size(storm_num))<1)
                        storm_num=input(['Enter tornado number: ']);
                    end
                    save_ctr_file=['20100510_' num2str(storm_num) '_' time_choice '_' num2str(ele_save(scan_num)) '.mat'];
                elseif(strcmp(date_save,'20100510'))
                    save_ctr_file=['20100510_' num2str(storm_num) '_' time_choice '_' num2str(ele_save(scan_num)) '.mat'];
                end
            else
                storm_num=input(['Enter tornado number: ']);
                save_ctr_file=['KOUN_' num2str(file_tmp(1:8)) '_' num2str(storm_num) '_' time_choice '_' num2str(ele_save(scan_num)) '.mat'];
            end

            [tds]=tds_analysis_edge(data,xx,yy,save_ctr_file,ele_save(scan_num),big_monitor,0,save_ctr_file,deb_fall_flag);
        elseif(options==15)
            %These options give some different options for filtering/performing image processing on the
            %data
            options3=1;
            while(options3~=4)
                figure(102)
                set(102,'Position',pos10+[0 0 50 0])
                h11 = uicontrol('Style', 'pushbutton', 'String', 'Clutter Filt.','Position', [2 62 100 50], 'Callback', 'options3=1; uiresume');
                h21 = uicontrol('Style', 'pushbutton', 'String', 'Median Filt.','Position', [2 2 100 50], 'Callback', 'options3=2; uiresume');
                h31 = uicontrol('Style', 'pushbutton', 'String', 'Avg. Filt.','Position', [102 2 100 50], 'Callback', 'options3=3; uiresume');
                h41 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [102 62 100 50], 'Callback', 'options3=4; uiresume');
    %             h31 = uicontrol('Style', 'pushbutton', 'String', 'Y+','Position', [102 62 100 50], 'Callback', 'options3=3; uiresume');
    %             h41 = uicontrol('Style', 'pushbutton', 'String', 'Y-','Position', [102 2 100 50], 'Callback', 'options3=4; uiresume');
    %             h51 = uicontrol('Style', 'pushbutton', 'String', 'X+','Position', [202 62 100 50], 'Callback', 'options3=5; uiresume');
    %             h61 = uicontrol('Style', 'pushbutton', 'String', 'X-','Position', [202 2 100 50], 'Callback', 'options3=6; uiresume');
    %             h71 = uicontrol('Style', 'pushbutton', 'String', 'Chg Zoom Factor','Position', [302 62 100 50], 'Callback', 'options3=7; uiresume');
    %             h81 = uicontrol('Style', 'pushbutton', 'String', 'Match Windows','Position', [302 2 100 50], 'Callback', 'options3=8; uiresume');
    %             h91 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [402 2 50 50], 'Callback', 'options3=9; uiresume');
                figure(102)
                uiwait;
                close(102)
                if(options3==1)
                    if(isfield(data,'Z_corr'))
                        yn6=data.Z_corr>45&abs(data.vr_h<2)&data.rho_hv<0.7;
                        data.Z_corr(yn6)=NaN;
                    else
                        yn6=data.Z_H>45&abs(data.vr_h<2)&data.rho_hv<0.7;
                        data.Z_H(yn6)=NaN;
                    end
                elseif(options3==2 || options3==3)
                    yn7=xx>lims(1)&xx<lims(2)&yy>lims(3)&yy<lims(4);
                    [row3,col3]=find(yn7==1); row_win=1; col_win=1; filter_win=2;
                    if(isfield(data,'Z_corr'))
                        Z_tmp=data.Z_corr;
                        Zdr_tmp=data.ZDR_corr;
                    else
                        Z_tmp=data.Z_H;
                        Zdr_tmp=data.Z_DR;
                    end
%                     vr_tmp=data.vr_h;
                    rhohv_tmp=data.rho_hv;
                    if(isfield(data,'SW'))
                        SW_tmp=data.SW;
                    else
                        SW_tmp=data.SW_h;
                    end
                    if(isfield(data,'phi_dp'))
                        phidp_tmp=data.phi_dp;
                    end
                    for idx=1:max(size(row3))
                        if(filter_win==1)
                            %Apply a n by n square window
                            tmp_ind_row=row3(idx)-row_win:row3(idx)+row_win; 
                            tmp_ind_col=col3(idx)-col_win:col3(idx)+col_win;
                        else
                            %Cross window (uses radially and azimuthally
                            %adjacent gates
                            tmp_ind_row=[row3(idx) row3(idx)-1 row3(idx) row3(idx)+1 row3(idx)];
                            tmp_ind_col=[col3(idx)-1 col3(idx) col3(idx) col3(idx) col3(idx)+1];
                        end
                        if(sum(tmp_ind_row<=0)>0)
                            [row5]=find(tmp_ind_row<=0);
                            tmp_ind_row(row5)=tmp_ind_row(row5)+size(data.rho_hv,1);
                        elseif(sum(tmp_ind_row>size(data.rho_hv,1))>0)
                            [row5]=find(tmp_ind_row>size(data.rho_hv,1));
                            tmp_ind_row(row5)=tmp_ind_row(row5)-size(data.rho_hv,1);
                        end
                        if(sum(tmp_ind_col<=0)>0)
                            [row5]=find(tmp_ind_col<=0);
                            tmp_ind_col(row5)=tmp_ind_col(row5)+size(data.rho_hv,2);
                        end
                        if(isfield(data,'Z_corr')&&att_corr_flag)
                            tmp=Z_tmp(tmp_ind_row,tmp_ind_col);
                            if(options3==2)
                                data.Z_corr(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                            else
                                data.Z_corr(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                            end
                            tmp=Zdr_tmp(tmp_ind_row,tmp_ind_col);
                            if(options3==2)
                                data.ZDR_corr(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                            else
                                data.ZDR_corr(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                            end
                        else
                            tmp=Z_tmp(tmp_ind_row,tmp_ind_col);
                            if(options3==2)
                                data.Z_H(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                            else
                                data.Z_H(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                            end
                            tmp=Zdr_tmp(tmp_ind_row,tmp_ind_col);
                            if(options3==2)
                                data.Z_DR(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                            else
                                data.Z_DR(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                            end
                        end
%                         tmp=vr_tmp(tmp_ind_row,tmp_ind_col);
%                         if(options3==2)
%                             data.vr_h(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
%                         else
%                             data.vr_h(row3(idx),col3(idx))=nanmean(nanmean(tmp));
%                         end
                        tmp=rhohv_tmp(tmp_ind_row,tmp_ind_col);
                        if(options3==2)
                            data.rho_hv(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                        else
                            data.rho_hv(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                        end
                        tmp=SW_tmp(tmp_ind_row,tmp_ind_col);
                        if(options3==2)
                            data.SW(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                        else
                            data.SW(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                        end
%                         tmp=phidp_tmp(tmp_ind_row,tmp_ind_col);
                        if(isfield(data,'phi_dp'))
                            tmp=phidp_tmp(tmp_ind_row,tmp_ind_col);
                            if(options3==2)
                                data.phi_dp(row3(idx),col3(idx))=median(tmp(isfinite(tmp)));
                            else
                                data.phi_dp(row3(idx),col3(idx))=nanmean(nanmean(tmp));
                            end
                        end
                        
                    end
                end
                %Apply filter to the figure
                for pdx=1:num_loops
                    figure(pdx)
                    if(plot_choice(pdx)~=7)
                        set(p1(pdx),'XData',xx)
                        set(p1(pdx),'YData',yy)
                    else
                        set(p1(pdx),'XData',xx2)
                        set(p1(pdx),'YData',yy2)
                    end
                    if(plot_choice(pdx)==1)
                        if(att_corr_flag)
                            set(p1(pdx),'CData',double(data.Z_corr))
                        else
                            set(p1(pdx),'CData',double(data.Z_H))
                        end
%                     elseif(plot_choice(pdx)==2)
%                         set(p1(pdx),'CData',double(data.vr_h))
                    elseif(plot_choice(pdx)==3)
                        set(p1(pdx),'CData',double(data.rho_hv))
                    elseif(plot_choice(pdx)==4)
                        if(att_corr_flag)
                            set(p1(pdx),'CData',double(data.ZDR_corr))
                        else
                            set(p1(pdx),'CData',double(data.Z_DR))
                        end
                    elseif(plot_choice(pdx)==6)
                        set(p1(pdx),'CData',double(data.phi_dp))
                    elseif(plot_choice(pdx)==9)
                        set(p1(pdx),'CData',double(data.SW))
                    elseif(plot_choice(pdx)==10)
                        set(p1(pdx),'CData',double(sd_phidp))
                    elseif(plot_choice(pdx)==11)
                        set(p1(pdx),'CData',double(sd_ZDR))
                    end
                end
            end
            options=8; options3=4; clear tmp vr_tmp SW_tmp rhohv_tmp Zdr_tmp Z_tmp
            
        elseif(options==16)
            %Options for dealiasing data or editting/censoring
            
            %Dealias data
            figure(102)
            set(102,'Position',pos10+[0 0 50 0])
            h11 = uicontrol('Style', 'pushbutton', 'String', 'Dealias Cir.','Position', [2 62 100 50], 'Callback', 'options3=1; uiresume');
            h31 = uicontrol('Style', 'pushbutton', 'String', 'Dealias Rect.','Position', [2 62 100 50], 'Callback', 'options3=5; uiresume');
            h21 = uicontrol('Style', 'pushbutton', 'String', 'Censor Data','Position', [2 2 100 50], 'Callback', 'options3=2; uiresume');
            h41 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [102 62 100 50], 'Callback', 'options3=4; uiresume');
            figure(102)
            uiwait;
            close(102)
            if(options3==1 || options3==5)
                if(~isfield(data,'fold_int_map'))
                    data.fold_int_map=int8(zeros(size(xx)));
                end
                %Dealias commands
                if(isfield(data,'va'))
                    va_flag=input(['Is ' num2str(round(data.va*10)/10) ' the correct va? (1/0): ']); 
                else
                    data.va=nanmax(nanmax(data.vr_h));
                    va_flag=input(['Is ' num2str(round(data.va*10)/10) ' the correct va? (1/0): ']);
                    if(va_flag==0)
                        data.va=input(['Input correct va: ']);
                    end
                end
                %dealias_v1;
                dealias_flag=true;
                while(dealias_flag)
                    %data.va=nanmax(nanmax(data.vr_h));
                    figure(2)
                    if(options3==1)
                        figure(104)
                        set(104,'Position',pos10+[0 0 100 0])
                        h16 = uicontrol('Style', 'pushbutton', 'String', 'R=0.25','Position', [2 62 100 50], 'Callback', 'rad_dealias=0.25; uiresume');
                        h26 = uicontrol('Style', 'pushbutton', 'String', 'R=0.50','Position', [2 2 100 50], 'Callback', 'rad_dealias=0.50; uiresume');
                        h36 = uicontrol('Style', 'pushbutton', 'String', 'R=0.75','Position', [102 2 100 50], 'Callback', 'rad_dealias=0.75; uiresume');
                        h46 = uicontrol('Style', 'pushbutton', 'String', 'R=1','Position', [102 62 100 50], 'Callback', 'rad_dealias=1; uiresume');
                        h56 = uicontrol('Style', 'pushbutton', 'String', 'R=1.5','Position', [202 2 100 50], 'Callback', 'rad_dealias=1.5; uiresume');
                        h66 = uicontrol('Style', 'pushbutton', 'String', 'R=2','Position', [202 62 100 50], 'Callback', 'rad_dealias=2; uiresume');
                        h76 = uicontrol('Style','pushbutton', 'String', 'End','Position', [302 2 100 50], 'Callback', 'rad_dealias=-1; uiresume');
                        uiwait;
                        close(104)
                        if(rad_dealias~=-1)
                            figure(105)
                            set(105,'Position',pos10+[600 0 200 0])
                            h76 = uicontrol('Style', 'pushbutton', 'String', 'Int 1','Position', [2 2 100 50], 'Callback', 'fold_int=1; uiresume');
                            h86 = uicontrol('Style', 'pushbutton', 'String', 'Int 2','Position', [2 62 100 50], 'Callback', 'fold_int=2; uiresume');
                            h96 = uicontrol('Style', 'pushbutton', 'String', 'Int -1','Position', [102 2 100 50], 'Callback', 'fold_int=-1; uiresume');
                            h106 = uicontrol('Style', 'pushbutton', 'String', 'Int -2','Position', [102 62 100 50], 'Callback', 'fold_int=-2; uiresume');
                            h116 = uicontrol('Style', 'pushbutton', 'String', 'Int 0','Position', [202 2 100 50], 'Callback', 'fold_int=0; uiresume');
                            h86 = uicontrol('Style', 'pushbutton', 'String', 'Int 3','Position', [202 62 100 50], 'Callback', 'fold_int=3; uiresume');
                            uiwait;
                            close(105)
                        end
                        if(rad_dealias~=-1)
                            figure(2)
                            [xtmp,ytmp]=ginput;
                        else
                            dealias_flag=false;
                        end

                        while(min(size(rad_dealias))==0)
                            rad_dealias=input('Enter radius: ');
                        end
                        if(rad_dealias~=-1)
                            for ldx=1:max(size(xtmp))
                                if(rad_dealias~=-1)
                                    if(exist('xx3','var'))
                                        pos=sqrt((xx3-xtmp(ldx)).^2+(yy3-ytmp(ldx)).^2);
                                    elseif(exist('xx2','var'))
                                        pos=sqrt((xx2-xtmp(ldx)).^2+(yy2-ytmp(ldx)).^2);
                                    else
                                        pos=sqrt((xx-xtmp(ldx)).^2+(yy-ytmp(ldx)).^2);
                                    end
                                    [row3,col3]=find(pos<rad_dealias);
                                    if(~exist('fold_int','var'))
                                        fold_int=input('Enter unfolding interval (e.g., +va is 1, 2*va is 2, -va is -1: ');
                                    end
                                    if(fold_int==-99)
                                        dealias_flag=false;
                                    end
                                else
                                    dealias_flag=false; row3=[];
                                end
                                if(dealias_flag && min(size(row3))>0)
                                    if(fold_int>0)
                                        for pdx=1:max(size(row3))
                                           if(data.vr_h(row3(pdx),col3(pdx))+double(data.fold_int_map(row3(pdx),col3(pdx)))*data.va<(fold_int-1)*data.va)
                                               data.fold_int_map(row3(pdx),col3(pdx))=int8(data.fold_int_map(row3(pdx),col3(pdx))+2);
                                           end
                                        end
                                    elseif(fold_int<0)
                                        for pdx=1:max(size(row3))
                                            if((data.vr_h(row3(pdx),col3(pdx))+double(data.fold_int_map(row3(pdx),col3(pdx)))*data.va)>(fold_int+1)*data.va)
                                                data.fold_int_map(row3(pdx),col3(pdx))=int8(data.fold_int_map(row3(pdx),col3(pdx))-2);
                                            end
                                        end
                                    else
                                        for pdx=1:max(size(row3))
                                            data.fold_int_map(row3(pdx),col3(pdx))=int8(fold_int);
                                        end
                                    end
                                    figure(2)
                                    set(p1(2),'CData',double(data.vr_h+double(data.fold_int_map)*data.va));
                                    caxis([-double(max(max(abs(data.fold_int_map))))*data.va double(max(max(abs(data.fold_int_map))))*data.va])
                                    clear row3 col3 pos
                                    cd(dir_name)
                                        save(files(row).name,'data')
                                    cd(base_dir)
                                end
                            end
                        end
                    
                    else
                        figure(105)
                        set(105,'Position',pos10+[600 0 0 0])
                        h76 = uicontrol('Style', 'pushbutton', 'String', 'Int 1','Position', [2 2 100 50], 'Callback', 'fold_int=1; uiresume');
                        h86 = uicontrol('Style', 'pushbutton', 'String', 'Int 2','Position', [2 62 100 50], 'Callback', 'fold_int=2; uiresume');
                        h96 = uicontrol('Style', 'pushbutton', 'String', 'Int -1','Position', [102 2 100 50], 'Callback', 'fold_int=-1; uiresume');
                        h106 = uicontrol('Style', 'pushbutton', 'String', 'Int -2','Position', [102 62 100 50], 'Callback', 'fold_int=-2; uiresume');
                        h116 = uicontrol('Style', 'pushbutton', 'String', 'Int 0','Position', [202 2 100 50], 'Callback', 'fold_int=0; uiresume');
                        h86 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [202 62 100 50], 'Callback', 'fold_int=-99; uiresume');
                        uiwait;
                        close(105)
                       if(fold_int==-99)
                           dealias_flag=false;
                       else
                           dealias_flag=true;
                       end
                       if(dealias_flag)
                           [xtmp,ytmp]=ginput(4);
                           %Find closest points in range, azimuth
                           xtmp=sort(xtmp); ytmp=sort(ytmp);
                           if(exist('xx3','var'))
                                yn=xx3>xtmp(2)&xx3<xtmp(3)&yy3>ytmp(2)&yy3<ytmp(3);
                                [row3,col3]=find(yn==1);
                           elseif(exist('xx2','var'))
                                yn=xx2>xtmp(2)&xx2<xtmp(3)&yy2>ytmp(2)&yy2<ytmp(3);
                                [row3,col3]=find(yn==1);
                           else
                                yn=xx>xtmp(2)&xx<xtmp(3)&yy>ytmp(2)&yy<ytmp(3);
                                [row3,col3]=find(yn==1);
                           end
                           if(~exist('fold_int','var'))
                                fold_int=input('Enter unfolding interval (e.g., +va is 1, 2*va is 2, -va is -1: ');
                           end
                       end
                       
                       if(dealias_flag && min(size(row3))>0)
                            if(fold_int>0)
                                for pdx=1:max(size(row3))
                                   if(data.vr_h(row3(pdx),col3(pdx))+double(data.fold_int_map(row3(pdx),col3(pdx)))*data.va<(fold_int-1)*data.va)
                                       data.fold_int_map(row3(pdx),col3(pdx))=int8(data.fold_int_map(row3(pdx),col3(pdx))+2);
                                   end
                                end
                            elseif(fold_int<0)
                                for pdx=1:max(size(row3))
                                    if((data.vr_h(row3(pdx),col3(pdx))+double(data.fold_int_map(row3(pdx),col3(pdx)))*data.va)>(fold_int+1)*data.va)
                                        data.fold_int_map(row3(pdx),col3(pdx))=int8(data.fold_int_map(row3(pdx),col3(pdx))-2);
                                    end
                                end
                            else
                                for pdx=1:max(size(row3))
                                    data.fold_int_map(row3(pdx),col3(pdx))=int8(fold_int);
                                end
                            end
                            figure(2)
                            set(p1(2),'CData',double(data.vr_h+double(data.fold_int_map)*data.va));
                            caxis([-double(max(max(abs(data.fold_int_map))))*data.va double(max(max(abs(data.fold_int_map))))*data.va])
                            clear row3 col3 pos
                            
                            if(data_type==1)
                                cd(dir_name)
                                    save(files(row).name,'data')
                                cd(base_dir)
                            else
                                write_vout;
                            end
                        end
                    end
                end 
                options3=2; options=8;
            end
           
        elseif(options==17)
            %[vad_data]=read_vel_vol(
        elseif(options==1)
            scan_num=scan_num+1; clear data xx2 yy2;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        elseif(options==2)
            scan_num=scan_num-1; clear data xx2 yy2;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        elseif(options==3)
            ele_current=ele_save(scan_num); ele_tmp=ele_save(scan_num); new_scan_num=scan_num+1; ele_current=ele_save(new_scan_num);
            while(ele_tmp~=ele_current)
                new_scan_num=new_scan_num+1;
                ele_current=ele_save(new_scan_num);
                if(new_scan_num>10000)
                    fprintf('Error: no tilt exists in next volume');
                    keyboard;
                end
            end
            scan_num=new_scan_num; clear data xx2 yy2;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        elseif(options==4)
            ele_current=ele_save(scan_num); ele_tmp=ele_save(scan_num); new_scan_num=scan_num-1; ele_current=ele_save(new_scan_num);
            while(ele_tmp~=ele_current)
                new_scan_num=new_scan_num-1;
                ele_current=ele_save(new_scan_num);
                if(abs(new_scan_num)>10000)
                    fprintf('Error: no tilt exists in next volume');
                    keyboard;
                end
            end
            scan_num=new_scan_num; clear data xx2 yy2;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        elseif(options==5)
%            fprintf(['Choose a time: \n']);
%             for idx=1:max(size(file_time_save))
%                 fprintf([num2str(idx) ': ' num2str(roundn(ele_save(idx),-1)) ' deg ' file_time_save{idx}(1:2) ':' file_time_save{idx}(3:4) ':' file_time_save{idx}(5:6) '\n']);
%             end
%             scan_num=input('Choose scan number: ');
%             [row]=find(scan_number==scan_num);
%             time_choice=file_time_save{scan_num}(1:6);
%             clear data;
        latcoords=input('Enter latitude: ');
        loncoords=input('Enter longitude: ');
        marker_name=input('Enter Marker Name','s');
        [x_mark,y_mark]=latlon2cart(latcoords,loncoords,lat,lon);
        fprintf(num2str(num_loops))
        for idx=1:num_loops
            figure(idx)
            hold on
            plot(x_mark,y_mark,'*','MarkerSize',8)
            text(x_mark,y_mark+1,marker_name);
            hold off
        end

        elseif(options==6)
            %Change zoom, x and y centering, and adjust limits
            lims(1:2)=get(h1,'xlim');
            lims(3:4)=get(h1,'ylim');

            while(options2~=9)
                zoom_factor_x=(lims(2)-lims(1))/2*0.2;
                zoom_factor_y=(lims(4)-lims(3))/2*0.2;
                if(options2==1)
                    lims(1)=lims(1)-zoom_factor_x; lims(2)=lims(2)+zoom_factor_x; lims(3)=lims(3)-zoom_factor_y; lims(4)=lims(4)+zoom_factor_y;
                elseif(options2==2)
                    lims(1)=lims(1)+zoom_factor_x; lims(2)=lims(2)-zoom_factor_x; lims(3)=lims(3)+zoom_factor_y; lims(4)=lims(4)-zoom_factor_y;
                elseif(options2==3)
                    lims(3:4)=lims(3:4)+zoom_factor_y;
                elseif(options2==4)
                    lims(3:4)=lims(3:4)-zoom_factor_y;
                elseif(options2==5)
                    lims(1:2)=lims(1:2)+zoom_factor_x;
                elseif(options2==6)
                    lims(1:2)=lims(1:2)-zoom_factor_x;
                elseif(options2==7)
                    zoom_factor_x=input(['Input new zoom factor x (current value: ' num2str(zoom_factor_x) '): ']); 
                    zoom_factor_y=input(['Input new zoom factor y (current value: ' num2str(zoom_factor_y) '): ']); 
                elseif(options2==8)
                    pause; fprintf('Change figure zoom now');
                    if(exist('pos3','var') && big_monitor==0)
                        figure(gcf)
                        fignum=inputdlg('Figure to match limits: 1. Z, 2. Vr, 3. rhv 4. zdr'); fignum=str2double(fignum);
                        if(fignum==1)
                            fignum=h1; fignum1=1; %fignum1=115; 
                        elseif(fignum==2)
                            fignum=h2; fignum1=3;
                        elseif(fignum==3)
                            fignum=h3; fignum1=2; %fignum1=116;
                        else
                            fignum=h4; fignum1=4;
                        end     
                    else
                        figure(gcf)
                        fignum=str2double(inputdlg('Figure 1 or 2?'));
                        if(fignum==1)
                            fignum=h1;
                            fignum1=1;
                        else
                            fignum=h2;
                            fignum1=2;
                        end 
                    end
                    figure(fignum1)
                    lims(1:2)=get(fignum,'xlim');
                    lims(3:4)=get(fignum,'ylim');
                    zoom_factor_x=(lims(2)-lims(1))/2*0.2;
                    zoom_factor_y=(lims(4)-lims(3))/2*0.2;
                end
                if(exist('pos3','var'))
                    figure(1)
                    axis(lims)
                    figure(3)
                    axis(lims)
                    figure(2)
                    axis(lims)
                    figure(4)
                    axis(lims)
                end
                options2=9;
            end
            options=8;
        elseif(options==7)
            %Data cursor: obtain values at a specific point by clicking
            getvalue_flag=1; pdx=1; qdx=1;
            while(getvalue_flag)
                if(pdx==1)
                    pause; pdx=pdx+1;
                end
                if(data_type==1)
                    if(az(5)-az(4)>0)
                        switch_factor=1;
                    else
                        switch_factor=-1;
                    end
                end
                [x1,y1]=ginput(1);
                ran1=sqrt(x1^2+y1^2)*1000; gat1=ceil((ran1-data.detr)/data.detr); azi1=atan2(y1,x1)*180/pi;
                if(data_type==1)
                    azi1=450-azi1-0.25*switch_factor;
                else
                    azi1=450-azi1;
                end
                ran1=ran1/1000;
                height_km=sqrt(ran1^2+(ae)^2+2*ran1*ae*sin(ele*pi/180))-ae+0.020;
                if(azi1>360)
                    azi1=azi1-360;
                end
                
                [c1,azi_id]=min(abs(azi1-az));
                if(~ouprime_flag && data_type==1)
                    gat_tmp=gat1; gat1=azi_id; azi_id=gat_tmp; 
                else
                    [c1,gat1]=min(abs(r-ran1));
                end
                fprintf('\n');
                fprintf(['Range: ' num2str(round(ran1*100)/100) ' km,' ' '])
                fprintf(['Azi.: ' num2str(round(azi1*10)/10) ' deg.,' ' '])
                fprintf(['Height: ' num2str(round(height_km*1000)/1000) ' km' '\n']);
                if(data_type==1)
                    fprintf(['Z: ' num2str(data.Z_H(azi_id,gat1)) ' dBZ, ']);
                    fprintf(['Vr: ' num2str(data.vr_h(azi_id,gat1)) ' m/s, ']);
                    fprintf(['rhv: ' num2str(round(data.rho_hv(azi_id,gat1)*1000)/1000) ', ']);
                    fprintf(['Zdr: ' num2str(data.Z_DR(azi_id,gat1)) ' dB, ' '\n']); 
                else
                    fprintf(['Z: ' num2str(data.Z_H(gat1,azi_id)) ' dBZ, ']);
                    fprintf(['Vr: ' num2str(data.vr_h(gat1,azi_id)) ' m/s, ']);
                    fprintf(['rhv: ' num2str(round(data.rho_hv(gat1,azi_id)*1000)/1000) ', ']);
                    fprintf(['Zdr: ' num2str(data.Z_DR(gat1,azi_id)) ' dB, ' '\n']);
                end
                
                range_profile=false;
                if(range_profile)
                    filt_leng=5;
                    if(~ouprime_flag)
                        max_loop_size=round(azi_id*1.5)-filt_leng;
                    else
                        max_loop_size=round(gat1*1.5)-filt_leng;
                    end
                    for idx=1:max_loop_size
                        r_avg(idx)=nanmean(r(idx:idx+filt_leng));
                        if(~ouprime_flag)

                            phidp_avg(idx)=atan2(nanmean(sind(data.phi_dp(idx:idx+filt_leng,gat1))),nanmean(cosd(data.phi_dp(idx:idx+filt_leng,gat1))))*180/pi;
                        else
                            phidp_avg(idx)=atan2(nanmean(sind(data.phi_dp(azi_id,idx:idx+filt_leng))),nanmean(cosd(data.phi_dp(azi_id,idx:idx+filt_leng))))*180/pi;
                        end
                    end
                    curr_fig=gcf;
                    figure(1003)
                    subplot(2,2,1)
                    if(~ouprime_flag)
                        plot(r(1:max_loop_size),squeeze(data.phi_dp(1:max_loop_size,gat1)),r_avg,phidp_avg-phidp_avg(1));
                    else
                        plot(r(1:max_loop_size),squeeze(data.phi_dp(azi_id,1:max_loop_size)),r_avg,phidp_avg-phidp_avg(1));
                    end
                    xlabel('r (km)','FontSize',12)
                    ylabel('\Phi_{DP}','FontSize',12)
                    title('\Phi_{DP} vs. Range','FontSize',14)
                    legend('Raw','Avg.','Location','Best')
                    grid on
                    subplot(2,2,2)
                    if(~ouprime_flag)
                        plot(r(1:max_loop_size),squeeze(data.Z_H(1:round(azi_id*1.5),gat1)));
                    else
                        plot(r(1:max_loop_size),squeeze(data.Z_H(azi_id,1:max_loop_size)));
                    end
                    xlabel('r (km)','FontSize',12)
                    ylabel('Z_{HH}','FontSize',12)
                    title('Z_{HH} vs. Range','FontSize',14)
                    grid on
                    subplot(2,2,3)
                    if(~ouprime_flag)
                        plot(r(1:max_loop_size),squeeze(data.Z_DR(1:round(azi_id*1.5),gat1)));
                    else
                        plot(r(1:max_loop_size),squeeze(data.Z_DR(azi_id,1:max_loop_size)));                    
                    end
                    xlabel('r (km)','FontSize',12)
                    ylabel('Z_{DR}','FontSize',12)
                    title('Z_{DR} vs. Range','FontSize',14)
                    grid on
                    subplot(2,2,4)
                    if(~ouprime_flag)
                        plot(r(1:max_loop_size),squeeze(data.rho_hv(1:round(azi_id*1.5),gat1)));
                    else
                        plot(r(1:max_loop_size),squeeze(data.rho_hv(azi_id,1:max_loop_size)));
                    end
                    xlabel('r (km)','FontSize',12)
                    ylabel('\rho_{HV}','FontSize',12)
                    title('\rho_{HV} vs. Range','FontSize',14)
                    grid on
                    close(1003)
                end
                    curr_fig=gcf;
                    figure(curr_fig)
                    getvalue_flag=input('Get more values? (1/0)');
                    clear phidp_avg;
            end
            options=8;
        elseif(options==9)
            is_running=false;
        end
    end