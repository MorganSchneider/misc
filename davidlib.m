function [varargout] = davidlib(fn,varargin)
if ~exist('fn','var')
	eval('help davidlib')
    if isempty(whos), feval('demo'); return; end
else
    if nargout == 0; feval(fn, varargin{:});
    elseif nargout == 1, x = feval(fn,varargin{:}); varargout(1) = {x};
    end
end
return

%Lists all files
function [] = flist(files)
    for idx = 1:max(size(files))
        fprintf([num2str(idx) ': ' files(idx).name '\n'])
    end
return

%Calculates the height of the radar beam
function [bh] = bh_calc(r_km, ele)
    a = 6.371e3; ae = 4/3*a;
    bh = (r_km .^ 2 + ae .^ 2 + 2 * r_km * ae .* sind(ele)) .^ (1/2) - ae;
return

%Calculates size of the resolution volume
function [res_size] = res_size_calc(r_km, bw_deg, dr_m)
    bw = bw_deg * pi / 180; r = r_km * 1000;
    res_size = (bw * r) ^ 2 * dr_m;
return

%This function calculates kdp using a linear fit. The filt_length is the
%number of range gates over which the filter is applied. An additional
%median filter can be applied to further smooth the data. A rhohv threshold
%is also applied to avoid calculating kdp in areas of non-meteorological
%scatterers (since the primary use here is attenuation correction).
function [kdp_data] = kdp_calc(phidp, r, rhohv, filt_flag, filt_length, min_pts)
    %Filt flag: apply median filter (1: Yes, 0: No)
    %filt_length=60; rho_thres=0; min_pts=round(filt_length*0.9);
    rho_thres = 0;
    std_phi_thres = 5;
    phidp(rhohv < rho_thres) = NaN;
    if filt_flag
        b = [1.625807356e-2 2.230852545e-2 2.896372364e-2 3.595993808e-2 4.298744446e-2 4.971005447e-2...
        5.578764970e-2 6.089991897e-2 6.476934523e-2 6.718151185e-2 6.800100000e-2]; 
        b = [b, fliplr(b(1:max(size(b))-1))];
        b = b / sum(b);
        a = 1;
    end
    dr = r(2) - r(1);
    kdp = nan(size(phidp)); %phif = kdp; 
    phif = phidp;
    
    if max(size(r)) == size(phidp,1)
        for idx = 1:size(phidp,2)
%             phif = phidp(:,idx);
            %Apply smoothing (filter cosine and sine and recombine)
            %Remove nonmeteorological scatterers and compute backscatter
            %differential phase
%             phi_cos = cosd(phi); phi_sin = sind(phi);
%             yn = isfinite(phi_cos) & isfinite(phi_sin);
%             %apply a polyfit over different segments
%             if sum(yn) > 0
%                 phi_cos_filt = interp1(r(yn),phi_cos(yn),r,'linear');
% %             phi_cos_filt = filter(b,a,phi_cos);
% %             phi_sin_filt = filter(b,a,phi_sin);
%                 phi_sin_filt = interp1(r(yn),phi_sin(yn),r,'linear');
%                 phif(:,idx) = atan2(phi_sin_filt,phi_cos_filt)*180/pi;
%             else
%                 phif(:,idx) = NaN;
%             end
            if filt_flag
                tmp = phif(:,idx);
                filt_phi = filter(b,1,tmp);
                num_it = 1;
                while std(abs(tmp-filt_phi), 'omitnan') > std_phi_thres
                    tmp = filt_phi;
                    filt_phi = filter(b,1,tmp);
                    num_it = num_it+1;
                    if num_it > 100
                        keyboard;
                    end
                end
                phif(:,idx) = filt_phi;
            end
            
            for jdx = ceil(filt_length/2): size(phidp,1)-ceil(filt_length/2)
                tmp = jdx - ceil(filt_length/2) + 1 : jdx + ceil(filt_length/2) - 1;
                xtmp = squeeze(r(tmp)); ytmp = squeeze(phif(tmp,idx))';
                if sum(~isnan(ytmp)) > min_pts
                    xtmp = xtmp(isfinite(ytmp)); ytmp = ytmp(isfinite(ytmp));
                    if max(size(xtmp)) ~= max(size(ytmp))
                        keyboard;
                    end
                    coeff = polyfit(xtmp,ytmp,1);
                    kdp(jdx,idx) = coeff(1)/2;
%                     figure(135)
%                     plot(xtmp,ytmp,'-k')
%                     title(['K_{dp}: ' num2str(kdp(jdx,idx)) '\n'])
%                     pause;
                else
                    kdp(jdx,idx) = NaN;
                end
            end
        end
        
        if filt_flag
            npts = 5; conv_type = 1;
            for idx = 1:size(phidp,1)
                [data_filt] = median_filt(kdp(idx,:),npts,conv_type);
                kdpf(idx,:) = data_filt;
            end
            kdp_data.kdpf = kdpf;
        end
        kdp_data.kdp = kdp;
        kdp_data.phif = phif;
%         kdp(2:size(phif,2)-1,:)=(phif(3:size(phif,2),:)-phif(1:size(phif,2)-2,:))./(2*dr);
    else
        keyboard;
        for idx = 1:size(phidp,1)
            phi = phidp(idx,:);
            phi_cos = cosd(phi); phi_sin = sind(phi);
            yn = isfinite(phi_cos) & isfinite(phi_sin);
            if sum(yn) > 0
                phi_cos_filt = interp1(r(yn),phi_cos(yn),r,'linear');
    %             phi_cos_filt = filter(b,a,phi_cos);
    %             phi_sin_filt = filter(b,a,phi_sin);
                phi_sin_filt = interp1(r(yn),phi_sin(yn),r,'linear');
                keyboard;
                phif(idx,:) = atan2(phi_sin_filt,phi_cos_filt)*180/pi;
            else
                phif(idx,:) = NaN;
            end
            for jdx = ceil(filt_length/2) : size(phidp,2) - ceil(filt_length/2)
                tmp = jdx - ceil(filt_length/2) + 1 : jdx + ceil(filt_length/2) - 1;
                xtmp = r(tmp); ytmp = phif(idx,tmp);
                if sum(~isnan(ytmp)) > round(filt_length/5)
                    xtmp = xtmp(isfinite(ytmp)); ytmp = ytmp(isfinite(ytmp));
                    coeff = polyfit(xtmp,ytmp,1);
                    kdp(jdx,idx) = coeff(1)/2;
                    keyboard;
                else
                    kdp(jdx,idx) = NaN;
                end
            end
            %Apply smoothing
        end
        
        if filt_flag
            npts = 5; conv_type = 1;
            for idx = 1:size(phidp,2)
                [data_filt] = median_filt(kdp(:,idx),npts,conv_type);
                kdpf(:,idx) = data_filt;
            end
            kdp_data.kdpf = kdpf;
        end
        kdp_data.kdp = kdp;
        kdp_data.phif = phif;
%         kdp(:,2:size(phif,1)-1)=(phif(:,3:size(phif,1))-phif(:,1:size(phif,1)-2))./(2*dr);
    end
return

%This function applies attenuation and differential attenuation correction
%to data. The cumulative differential phase is calculated by summing KDP
%over each range gate. Then, the cumulative differential phase is
%multiplied by the linear coefficients alpha and beta to compute the
%attenuation correction, which is added to the input Z and ZDR to provide
%the attenuation corrected data in the structure data_corr.
function [data_corr] = att_corr(Z, ZDR, kdp, r, alpha, beta, gate_start)
    kdp(~isfinite(kdp)) = 0; phi_sum = zeros(size(ZDR));
    if max(size(r)) == size(ZDR,1)
        kdp(1:gate_start-1,:) = 0;
        for idx = gate_start:size(ZDR,1)
            phi_sum(idx,:) = nanmean(kdp(gate_start:idx,:),1) * 2 * r(idx);
        end
    else
        kdp(:,1:gate_start-1) = 0;
        for idx = gate_start:size(ZDR,2)
            phi_sum(:,idx) = nanmean(kdp(:,gate_start:idx),2) * 2 * r(idx);
        end
    end
    phi_sum(phi_sum < 0) = 0;
    data_corr.Z = Z + phi_sum * alpha;
    data_corr.ZDR = ZDR + phi_sum * beta;
    data_corr.phi_sum = phi_sum;
return

%1D N-point median filter 
function [data_filt] = median_filt(data, npts, conv_type)
%conv_type=1: circular convolution, otherwise not implemented
    data_filt = nan(size(data));
    if conv_type == 1
        dat = [data(max(size(data))-floor(npts/2)) : data(max(size(npts))), data, data(1:npts+floor(npts/2))];
        for idx = 1:max(size(data))
            tmp = dat(idx:idx+npts-1);
            if min(size(median(tmp(isfinite(tmp))))) > 0
                data_filt(idx) = median(tmp(isfinite(tmp)));
            else
                data_filt(idx) = NaN;
            end
        end
    else
        keyboard;
    end
return

%Calculates x and y from r_km and az
function [xx,yy] = radpol2cart_ppi(r,az)
    if size(r,1) == 1
        rr = repmat(r, [max(size(az)) 1]);
    else
        rr = repmat(r', [max(size(az)) 1]);
    end
    if size(az,2) == 1
        azz = repmat(az, [1 max(size(r))]);
    else
        azz = repmat(az', [1 max(size(r))]);
    end
    xx = rr .* sind(azz);
    yy = rr .* cosd(azz);
return

%Calculates x and z from r_km and az
function [coords] = radpol2cart_rhi(r,el)
    if size(r,1) == 1
        rr = repmat(r, [max(size(el)) 1]);
    else
        rr = repmat(r', [max(size(el)) 1]);
    end
    if size(el,2) == 1
        ell = repmat(el, [1 max(size(r))]);
    else
        ell = repmat(el', [1 max(size(r))]);
    end
    coords.xx = rr .* cosd(ell);
    coords.zz = bh_calc(rr,ell);
return

% The engine behind all the color maps
function [cmap] = fleximap(num,pt)
if nargin < 1
    num = 15;
    pt = [0    0.5 0.0 0.0;...
          0.20 1.0 0.0 0.0;...
          0.50 1.0 1.0 1.0;...
          0.80 0.0 0.0 1.0;...
          1.00 0.0 0.0 0.5];
end
if any(pt > 1), fprintf('PT Matrix cannot have value > 1.0.\n'); return; end
x = linspace(0,1,num); x = x(:);
cmap = interp1(pt(:,1), pt(:,2:4), x, 'linear');
cmap(cmap < 0) = 0;
cmap(cmap > 1) = 1;
return

function [x] = bjetmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.00 0.00 0.50;...
      0.13 0.00 0.10 1.00;...
      0.37 0.00 1.00 1.00;...
      0.50 1.00 1.00 1.00;...
      0.57 1.00 1.00 0.00;...
      0.88 1.00 0.00 0.00;...
      1.00 0.50 0.00 0.00];
x = feval('fleximap',num,pt);
return