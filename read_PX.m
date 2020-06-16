
%Find other moment files

%Get global attributes
cd(dir_name)
ncd = files(row).name;
data.time = ncreadatt(ncd, '/', 'Time');
data.ht = ncreadatt(ncd, '/', 'Height');
if data_type == 2
    data.va = 15.7;
else
    data.va = ncreadatt(ncd, '/', 'Nyquist_Vel-value');
end
data.lon = ncreadatt(ncd, '/', 'Longitude');
data.lat = ncreadatt(ncd, '/', 'Latitude');
data.Station = ncreadatt(ncd, '/', 'radarName-value');
if strcmp(data.Station, 'PX-1000')
    data.Station = 'PX1000';
end

data.missing_data = ncreadatt(ncd, '/', 'MissingData');
data.prf = ncreadatt(ncd, '/', 'PRF-value');
data.pw = ncreadatt(ncd, '/', 'PulseWidth-value');
data.az = ncread(ncd, 'Azimuth');
data.el = ncread(ncd, 'Elevation');
data.gw = ncread(ncd, 'GateWidth');
data.Z_H = ncread(ncd, 'Corrected_Intensity');
ncd = files_D(row).name;
data.Z_DR = ncread(ncd, 'Differential_Reflectivity');
ncd = files_P(row).name;
data.phi_dp = ncread(ncd, 'PhiDP');
ncd = files_R(row).name;
data.rho_hv = ncread(ncd, 'RhoHV');
ncd = files_V(row).name;
data.vr_h = ncread(ncd, 'Radial_Velocity');

%Check for editted radial velocity data
str = 'good';
var = 'Radial_Velocity_ed';
ncid = netcdf.open(ncd);
try
    ID = netcdf.inqVarID(ncid, var);
catch exception
    %if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
    if strcmp(exception.identifier, 'MATLAB:imagesci:netcdf:libraryFailure')
        str = 'bad';
    end
end
disp(str)
netcdf.close(ncid);
if strcmp(str, 'good')
    data.vr_h_old = data.vr_h;
    data.vr_h = ncread(ncd, 'Radial_Velocity_ed');
end

tmp_bypass = false;
if tmp_bypass
    load('NewV-PX-20130520-202018-E2.6.mat');
    data.vr_h = V(1:2010, 1:360);
end
ncd = files_W(row).name;
data.SW = ncread(ncd, 'Width');
data.Z_H(data.Z_H == data.missing_data) = NaN;
data.Z_DR(data.Z_DR == data.missing_data) = NaN;
data.vr_h(data.vr_h == data.missing_data) = NaN;
data.SW(data.SW == data.missing_data) = NaN;
data.rho_hv(data.rho_hv == data.missing_data) = NaN;
data.phi_dp(data.phi_dp == data.missing_data) = NaN;

vr_lims = [-nanmax(abs(data.vr_h(:))) nanmax(abs(data.vr_h(:)))];
cd(base_dir)
