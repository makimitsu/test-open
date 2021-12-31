function [B_z,r_probe,z_probe,ch_dist,data_return] = get_B_z(date,TF_shot,shot,offset_TF,offset_EF)
% input:
%   integer: date of experiment. Example:(2019 Aug. 01->190801)
%   integer: TF shot number. Example:(2)
%   integer: shot number. Example:(38)
%   boolean: option for offseting TF signal
%   boolean: option for offseting EF signal
% output:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   2d array of int: ch_dist, keep track of channel locations

%************** Physical Constants ***************
mu0   = 4*pi*10^(-7);   % vacuum permeability (H/m)

%************** EF Related ***********************
z1_EF = 0.680;  % EF z1 location (m)
z2_EF = -0.680; % EF z2 location (m)
r_EF = 0.5;     % EF radius (m)
n_EF = 234;     % EF turns
i_EF = 150;     % EF current (A)

%************** B_z Probe Properties**************
ch_per_module = 4;                              % 4 channels per module
module_count = 69;                              % total number of modules
time_count = 1024;                              % 1024 points in time space
max_ch = 232;                                   % number of channels deployed
ch_r_count = 8;                                 % 8 ch along r
ch_zL_count = 14;                               % 14 ch along z on left(oku)
ch_zR_count = 15;                               % 15 ch along z on right(temae)
ch_z_count = ch_zL_count+ch_zR_count;           % total ch along z

%************** Read Files ***********************
% read: 
%   data (channel_count x time_count)
%   data_TF (channel_count x time_count)
%   pickup_coil_coeff 
%   pickup_coil_direction (+-)
%   pickup coil location (r,z)
[data_dir,data_TF_dir] = directory_generation(date,TF_shot,shot);
data = read_mrd_file(data_dir,ch_per_module,module_count,time_count);
if date < 191218 && date >= 191209
    r_z = readcell('coil_calibration.xlsx','Sheet','Bz_191209');
elseif date >= 191218 && date < 200121
    r_z = readcell('coil_calibration.xlsx','Sheet','Bz_191218');
elseif date >= 200121
    r_z = readcell('coil_calibration.xlsx','Sheet','Bz_200121');
end
ch_dist_file = get_ch_dist('channel distributions.xlsx',date);
data_return = data;

% Get rid of TF noise
if offset_TF
    data_TF = read_mrd_file(data_TF_dir,ch_per_module,module_count,time_count);
    data = data - data_TF;
end

%**************** get probe_B_z ****************

% get r and z coordinates of pickup coils
r_probe = rot90(cell2mat(ch_dist_file(3:end,2)));
z_probe = cell2mat(ch_dist_file(2,3:end));
ch = cell2mat(ch_dist_file(3:end,3:end));

% get information from 'coil_coordinate.xlsx'
r = cell2mat(r_z(2:end,4));
z = cell2mat(r_z(2:end,5));
dir = cell2mat(r_z(2:end,6));
coefficient = cell2mat(r_z(2:end,7));
B_z = zeros(ch_r_count,ch_z_count,time_count);
ch_dist = zeros(ch_r_count,ch_z_count);

% calculate probe_B_z from data and coefficients; mT->T
for i = 1:length(r_probe)
    for j = 1:length(z_probe)
        r_ok = logical(r == r_probe(i));
        z_ok = logical(z == z_probe(j));
        ok = r_ok.*z_ok;
        index = find(ok,1,'first');
        data(ch(i,j),:) = data(ch(i,j),:)*dir(index)*coefficient(index)*0.001;
    end
end

% assign value to dead channels
data = corrections(data,date);

% rearrange channels according to their pickup coil locations
for i = 1:length(r_probe)
    for j = 1:length(z_probe)
        B_z(i,j,:) = data(ch(i,j),:);
        ch_dist(i,j) = ch(i,j);
    end
end

% ********** Offset, Smoothing and Others **************
% set offset to 0 (signal should start at zero offset)
% smooth data

for i = 1:ch_r_count
    for j = 1:ch_z_count
        B_z(i,j,:) = B_z(i,j,:) - mean(B_z(i,j,100:300));
        B_z(i,j,:) = rot90(smoothdata(B_z(i,j,:),'rloess',5));
    end
end

% EF Calculation
% 2d mesh indicating probe location ([r][z])

if offset_EF
    [probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);
    % Bz and Br due to EF current at B_z probe positions
    [Bz_EF,Br_EF] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,probe_mesh_r,probe_mesh_z,false);
    
    % total B_z = probe_B_z + Bz_EF
    B_z = -Bz_EF+B_z;
end

clearvars -except B_z r_probe z_probe ch_dist data_return;
end