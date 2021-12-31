function [B_z,r_probe,z_probe,ch_dist] = get_B_z(data_dir,data_TF_dir)
% input:
%   string: directory of the interested raw data mrd file
%   string: directory of TF discharge only raw data mrd file
% output:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   2d array of int: ch_dist, keep track of channel locations

%************** Directories **********************
pickup_coil_coeffs_dir = 'pickupcoil_coefficient.xlsx';

%************** Physical Constants ***************
mu0   = 4*pi*10^(-7);   % vacuum permeability (H/m)

%************** EF Related ***********************
z_EF = 0.702;   % EF z location (m)
r_EF = 0.5;     % EF radius (m)
n_EF = 234;     % EF turns
i_EF = 140;     % EF current (A)

%************** B_z Probe Properties**************
ch_per_module = 4;                              % 4 channels per module
module_count = 70;                              % total number of modules
time_count = 1024;                              % 1024 points in time space
max_ch = 184;                                   % number of channels deployed
broken_ch = [57,64,72,78,81,120,134,138,...
            153,157,166,169,171];               % channels broken
unused_ch = [92,93,94,95,96,97,98,99,100];      % channels unused
                                                % probe locations in r and z
r_probe = [ 0.2750, 0.2500, 0.2200, 0.1900, 0.1670, 0.1210, 0.0750];
z_probe = [ 0.2085, 0.1835, 0.1585, 0.1335, 0.1085,...
            0.0835, 0.0585, 0.0485, 0.0385, 0.0285,...
            0.0185, 0.0085, 0.0000,-0.0080,-0.0180,...
           -0.0280,-0.0380,-0.0480,-0.0580,-0.0830,...
           -0.1080,-0.1330,-0.1580,-0.1830,-0.2080];
ch_r_count = 7;                                 % 7 ch along r
ch_zL_count = 12;                               % 12 ch along z on left
ch_zR_count = 13;                               % 13 ch along z on right
ch_z_count = ch_zL_count+ch_zR_count;           % total ch along z

%************** Read Files ***********************
% read: 
%   data (channel_count x time_count)
%   data_TF (channel_count x time_count)
%   pickup_coil_coeff 
%   pickup_coil_direction (+-)
%   pickup coil location (r,z)
data = read_mrd_file(data_dir,ch_per_module,module_count,time_count);
data_TF = read_mrd_file(data_TF_dir,ch_per_module,module_count,time_count);
pickup_coil_coeff = xlsread(pickup_coil_coeffs_dir,'B2:B281');
pickup_coil_direction = xlsread(pickup_coil_coeffs_dir,'C2:C281');

% Get rid of TF noise
data = data - data_TF;

%**************** get probe_B_z ****************
% create channel_distri to keep track of channels
channel_distri = (1:184);
channel_distri(92:100) = -1; % chs not used
% replace probe insertions
channel_distri(173) = 157; 
channel_distri(177:184) = (189:196);
channel_distri(39) = 92; 
channel_distri(174) = 171;

% replace probe insertions
% 173->157,177:184->189:196,39->92,174->171
data(173,:) = data(157,:);
data(177:184,:) = data(189:196,:);
data(39,:) = data(92,:);
data(174,:) = data(171,:);

% calculate probe_B_z from data and coefficients; mT->T
probe_B_z = zeros(max_ch,time_count);
for i = 1:max_ch
    probe_B_z(i,:) = data(i,:)*pickup_coil_coeff(i)*pickup_coil_direction(i)*0.001;
end

% dead inserted probes:
% 57,64,72,78,81,120,134,138,153,157,166,169,171
probe_B_z(57,:)  =  (probe_B_z(56,:)+probe_B_z(58,:))/2;
probe_B_z(64,:)  =  (probe_B_z(63,:)+probe_B_z(65,:))/2;
probe_B_z(72,:)  =  2*probe_B_z(71,:)-probe_B_z(70,:);
probe_B_z(78,:)  =  2*probe_B_z(77,:)-probe_B_z(76,:);
probe_B_z(81,:) =  (probe_B_z(80,:)+probe_B_z(82,:))/2;
probe_B_z(120,:) =  (probe_B_z(119,:)+probe_B_z(121,:))/2;
probe_B_z(134,:) =  (probe_B_z(133,:)+probe_B_z(135,:))/2;  
probe_B_z(138,:) =  (probe_B_z(137,:)+probe_B_z(139,:))/2;
probe_B_z(153,:) =  (probe_B_z(152,:)+probe_B_z(154,:))/2;
probe_B_z(157,:) =  (probe_B_z(156,:)+probe_B_z(158,:))/2;
probe_B_z(166,:) =  2*probe_B_z(165,:)-probe_B_z(164,:);
probe_B_z(169,:) =  (probe_B_z(168,:)+probe_B_z(170,:))/2;    
probe_B_z(171,:) =  (probe_B_z(170,:)+probe_B_z(172,:))/2;
% mark dead probes
for i = broken_ch
    channel_distri(i) = 0;
end

% rearrange probe index according to inserted location
% B_z([r][z][t])
B_z = zeros(ch_r_count,ch_z_count,time_count);
B_z(1,1:13,:) = flipud(probe_B_z(1:13,:));
B_z(2,1:13,:) = flipud(probe_B_z(14:26,:));
B_z(3,1:13,:) = flipud(probe_B_z(79:91,:));
B_z(4,1:13,:) = flipud(probe_B_z(27:39,:));
B_z(5,1:13,:) = flipud(probe_B_z(40:52,:));
B_z(6,1:13,:) = flipud(probe_B_z(53:65,:));
B_z(7,1:13,:) = flipud(probe_B_z(66:78,:));
B_z(1,14:25,:) = probe_B_z(101:112,:);
B_z(2,14:25,:) = probe_B_z(113:124,:);
B_z(3,14:25,:) = probe_B_z(173:184,:);
B_z(4,14:25,:) = probe_B_z(125:136,:);
B_z(5,14:25,:) = probe_B_z(137:148,:);
B_z(6,14:25,:) = probe_B_z(149:160,:);
B_z(7,14:25,:) = probe_B_z(161:172,:);
% reshape channel tracking to 2d
ch_dist = zeros(ch_r_count,ch_z_count);
ch_dist(1,1:13) = flip(channel_distri(1:13));
ch_dist(2,1:13) = flip(channel_distri(14:26));
ch_dist(3,1:13) = flip(channel_distri(79:91));
ch_dist(4,1:13) = flip(channel_distri(27:39));
ch_dist(5,1:13) = flip(channel_distri(40:52));
ch_dist(6,1:13) = flip(channel_distri(53:65));
ch_dist(7,1:13) = flip(channel_distri(66:78));
ch_dist(1,14:25) = channel_distri(101:112);
ch_dist(2,14:25) = channel_distri(113:124);
ch_dist(3,14:25) = channel_distri(173:184);
ch_dist(4,14:25) = channel_distri(125:136);
ch_dist(5,14:25) = channel_distri(137:148);
ch_dist(6,14:25) = channel_distri(149:160);
ch_dist(7,14:25) = channel_distri(161:172);

% EF Calculation
% 2d mesh indicating probe location ([r][z])
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);
% Bz and Br due to EF current at B_z probe positions
[Bz_EF,Br_EF] = B_EF(z_EF,r_EF,i_EF,n_EF,probe_mesh_r,probe_mesh_z);

% total B_z = probe_B_z + Bz_EF
B_z = B_z + Bz_EF;

% ********** Offset, Smoothing and Others **************
% set offset to 0 (signal should start at zero)
% smooth data
for i = 1:ch_r_count
    for j = 1:ch_z_count
        B_z(i,j,:) = B_z(i,j,:) - mean(B_z(i,j,1:10));
        B_z(i,j,:) = rot90(smoothdata(B_z(i,j,:),'rloess',5));
    end
end

clearvars -except B_z r_probe z_probe ch_dist;

end