function [B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] = get_B_z(date,TF_shot,shot,offset_TF,offset_EF)
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

%************** Folder directory *****************
%folder_directory = 'C:/Users/take_/OneDrive/デスクトップ/program/fluctuation/data/';
folder_directory_Bz = 'E:/experiment/results/ts-3u/';

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
[data_dir, data_TF_dir] = directory_generation(date,TF_shot,shot);
data = read_mrd_file(data_dir,ch_per_module,module_count,time_count);
data_return = data;
r_z = get_r_z(date);
ch_dist_file = get_ch_dist(date);

date = date;
shot_num = shot;

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
B_z_return = zeros(ch_r_count,ch_z_count,time_count);
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

% assign NaN to dead channels
data_return = data;
data = corrections(data,date);

% rearrange channels according to their pickup coil locations
for i = 1:length(r_probe)
    for j = 1:length(z_probe)
        B_z(i,j,:) = data(ch(i,j),:);
        B_z_return(i,j,:) = data_return(ch(i,j),:);
        ch_dist(i,j) = ch(i,j);
    end
end

% ********** Offset, Smoothing and Others **************
% set offset to 0 (signal should start at zero offset)
% smooth data

for i = 1:ch_r_count
    for j = 1:ch_z_count
        B_z(i,j,:) = B_z(i,j,:) - mean(B_z(i,j,10:40));
        B_z_return(i,j,:) = B_z_return(i,j,:) - mean(B_z_return(i,j,10:40));
        B_z(i,j,:) = rot90(smoothdata(B_z(i,j,:),'rloess',15));
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

    function [data_dir,data_TF_dir] = directory_generation(date,TF_shot,shot)
    if shot < 10
        data_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),'00',num2str(shot),'.mrd'];
    elseif shot < 100
        data_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),'0',num2str(shot),'.mrd'];
    elseif shot < 1000
        data_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),num2str(shot),'.mrd'];
    else
        disp('More than 999 shots! You need some rest!!!')
        return
    end


    
    if TF_shot < 10
        data_TF_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),'00',num2str(TF_shot),'.mrd'];
    elseif TF_shot < 100
        data_TF_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),'0',num2str(TF_shot),'.mrd'];
    elseif TF_shot < 1000
        data_TF_dir = [folder_directory_Bz,num2str(date),'/',num2str(date),num2str(TF_shot),'.mrd'];
    else
        disp('More than 999 shots! You need some rest!!!')
        return
    end
    end
    function data = read_mrd_file(dir,ch_per_module,module_count,time_count)
    % input: 
    %   string: directory to the mrd data file
    %   int: number of channels per module
    %   int: number of modules
    %   int: number of points in time space (1us interval)
    % output: 
    %   2d array of double: data stores the digitizer data (unit is mV) 
    %                       data(i,j): i iterates channel, j iterates time.
    %                       data is calibrated with range coefficient

    fileID = fopen(dir,'r');
    assert(module_count == fscanf(fileID,'%d',1));
    channel_count = ch_per_module * module_count;   % total number of channels

    % ******************** READ FILE **************************
    % create an array to store each module id
    module_id = zeros(module_count,1);
    % create an array to store digitizer calibration coefficients
    calibration_array = zeros(channel_count,1);

    % read module id and calibration coeff for all modules
    for i_module = 1:module_count
        module_id(i_module) = fscanf(fileID,'%d',1);
        calibr = fscanf(fileID,'%d',1);
        calibration_array(1+(i_module-1)*ch_per_module:i_module*ch_per_module) = str2double(regexp(num2str(calibr),'\d','match'));
    end

    % read 2 characters to get through newline to the start of data
    fscanf(fileID,'%c',2); 

    % read data to string -> split every 4 -> str2num -> reshape again
    % creating a 2d array that stores main data
    % data: channel_count x time_count
    temp_data = fscanf(fileID,'%c');
    data = reshape(str2num(reshape(temp_data,4,[])'),time_count,channel_count)';

    % ******************** CALIBRATION ************************
    for i_ch = 1:channel_count
        switch calibration_array(i_ch)
            case 0
                data(i_ch,:) = data(i_ch,:)*2.442;
            case 1
                data(i_ch,:) = data(i_ch,:)*1.221;
            case 2
                data(i_ch,:) = data(i_ch,:)*0.661;
            case 3
                data(i_ch,:) = data(i_ch,:)*0.305;
            case 4 
                data(i_ch,:) = data(i_ch,:)*0.153;
            case 5
                data(i_ch,:) = data(i_ch,:)*0.077;
            case 6
                data(i_ch,:) = data(i_ch,:)*48.85;
            case 7
                data(i_ch,:) = data(i_ch,:)*4.885;
            case 8 
                data(i_ch,:) = data(i_ch,:)*9.77;
            case 9
                data(i_ch,:) = data(i_ch,:)*0.977;
            otherwise
                data(i_ch,:) = data(i_ch,:)*0;
        end
    end

    % ************************************************
    fclose(fileID);
    end
    function cell = get_ch_dist(date)
    directory = 'channel distributions.xlsx';
    if date >= 191209 && date < 191218
        cell = readcell(directory,'Sheet','Bz_191209');
    elseif date >= 191218 && date < 200121
        cell = readcell(directory,'Sheet','Bz_191218');
    elseif date >= 200121
        cell = readcell(directory,'Sheet','Bz_200121');
    else
        disp('Wrong date!!!!');
        return;
    end

    end
    function r_z = get_r_z(date)
    if date < 191218 && date >= 191209
        r_z = readcell('coil_calibration.xlsx','Sheet','Bz_191209');
    elseif date >= 191218 && date < 200121
        r_z = readcell('coil_calibration.xlsx','Sheet','Bz_191218');
    elseif date >= 200121 && date < 200513
        r_z = readcell('coil_calibration.xlsx','Sheet','Bz_200121');
    elseif date >= 200513
        r_z = readcell('coil_calibration.xlsx','Sheet','Bz_200513');
    end
    end
    function data = corrections(data,date)
    % This function assigns values to dead channels
    % input: 
    %   raw data
    %   date of experiment(yymmdd)
    % output:
    %   data corrected

    if date == 191209
    data(39,:) = (data(38,:) + data(40,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(120,:) = (data(119,:) + data(121,:))/2;
    data(185,:) = (data(176,:) + data(186,:))/2;
    data(205,:) = (data(204,:) + data(206,:))/2;
    data(239,:) = 2*data(110,:) - data(109,:);
    data(240,:) = 2*data(239,:) - data(110,:);

    data(225:236,:) = data(99:110,:);
    data(237:238,:) = data(239:240,:);

    elseif date >= 191210 && date < 191218 
    data(39,:) = (data(38,:) + data(40,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(120,:) = (data(119,:) + data(121,:))/2;
    data(185,:) = (data(176,:) + data(186,:))/2;
    data(205,:) = (data(204,:) + data(206,:))/2;
    data(239,:) = 2*data(110,:) - data(109,:);
    data(240,:) = 2*data(239,:) - data(110,:);

    data(185:188,:) = data(62:65,:);
    data(225:236,:) = data(99:110,:);
    data(237:238,:) = data(239:240,:);

    elseif date >= 191218 && date < 200121
    data(39,:) = (data(38,:) + data(40,:))/2;
    data(42,:) = (data(28,:) + data(98,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(120,:) = (data(119,:) + data(121,:))/2;
    data(205,:) = (data(204,:) + data(206,:))/2;
    data(232,:) = (data(84,:) + data(98,:))/2;
    data(234,:) = (data(84,:) + data(56,:))/2;
    data(236,:) = 2*data(56,:) - data(234,:);
    data(70,:) = 2*data(236,:) - data(56,:);
    data(239,:) = 2*data(110,:) - data(109,:);

    elseif date >= 200121 && date < 200128
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);
    data(129,:) = (data(130,:) + data(128,:))/2;

    elseif date >= 200128 && date < 200121
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);

    elseif date >= 200130 && date < 200414
    %data(77,:) = (data(76,:) + data(78,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(236,:) = (data(235,:) + data(237,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);
    data(173,:) = (data(174,:) + data(172,:))/2;
    data(266,:) = (data(265,:) + data(267,:))/2;
    data(131,:) = (data(130,:) + data(132,:))/2;
    data(112,:) = (data(111,:) + data(113,:))/2;

    elseif date >= 200414 && date < 200424
    data(58,:) = (data(57,:) + data(59,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);
    data(259,:) = (data(242,:) + data(240,:))/2;
    data(173,:) = 2*data(172,:) - data(171,:);
    data(174,:) = (data(173,:) + data(175,:))/2;
    data(266,:) = (data(265,:) + data(267,:))/2;

    elseif date >= 200424 && date < 201016
    data(32,:) = (data(31,:) + data(33,:))/2;
    data(55,:) = (data(54,:) + data(56,:))/2;
    data(58,:) = (data(57,:) + data(59,:))/2;
    data(73,:) = (data(72,:) + data(74,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(112,:) = (data(111,:) + data(113,:))/2;
    data(121,:) = (data(122,:) + data(120,:))/2;
    data(123,:) = (data(124,:) + data(122,:))/2;
    data(131,:) = (data(132,:) + data(130,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);
    data(259,:) = (data(242,:) + data(240,:))/2;
    data(173,:) = 2*data(172,:) - data(171,:);
    data(174,:) = (data(173,:) + data(175,:))/2;
    
    data(266,:) = (data(265,:) + data(267,:))/2;
    data(270,:) = (data(269,:) + data(271,:))/2;
    
    %Temporary
    data(4,:) = (data(3,:) + data(5,:))/2;
    data(8,:) = (data(7,:) + data(9,:))/2;
    data(28,:) = 2*data(27,:) - data(26,:);
    data(36,:) = (data(35,:) + data(37,:))/2;
    data(38,:) = (data(177,:) + data(37,:))/2;
    data(64,:) = (data(63,:) + data(65,:))/2;
    data(97,:) = (data(96,:) + data(98,:))/2;
    data(222,:) = (data(221,:) + data(223,:))/2;
    data(227,:) = (data(228,:) + data(226,:))/2;

    elseif date >= 201016
    data(55,:) = (data(54,:) + data(56,:))/2;
    data(58,:) = (data(57,:) + data(59,:))/2;
    data(93,:) = (data(92,:) + data(94,:))/2;
    data(97,:) = (data(96,:) + data(98,:))/2;
    data(112,:) = (data(111,:) + data(113,:))/2;
    data(131,:) = (data(132,:) + data(130,:))/2;
    data(173,:) = 2*data(172,:) - data(171,:);
    data(174,:) = (data(173,:) + data(175,:))/2;
    data(219,:) = (data(110,:) + data(244,:))/2;
    data(243,:) = 2*data(242,:) - data(259,:);
    data(259,:) = (data(242,:) + data(240,:))/2;
    data(266,:) = (data(265,:) + data(267,:))/2;

    else
        disp('Wrong date!!!');
        return
    end
    end

clearvars -except B_z r_probe z_probe ch_dist B_z_return data_return shot_num date Bz_EF Br_EF;
end