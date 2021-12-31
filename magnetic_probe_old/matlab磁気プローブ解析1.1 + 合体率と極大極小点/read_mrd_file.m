function data = read_mrd_file(dir,ch_per_module,module_count,time_count)
% input: 
%   string: directory to the mrd data file
%   int: number of channels per module
%   int: number of modules
%   int: number of points in time space (1us interval)
% output: 
%   2d array of double: data stores the digitizer data (unit is mV) 
%                       data(i,j): i iterates channel, j iterates time
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
for i = 1:module_count
    module_id(i) = fscanf(fileID,'%d',1);
    calibr = fscanf(fileID,'%d',1);
    calibration_array(1+(i-1)*ch_per_module:i*ch_per_module) = str2double(regexp(num2str(calibr),'\d','match'));
end

% read 2 characters to get through newline to the start of data
fscanf(fileID,'%c',2); 

% read data to string -> split every 4 -> str2num -> reshape again
% creating a 2d array that stores main data
% data: channel_count x time_count
temp_data = fscanf(fileID,'%c');
data = reshape(str2num(reshape(temp_data,4,[])'),time_count,channel_count)';

% ******************** CALIBRATION ************************
for i = 1:channel_count
    switch calibration_array(i)
        case 0
            data(i,:) = data(i,:)*2.442;
        case 1
            data(i,:) = data(i,:)*1.221;
        case 2
            data(i,:) = data(i,:)*0.661;
        case 3
            data(i,:) = data(i,:)*0.305;
        case 4 
            data(i,:) = data(i,:)*0.153;
        case 5
            data(i,:) = data(i,:)*0.077;
        case 6
            data(i,:) = data(i,:)*48.85;
        case 7
            data(i,:) = data(i,:)*4.885;
        case 8 
            data(i,:) = data(i,:)*9.77;
        case 9
            data(i,:) = data(i,:)*0.977;
        otherwise
            data(i,:) = data(i,:)*0;
    end
end
% ************************************************
fclose(fileID);
end