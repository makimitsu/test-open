function Bt = get_B_troidal_test(date,shot,pathname)

%directory_rogo = strcat(pathname.fourier,'rogowski/');
% pathname =
% '/Users/shohgookazaki/Documents/UTokyo/OnoLab/koala/mnt/fourier/rogowski/';
directory_rogo = strcat(pathname);

aquisition_rate = 10;
offset = 0;
Bt = rogowski(date,aquisition_rate,offset,shot,directory_rogo);


end

function B_t = rogowski(date,aquisition_rate,offset,shot,directory_rogo)

% input:
%   integer: date, date of experiment. Example:(2019 Aug. 01->190801)
%   integer: shot, shot number.
%   aquisition_rate: MHz
%   offset: starting time offset in us
% output:
%   2d array of double: data, raw rogowski signal for all channels

% folder_directory_rogo = '/Users/shinjirotakeda/mountpoint/';
% folder_directory_rogo = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/probedata/';
% folder_directory_rogo = 'C:\Users\denjo\OneDrive - The University of Tokyo\data\mag_probe\';

% channels = [2,7,8]; % channels to display
channels = [1,7,8]; % channels to display
% channels =1:8; % channels to display

t_start = 1; % us
t_end = 1000; % us

t_start=t_start-offset;t_end=t_end-offset; % set offset
time_step = 1;%0.2; % us; time step of plot; must be an integer times of time step of raw data; larger time step gives faster plotting.
upper = 0.6; % upper limit of y axis
lower = -0.4; % lower limit of y axis
perc = 0.04; % span used in smoothing
windowSize = 5; % window size in filtering
%calibration = [1, -1516.4, 1, 1, 1533.2, 80, 140, 480, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
% calibration = [1, -1, 1, 1, 1, 1, 1, 1, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
calibration = [116.6647, -1, 1, 1, 1, 1, 63.9568, 1, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
plot_columns = 3;
plot_rows = fix(length(channels)/plot_columns)+(rem(length(channels),plot_columns)~=0)+1;

[date_str,shot_str,path] = directory_generation_Rogowski(date,shot,directory_rogo);

% transfer all rgw file to txt file; do nothing if there is no rgw file in the folder
% rgw2txt(date_str,'mag_probe');
% rgw2txt(date_str);

rgw2txt_shot(date_str,shot_str,directory_rogo);

desktop_path = '/Users/shohgookazaki/Desktop';
path = strcat(desktop_path, '/', date_str, shot_str, '.txt');


if ~isfile(path)
    disp(strcat('No such file: ',path));
    return
elseif (t_start<0)
    disp('offset should be smaller than t_start!');
    return
elseif aquisition_rate * time_step < 1
    disp('time resolution should be < aquisition rate!');
    return
end
    
data = readmatrix(path);
step = aquisition_rate * time_step;
x = t_start * aquisition_rate : step : t_end * aquisition_rate;
smoothed_and_filtered = zeros(length(x),length(channels));
I_TF = data(x,1+2)*calibration(1);
timing = x/aquisition_rate==475;
m0 = 4*pi*10^(-7);
r = 0.23;
B_t = m0*I_TF(timing)*1e3*12/(2*pi()*r);

% figure('Position', [0 0 1500 1500],'visible','on');
% subplot(plot_rows,plot_columns,[1,2,3]);
% hold on
% lengend_str_all = {};
% count = 1;
% for i = channels
%     smoothed = smooth(data(x,i+2),perc,'rloess')*calibration(i);
% 
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     smoothed_and_filtered(:,count) = filter(b,a,smoothed);
% 
%     %plot(x/aquisition_rate+offset,smoothed_and_filtered(:,count),'LineWidth',4);
%     plot(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'LineWidth',2);
% 
%     legend_str = ['ch',num2str(i)];
%     lengend_str_all = [lengend_str_all,legend_str];
%     count = count + 1;
% end
% 
% legend(lengend_str_all);
% ylim([lower upper]);
% xlim([t_start+offset t_end+offset]);
% xlabel('time (us)','FontSize',18);
% ylabel('kA','FontSize',18);
% ax = gca;
% ax.FontSize = 18; 
% 
% hold off
% 
% count = 1;
% for i = channels
%     legend_str = ['ch',num2str(i)];
%     subplot(plot_rows,plot_columns,count+plot_columns);
%     hold on;
%     %plot(x/aquisition_rate+offset,smoothed_and_filtered(:,count),'k','LineWidth',1);
%     %scatter(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'r');
%     plot(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'r');
%     legend(legend_str);
%     hold off;
%     ylim([lower upper]);
%     xlim([t_start+offset t_end+offset]);
%     %return_data(:,i) = smoothed_and_filtered(:,count);
%     return_data(:,i) = data(x,i+2)*calibration(i);
%     count = count+1;
% end
% 
% saveas(figure(1),[path(1:end-3),'png']);
% %close(1);
% 
% figure;
% plot(x/aquisition_rate+offset,data(x,1+2)*calibration(1),'LineWidth',2);
% hold on
% plot(x/aquisition_rate+offset,data(x,7+2)*calibration(7),'LineWidth',2);
% legend('TF coil','PF coil');
% % ylim([lower upper]);
% xlim([0 t_end+offset]);
% xlabel('time (us)','FontSize',18);
% ylabel('kA','FontSize',18);
% ax = gca;
% ax.FontSize = 18; 

end

function [date_str,shot_str,data_dir] = directory_generation_Rogowski(date,shot,directory_rogo)
    date_str = num2str(date);
    if shot < 10
        data_dir = [directory_rogo,date_str,'/',date_str,'00',num2str(shot),'.txt'];
        shot_str = ['00',num2str(shot)];
    elseif shot < 100
        data_dir = [directory_rogo,date_str,'/',date_str,'0',num2str(shot),'.txt'];
        shot_str = ['0',num2str(shot)];
    elseif shot < 1000
        data_dir = [directory_rogo,date_str,'/',date_str,num2str(shot),'.txt'];
        shot_str = num2str(shot);
    else
        disp('More than 999 shots! You need some rest!!!');
        return
    end
end

function [] = rgw2txt(date)

% Get all rgw files in the current folder
% current_folder = strcat('/Users/shinjirotakeda/mountpoint/',date,'/');
date = num2str(date);
current_folder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/probedata/',date,'/');
files = dir(strcat(current_folder,'*.rgw'));
% Loop through each
% disp(length(files));
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~,name,~] = fileparts(files(id).name);
    rename = strcat(current_folder,name,'.txt');
    old_name = strcat(current_folder,files(id).name);
    copyfile(old_name,rename);
%     disp(id);
end

end

% function [] = rgw2txt_shot(date_str,shot_str,directory_rogo)
% 
% % current_folder = strcat('/Users/shinjirotakeda/mountpoint/',date,'/');
% current_folder = strcat(directory_rogo,date_str,'/');
% filename = strcat(current_folder,date_str,shot_str,'.rgw');
% rename = strcat(current_folder,date_str,shot_str,'.txt');
% if isfile(rename)
%     return
% end
% copyfile(filename,rename, 'f');
% 
% end

function [] = rgw2txt_shot(date_str, shot_str, directory_rogo)

    % Define the current folder based on the provided directory and date
    current_folder = strcat(directory_rogo, date_str, '/');
    filename = strcat(current_folder, date_str, shot_str, '.rgw');

    % Set the destination to your desktop
    desktop_path = '/Users/shohgookazaki/Desktop';
    rename = strcat(desktop_path, '/', date_str, shot_str, '.txt');

    % Check if the file already exists on the desktop
    if isfile(rename)
        return; % If the file exists, do nothing
    end

    % Copy the file to the desktop
    try
        copyfile(filename, rename, 'f'); % Use 'f' to force overwrite
    catch ME
        disp(['Error copying file: ', ME.message]);
    end

end