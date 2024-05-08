function return_data = rogowski(date,aquisition_rate,offset,shot)
% input:
%  integer: date, date of experiment. Example:(2019 Aug. 01->190801)
%  integer: shot, shot number.
%  aquisition_rate: MHz
%  offset: starting time offset in us
% output:
%  2d array of double: data, raw rogowski signal for all channels
% folder_directory_rogo = '/Users/shinjirotakeda/mountpoint/';
folder_directory_rogo = '/Users/itsuki/Documents/東大/研究資料/Doppler_probe_data/';
% folder_directory_rogo = 'C:\Users\denjo\OneDrive - The University of Tokyo\data\mag_probe\';
channels = [2,2,2]; % channels to display
t_start = 1; % us
t_end = 700; % us
t_start=t_start-offset;t_end=t_end-offset; % set offset
time_step = 1;%0.2; % us; time step of plot; must be an integer times of time step of raw data; larger time step gives faster plotting.
upper = 0.4; % upper limit of y axis
lower = -0.6; % lower limit of y axis
perc = 0.04; % span used in smoothing
windowSize = 5; % window size in filtering
%calibration = [1, -1516.4, 1, 1, 1533.2, 80, 140, 480, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
% calibration = [116.6647, -1, 1, 1, 1, 1, 1, 1, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
% calibration = [1, 12468.6, 1, 1, 63.9568, 1, 223.2319, 1, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
calibration = [1, 12468.6, 1, 1, 1, 1, 1, 1, 1, 1,]; % calibration for each channel. I'm not sure about the exact calibration coefficient.
plot_columns = 3;
plot_rows = fix(length(channels)/plot_columns)+(rem(length(channels),plot_columns)~=0)+1;
[date_str,shot_str,path] = directory_generation_Rogowski(date,shot);
% transfer all rgw file to txt file; do nothing if there is no rgw file in the folder
% rgw2txt(date_str,'mag_probe');
rgw2txt_shot(date_str,shot_str);
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
figure('Position', [0 0 1500 1500],'visible','on');
subplot(plot_rows,plot_columns,[1,2,3]);
hold on
lengend_str_all = {};
count = 1;
for i = channels
  smoothed = smooth(data(x,i+2),perc,'rloess')*calibration(i);
  b = (1/windowSize)*ones(1,windowSize);
  a = 1;
  smoothed_and_filtered(:,count) = filter(b,a,smoothed);
  %plot(x/aquisition_rate+offset,smoothed_and_filtered(:,count),'LineWidth',4);
  plot(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'LineWidth',2);
  legend_str = ['ch',num2str(i)];
  lengend_str_all = [lengend_str_all,legend_str];
  count = count + 1;
end
legend(lengend_str_all);
% ylim([lower upper]);
xlim([t_start+offset t_end+offset]);
xlabel('time (us)','FontSize',18);
ylabel('kA','FontSize',18);
ax = gca;
ax.FontSize = 18;
hold off
count = 1;
for i = channels
  legend_str = ['ch',num2str(i)];
  subplot(plot_rows,plot_columns,count+plot_columns);
  hold on;
  %plot(x/aquisition_rate+offset,smoothed_and_filtered(:,count),'k','LineWidth',1);
  %scatter(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'r');
  plot(x/aquisition_rate+offset,data(x,i+2)*calibration(i),'r');
  legend(legend_str);
  hold off;
  ylim([lower upper]);
  xlim([t_start+offset t_end+offset]);
  %return_data(:,i) = smoothed_and_filtered(:,count);
  return_data(:,i) = data(x,i+2)*calibration(i);
  count = count+1;
end
saveas(figure(1),[path(1:end-3),'png']);
%close(1);
figure;
plot(x/aquisition_rate+offset,data(x,5+2)*calibration(5),'LineWidth',2);
hold on
plot(x/aquisition_rate+offset,data(x,7+2)*calibration(7),'LineWidth',2);
legend('TF coil','PF coil');
% ylim([lower upper]);
xlim([0 t_end+offset]);
xlabel('time (us)','FontSize',18);
ylabel('kA','FontSize',18);
ax = gca;
ax.FontSize = 18;
Iave = ( max(data(x(410:end),2+2)*calibration(2)) + max(data(x(410:end),2+2)*calibration(2)) ) / 2;
disp(Iave)
function [date_str,shot_str,data_dir] = directory_generation_Rogowski(date,shot)
  date_str = num2str(date);
  if shot < 10
    data_dir = [folder_directory_rogo,date_str,'/',date_str,'00',num2str(shot),'.txt'];
    shot_str = ['00',num2str(shot)];
  elseif shot < 100
    data_dir = [folder_directory_rogo,date_str,'/',date_str,'0',num2str(shot),'.txt'];
    shot_str = ['0',num2str(shot)];
  elseif shot < 1000
    data_dir = [folder_directory_rogo,date_str,'/',date_str,num2str(shot),'.txt'];
    shot_str = num2str(shot);
  else
    disp('More than 999 shots! You need some rest!!!');
    return
  end
end
end