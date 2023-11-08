function [] = rogowski_plot(date,shot,ch,calib)

% ch_cal = 1;
% ch_rogo = 1;
discharge = 1;

if discharge == 1
    t_start = 3500;
    t_end = 6000;
    p_start = 680;
    p_end = 700;
end

if discharge == 2
    t_start = 1;
    t_end = 10000;
    p_start = 4000;
    p_end = 6000;
end

% folder_directory_rogo = '/Users/Ryo/mountpoint/';
folder_directory_rogo = '/Volumes/md0/rogowski/';
date_str = num2str(date);

if shot < 10
    data_dir = [folder_directory_rogo,date_str,'/',date_str,'00',num2str(shot)];
    shot_str = ['00',num2str(shot)];
elseif shot < 100
    data_dir = [folder_directory_rogo,date_str,'/',date_str,'0',num2str(shot)];
    shot_str = ['0',num2str(shot)];
elseif shot < 1000
    data_dir = [folder_directory_rogo,date_str,'/',date_str,num2str(shot)];
    shot_str = num2str(shot);
else
    disp('More than 999 shots! You need some rest!!!')
    return
end

txt_name = strcat(data_dir,'.txt');
rgw = importdata(txt_name);
t_ax = rgw.data(t_start:t_end,2);
% data_cal = rgw.data(t_start:t_end,ch_cal+2);
data_rogo = rgw.data(t_start:t_end,3:end);
[~,n_ch] = size(data_rogo);

% offset_cal = mean(data_cal(1:50));
% data_cal = data_cal - offset_cal;
% data_rogo = data_rogo - mean(data_rogo(1:50,:));

figure('Position', [0 0 800 800],'visible','on');
for i = ch
    plot(t_ax,data_rogo(:,i)*calib(i))
    hold on
end
% title('TF: Before calibration')
% legend('Calibrator','Target Rogowski')
legendStrings = "CH" + string(ch);
legend(legendStrings)
xlabel('us')
ylabel('V')
% xlim([400 500])
% ylim([0 40])
% ylim([-30 50])
% savefig('TF_Before')
% close(f1)

% peak_cal = mean(data_cal(p_start:p_end));
% peak_rogo = mean(data_rogo(p_start:p_end));
%
% cal_value = peak_cal*10/peak_rogo%Z³ŒW”kA/V
%
% data_cal = data_cal*10;
% data_rogo = data_rogo*cal_value;
%
% f2 = figure;
% plot(t_ax,data_cal,t_ax,data_rogo)
% title('TF: After calibration')
% legend('Calibrator','Target Rogowski')
% xlabel('us')
% ylabel('kA')
% savefig('TF_After')
% close(f2)

