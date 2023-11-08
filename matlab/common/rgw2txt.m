function [] = rgw2txt(date,shot)

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

rgw_name = strcat(data_dir,'.rgw');
txt_name= strcat(data_dir,'.txt');

if exist(rgw_name,'file')
    if not(isfile(txt_name))
        copyfile(rgw_name,txt_name);
    end
else
    warning('Rogowski file does not exist. External Bt cannot be calculated.')
end

end
