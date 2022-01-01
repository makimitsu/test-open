function [low_n_data] = low_n_mode(date,shot,TF_shot,TF_offset)

%## prepare folder directory ##
%folder_directory_Bz = 'C:/Users/take_/OneDrive/デスクトップ/program/モード計測/data/Bz/';
%folder_directory_Bz= 'D:\pub\mnt\old-koala\experiment\results\ts-3u/';
folder_directory_Bz = 'E:/experiment/results/ts-3u/';

%## make directory generation ##
[data_dir_oscilloscope,TF_dir_oscilloscope] = directory_generation_rgw(date,shot,TF_shot);

%## get struct data ##
struct_data = tdfread(data_dir_oscilloscope);
TF_struct_data = tdfread(TF_dir_oscilloscope);

%## convert struct to table ##
low_n_data = struct2table(struct_data);
TF_low_n_data = struct2table(TF_struct_data);

%## convert table to array ##
low_n_data = table2array(low_n_data);
TF_low_n_data = table2array(TF_low_n_data);


%## offset TF data ##
if TF_offset
    low_n_data = low_n_data - TF_low_n_data;
end


    function [data_dir_oscilloscope,TF_dir_oscilloscope] = directory_generation_rgw(date,shot,TF_shot)    
     if shot < 10
        data_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),'00',num2str(shot),'.rgw'];
    elseif shot < 100
        data_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),'0',num2str(shot),'.rgw'];
    elseif shot < 1000
        data_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),num2str(shot),'.rgw'];
    else
        disp('More than 999 shots! You need some rest!!!')
        return
     end
    
      if TF_shot < 10
        TF_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),'00',num2str(TF_shot),'.rgw'];
    elseif shot < 100
        TF_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),'0',num2str(TF_shot),'.rgw'];
    elseif shot < 1000
        TF_dir_oscilloscope = [folder_directory_Bz,num2str(date),'/',num2str(date),num2str(TF_shot),'.rgw'];
    else
        disp('More than 999 shots! You need some rest!!!')
        return
      end
    
    
    end
    
    end