function [data_dir,data_TF_dir] = directory_generation(date,TF_shot,shot)

if shot < 10
    data_dir = ['Data/Bz/',num2str(date),'/',num2str(date),'00',num2str(shot),'.mrd'];
elseif shot < 100
    data_dir = ['Data/Bz/',num2str(date),'/',num2str(date),'0',num2str(shot),'.mrd'];
elseif shot < 1000
    data_dir = ['Data/Bz/',num2str(date),'/',num2str(date),num2str(shot),'.mrd'];
else
    disp('More than 999 shots! You need some rest!!!')
    return
end

if TF_shot < 10
    data_TF_dir = ['Data/Bz/',num2str(date),'/',num2str(date),'00',num2str(TF_shot),'.mrd'];
elseif TF_shot < 100
    data_TF_dir = ['Data/Bz/',num2str(date),'/',num2str(date),'0',num2str(TF_shot),'.mrd'];
elseif TF_shot < 1000
    data_TF_dir = ['Data/Bz/',num2str(date),'/',num2str(date),num2str(TF_shot),'.mrd'];
else
    disp('More than 999 shots! You need some rest!!!')
    return
end

end