%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata='/Users/yunhancai/Google Drive/Data/pcb_probe/raw'; %保存先

a039_shot = 2024;
a039_tfshot = 2032 * ones(size(a039_shot));
a040_shot = a039_shot - 1521;
a040_tfshot = a039_tfshot - 1521;

dtacqlist=39;
shotlist=a039_shot;%【input】dtacqの保存番号
tfshotlist = a039_tfshot;
n=numel(shotlist);%計測データ数

for i=1:n
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
    if tfshot < 0
        tfshot = 0;
    end
    save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    disp(['completed ',num2str(i),'/',num2str(n)]);
    if tfshot>0
        [rawdata0]=getMDSdata(dtacq_num,shot,0);
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
    end
end

dtacqlist=40;
shotlist=a040_shot;%【input】dtacqの保存番号
tfshotlist = a040_tfshot;
n=numel(shotlist);%計測データ数

for i=1:n
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
    if tfshot < 0
        tfshot = 0;
    end
    save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    disp(['completed ',num2str(i),'/',num2str(n)]);
    if tfshot>0
        [rawdata0]=getMDSdata(dtacq_num,shot,0);
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
    end
end

