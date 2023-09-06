%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
%setenv('MDSPLUS_DIR','/usr/local/mdsplus');
%addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata=getenv('rawdata_path'); %保存先

%dtacqlist=[38 39 40];
%dtacqlist=[38 39];
dtacqlist=38;

%shotlist=[11725 2546 1025];%【input】dtacqの保存番号
shotlist=11797%[10650:10692];
%shotlist=[10647 585];

%tfshotlist = [11761 2582 1061];
%tfshotlist=10646*ones(size(shotlist));
tfshotlist=zeros(size(shotlist));
n=numel(shotlist);%計測データ数


for i=1:n
    dtacq_num=dtacqlist;
    %dtacq_num=dtacqlist(i);
    shot=shotlist(i)
    tfshot=tfshotlist(i);
    [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
    save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    if tfshot>0
        [rawdata0]=getMDSdata(dtacq_num,shot,0);
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
    end
end
