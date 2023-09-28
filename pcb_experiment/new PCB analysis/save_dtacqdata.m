%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
%setenv('MDSPLUS_DIR','/usr/local/mdsplus');
%addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata=getenv('rawdata_path'); %保存先


%% 自動入力の場合
prompt = {'Date:','Shot number:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2double(cell2mat(answer(2)));
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);%計測データ数
dtacqlist=[38 39 40].*ones(n_data,1);
shotlist_a038 =T.a038(IDXlist);
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a038, shotlist_a039, shotlist_a040];
tfshotlist_a038 =T.a038_TF(IDXlist);
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a038, tfshotlist_a039, tfshotlist_a040];


%% 直接入力の場合
%{
%dtacqlist=[38 39];
dtacqlist=38;

%shotlist=[11725 2546 1025];%【input】dtacqの保存番号
shotlist=11839%[10650:10692];
%shotlist=[10647 585];

%tfshotlist = [11761 2582 1061];
%tfshotlist=10646*ones(size(shotlist));
tfshotlist=zeros(size(shotlist));
n=numel(shotlist);%計測データ数
%}

%% 保存計算パート

for i=1:n_data

    dtacq_num_col=dtacqlist(i,:);
    shot_col=shotlist(i,:);
    tfshot_col=tfshotlist(i,:);
    if shot_col == tfshot_col
        tfshot_col = [0,0,0];
    end
   
    
    for j = 1:length(dtacq_num_col)
        disp(dtacq_num_col(j));
        dtacq_num = dtacq_num_col(j);
        shot = shot_col(j);
        tfshot = tfshot_col(j);

        [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
        if tfshot>0
            [rawdata0]=getMDSdata(dtacq_num,shot,0);
            save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
        end
    end

end
