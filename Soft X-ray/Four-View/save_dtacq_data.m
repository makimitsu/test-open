function save_dtacq_data(dtacq_num,shot,tfshot,rawdataPath)

% rawdataPath=getenv('rawdata_path'); %保存先

% %%%%実験オペレーションの取得
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
% node='date';
% % pat=230707;
% pat=date;
% T=searchlog(T,node,pat);
% IDXlist=6;%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
% IDXlist = shot;
% n_data=numel(IDXlist);%計測データ数
% shotlist=T.a039(IDXlist);
% tfshotlist=T.a039_TF(IDXlist);
% EFlist=T.EF_A_(IDXlist);
% TFlist=T.TF_kV_(IDXlist);
% dtacqlist=39.*ones(n_data,1);

% dtacqlist=39;
% shotlist= 1819;%【input】dtacqの保存番号
% tfshotlist = 1817;

% dtacqlist=40;
% shotlist= 299;%【input】dtacqの保存番号
% tfshotlist = 297;

% date = 230706;%【input】計測日
% n=numel(shotlist);%計測データ数

%RC係数読み込み

% for i=1:n
%     dtacq_num=dtacqlist(i);
%     shot=shotlist(i);
%     tfshot=tfshotlist(i);
%     [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
%     save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
%     if tfshot>0
%         [rawdata]=getMDSdata(dtacq_num,shot,0);
%         % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
%         save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
%     end
% end

% dtacq_num=dtacqlist(i);
% shot=shotlist(i);
% tfshot=tfshotlist(i);

[rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
if tfshot==0
    [rawdata]=getMDSdata(dtacq_num,shot,0);
    % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
    save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
end

end