%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
pathname.rawdata=getenv('rawdata_path'); %保存先

%%%%実験オペレーションの取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
% --ここいじる
pat=230929;
% --
T=searchlog(T,node,pat);
IDXlist = T.shot;% すべてのshotについて保存している，選択したければ配列で渡せばOK．
n_data=numel(IDXlist);% 計測データ数
shotlist_38=T.a038(IDXlist);% DTACQ番号リスト
shotlist_39=T.a039(IDXlist);% DTACQ番号リスト
shotlist_40=T.a040(IDXlist);
tfshotlist_38=T.a038_TF(IDXlist);
tfshotlist_39=T.a039_TF(IDXlist);
tfshotlist_40=T.a040_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);

%RC係数読み込み

for i=1:n_data
    disp(strcat(num2str(i), '/', num2str(n_data)));
    % dtacq38
    shot=shotlist_38(i);
    tfshot=tfshotlist_38(i);
    [rawdata]=getMDSdata(38,shot,tfshot);% 測定した生信号
    save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(38),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    % tfshotが0でない，つまり通常のshotである場合は，TFを差し引かない信号も保存する．
    if tfshot>0
        [rawdata0]=getMDSdata(38,shot,0);
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(38),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
    end
%     % dtacq39
%     shot=shotlist_39(i);
%     tfshot=tfshotlist_39(i);
%     [rawdata]=getMDSdata(39,shot,tfshot);% 測定した生信号
%     save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
%     % tfshotが0でない，つまり通常のshotである場合は，TFを差し引かない信号も保存する．
% %     if tfshot>0
% %         [rawdata0]=getMDSdata(39,shot,0);
% %         save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
% %     end
%     % dtacq40
%     shot=shotlist_40(i);
%     tfshot=tfshotlist_40(i);
%     [rawdata]=getMDSdata(40,shot,tfshot);%測定した生信号
%     save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
%     % tfshotが0でない，つまり通常のshotである場合は，TFを差し引かない信号も保存する．
%     %     if tfshot>0
% %         [rawdata0]=getMDSdata(40,shot,0);
% %         save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
% %     end
end