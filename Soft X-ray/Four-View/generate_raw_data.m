function generate_raw_data(date)

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）


%%%%実験オペレーションの取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
% date=230714;
T=searchlog(T,node,date);

n_data=numel(T.date);%計測データ数
shotlist_a039 =T.a039;
shotlist_a040 = T.a040;
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF;
tfshotlist_a040 =T.a040_TF;
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
% EFlist=T.EF_A_;
% TFlist=T.TF_kV_;
% dtacqlist=39.*ones(n_data,1);

sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
% r_shift = 0.00;
% ok = logical(C(:,14));
dtacq_num_list = C(:,1);

for i = 1:n_data
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    % i_EF=EFlist(i);
    % TF=TFlist(i);
    if ismember(39,dtacq_num_list)
        filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
        if exist(filename1,"file")==0
            % disp('No rawdata file of a039 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
            % disp(['File:',filename1,' does not exit']);
            % return
        end
    end
    if ismember(40,dtacq_num_list)
        filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
        if exist(filename2,"file")==0
            % disp('No rawdata file of a040 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
            % disp(['File:',filename2,' does not exit']);
            % return
        end
    end
end