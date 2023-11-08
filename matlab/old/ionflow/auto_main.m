function [] = auto_main(gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor)
%ガス種(Ar:1,H:2)/ファイバーCH数/z方向データ数(数値)/フィッティングを表示(TorF)/流速を計算(TorF)/流速を保存(TorF)/流速を表示(TorF)/figを保存(TorF)/矢印サイズ(数値:0.05など)

%入力変数
date = 230314;%実験日
begin_flow = 2;%フロー計算始めshot番号
end_flow = 2;%フロー計算終わりshot番号(0にするとbegin_flow以降全shot計算)
min_r = 12.5;%プローブ計測点最小r座標[mm]
int_r = 2.5;%プローブ計測点r方向間隔[mm]
min_z = -2.1;%プローブ計測点最小z座標[mm]
int_z = 4.2;%プローブ計測点z方向間隔[mm]

%計測点配列を生成
[r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z);

%時短用(固定値で良い)
load_s = 5300;%実験ログ読み始め行番号(230309~)
load_f = 10000;%実験ログ読み終わり行番号(固定)

%実験ログがcdになければ最新の実験ログをwebから取得
if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) > date
        disp(append(string(FileDate), '更新の実験ログを使用します。'))
    else
        run save_log.m
        disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
    end
else
    run save_log.m
    disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
end

%実験ログ中の実験日に対応する範囲を特定
exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':AR' num2str(load_f)]);
[n_row,n_col] = size(exp_log);
begin_row = find(exp_log(:,3) == date);%実験日の最初のshotの行番号を取得
if isempty(begin_row)
    disp('実験日が実験ログ中に存在しません。')
    return
end
end_row = begin_row;
while end_row<n_row && isnan(exp_log(end_row+1,3)) && exp_log(end_row+1,4)%日付がNaN&&shot番号が記入済=実験日のshot
    end_row = end_row+1;
end%実験日の最後のshotの行番号を取得

%実験ログからICCDファイル名を特定し、フローを計算
if end_flow == 0
    %begin_flow以降全部計算
    for i = begin_row + begin_flow - 1:end_row
        shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,8);%a039ショット番号
        a039tfshot = exp_log(i,9);%a039TFショット番号
        trg = exp_log(i,42);%ICCDトリガ時間
        exp_w = exp_log(i,43);%ICCD露光時間
        gain = exp_log(i,44);%Andor gain
        manual_main(true,date,shot,trg,exp_w,gain,gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor,r_measured,z_measured)
    end
else
    %begin_flowからend_flowまで計算
    for i = begin_row + begin_flow - 1:begin_row + end_flow - 1
        shot = exp_log(i,4);%ショット番号
        trg = exp_log(i,42);%ICCDトリガ時間
        exp_w = exp_log(i,43);%ICCD露光時間
        gain = exp_log(i,44);%Andor gain
        manual_main(true,date,shot,trg,exp_w,gain,gas,NofCH,nz,plot_fitting,cal_flow,save_flow,plot_flow,save_fig,factor,r_measured,z_measured)
    end
end
