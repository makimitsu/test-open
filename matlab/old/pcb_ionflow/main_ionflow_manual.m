%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータを手入力して
%ドップラープローブによるイオン温度、フローをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath '/Users/rsomeya/Documents/lab/matlab/common';

%各PCのパスを定義
run define_path.m

%------【input】-------
date = 230817;%【input】実験日
expval.PF1 = 39;%PF1電圧(kv)
expval.PF2 = 39;%PF2電圧(kv)
expval.TF = 4;%PF2電圧(kv)
expval.EF = 150;%EF電流
ICCD.shot = 13;%【input】ショット番号
ICCD.trg = 484;%【input】ICCDトリガ時間
ICCD.exp_w = 4;%【input】ICCD露光時間
ICCD.gain = 4095;%【input】ICCD gain
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
min_r = 10;%【input】ドップラープローブ計測点最小r座標[cm]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm]
min_z = 0;%【input】ドップラープローブ計測点最小z座標[cm]
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm]
n_CH = 20;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
%------詳細設定【input】-------
cal_flow = true;%【input】流速を計算(true,false)
show_offset = true;%【input】分光offsetを表示(true,false)
plot_fit = true;%【input】ガウスフィッティングを表示(true,false)
save_fit = false;%【input】ガウスフィッティングpngを保存(true,false)
save_flow = false;%【input】流速データを保存(true,false)
load_flow = false;%【input】流速データを読み込む(true,false)
plot_flow = true;%【input】流速をプロット(true,false)
save_fig = false;%【input】流速pngを保存(true,false)
factor = 0.1;%【input】イオンフロー矢印サイズ(数値:0.1など)

%計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%イオン温度、フローを計算、プロット
if cal_flow
    %イオン温度、フローを計算
    [V_i,err_V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow);
elseif load_flow
    %保存済みイオン温度、フローを読み取り
    [V_i,err_V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
end

if plot_flow
    if isempty(V_i)
        return
    else
        plot_ionflow(V_i,err_V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionflow')
    end
end