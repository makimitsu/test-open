%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン速度分布関数、温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all

%各PCのパスを定義
run define_path.m

%------【input】---------------------------------------------------
date = 230524;%【input】実験日
begin_cal = 15;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 24;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm](2.5)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm](4.2)
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
n_CH = 20;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
%-----------------------詳細設定【input】----------------------------
cal_vdist = true;%【input】速度分布を計算(true,false)
save_vdist = true;%【input】速度分布データを保存(true,false)
load_vdist = false;%【input】速度分布データを読み込む(true,false)

cal_pcb = false;%【input】磁場を計算(true,false)
save_pcb = false;%【input】磁場データを保存(true,false)
load_pcb = true;%【input】磁場データを読み込む(true,false)

plot_spectra = false;%【input】スペクトルをプロット(true,false)
plot_analisis = false;%【input】逆変換解析をプロット(true,false)
plot_vdist = false;%【input】速度分布をプロット(true,false)
plot_compare = false;%【input】再構成比較をプロット(true,false)
plot_flow = true;%【input】流速をプロット(true,false)
plot_psi = true;%【input】磁気面をプロット(true,false)
overlay_plot = true;%【input】流速と磁気面を重ねる(true,false)
plot_Br = false;

save_fig = false;%【input】速度分布、フローpngを保存(true,false)

plot_type = 'contour';%【input】速度分布プロット種類('contour','surf')
Ti_type = 'dispersion';%【input】イオン温度計算法('dispersion')

show_offset = false;%【input】分光offsetを表示(true,false)
inversion_method = 5;%【input】速度分布逆変換手法(1~6)
factor = 0.1;%【input】イオンフロー矢印サイズ(数値:0.1など)
dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
mesh_rz = 100;%【input】磁気プローブrz方向のメッシュ数(50)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
%------------------------------------------------------------------

%実験ログ読み取り
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------磁気面&フローを計算------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        error('end_cal must <= begin_cal.')
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_calからend_calまで計算
    else
        error('end_cal must <= %d.', exp_log(end_row,4))
    end
    for i = start_i:end_i
        ICCD.shot = exp_log(i,index.shot);%ショット番号
        a039shot = exp_log(i,index.a039);%a039ショット番号
        a039tfshot = exp_log(i,index.a039_TF);%a039TFショット番号
        expval.PF1 = exp_log(i,index.PF1);%PF1電圧(kv)
        expval.PF2 = exp_log(i,index.PF2);%PF2電圧(kv)
        expval.TF = exp_log(i,index.TF);%PF2電圧(kv)
        expval.EF = exp_log(i,index.EF);%EF電流
        ICCD.trg = exp_log(i,index.ICCD_trg);%ICCDトリガ時間
        ICCD.exp_w = exp_log(i,index.ICCD_exp_w);%ICCD露光時間
        ICCD.gain = exp_log(i,index.ICCD_gain);%Andor gain
        time = round(ICCD.trg+ICCD.exp_w/2);%磁気面プロット時間

        min_r = exp_log(i,index.minR);%【input】ドップラープローブ計測点最小r座標
        min_z = exp_log(i,index.minZ);%【input】ドップラープローブ計測点最小z座標
        %ドップラープローブ計測点配列を生成
        mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_vdist
            %イオン速度分布を計算
            [V_i,absV,T_i] = cal_ionvdist(date,expval,ICCD,mpoints,pathname,show_offset,plot_spectra, ...
                inversion_method,plot_analisis,plot_vdist,plot_type,save_fig,plot_compare,save_vdist,Ti_type);
        elseif load_vdist
            %保存済みイオン温度、フローを読み取り
            [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,ICCD,pathname);
            if plot_vdist
                plot_ionvdist(Vx,Vy,F,date,expval,ICCD,pathname,mpoints,ppoints,plot_type,save_fig)
            end
            if plot_compare
                plot_inversion_compare(F,W,P,Lambda,mpoints,Angle,ICCD)
            end
        else
            V_i = char.empty;
            absV = char.empty;
            T_i = char.empty;
        end
        %磁場を計算
        if cal_pcb
            [grid2D,data2D,ok_z,ok_r] = cal_pcb200ch(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
        elseif load_pcb
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
        end
        if not(isempty(data2D))
            %磁気面をプロット
            if plot_psi
                plot_psi200ch_at_t(time,trange,grid2D,data2D,ok_z,ok_r);
            end
            if plot_Br
                plot_Br200ch(cut_z_Br,n_plot,t_start,dt,trange,grid2D,data2D);
            end
        end
        %イオン温度、フローをプロット
        if plot_flow
            if not(isempty(V_i))
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,'ionvdist')
                else
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionvdist')
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
