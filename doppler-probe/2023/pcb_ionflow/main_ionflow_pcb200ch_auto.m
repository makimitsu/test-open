%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%各PCのパスを定義
run define_path.m

%------【input】-------
date = 230524;%【input】実験日
begin_cal = 15;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 24;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm](2.5)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm](4.2)
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
%------詳細設定【input】-------
cal_flow = true;%【input】流速を計算(true,false)
save_flow = true;%【input】流速データを保存(true,false)
load_flow = false;%【input】流速データを読み込む(true,false)

cal_pcb = false;%【input】磁場を計算(true,false)
save_pcb = false;%【input】磁場データを保存(true,false)
load_pcb = true;%【input】磁場データを読み込む(true,false)

plot_fit = true;%【input】ガウスフィッティングを表示(true,false)
plot_flow = true;%【input】流速をプロット(true,false)
plot_psi = true;%【input】磁気面をプロット(true,false)
overlay_plot = true;%【input】流速と磁気面を重ねる(true,false)

save_fit = true;%【input】ガウスフィッティングpngを保存(true,false)
save_fig = true;%【input】流速pngを保存(true,false)

show_offset = true;%【input】分光offsetを表示(true,false)
factor = 0.1;%【input】イオンフロー矢印サイズ(数値:0.1など)
dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
mesh_rz = 100;%【input】磁気プローブrz方向のメッシュ数(50)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)

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
        IDSP.shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,index.a039);%a039ショット番号
        a039tfshot = exp_log(i,index.a039_TF);%a039TFショット番号
        expval.PF1 = exp_log(i,index.PF1);%PF1電圧(kv)
        expval.PF2 = exp_log(i,index.PF2);%PF2電圧(kv)
        expval.TF = exp_log(i,index.TF);%PF2電圧(kv)
        expval.EF = exp_log(i,index.EF);%EF電流
        IDSP.trg = exp_log(i,index.IDSP_trg);%IDSPトリガ時間
        IDSP.exp_w = exp_log(i,index.IDSP_exp_w);%IDSP露光時間
        IDSP.gain = exp_log(i,index.IDSP_gain);%Andor gain
        time = round(IDSP.trg+IDSP.exp_w/2);%磁気面プロット時間
        min_r = exp_log(i,index.IDSP_minR);%【input】ドップラープローブ計測点最小r座標
        min_z = exp_log(i,index.IDSP_minZ);%【input】ドップラープローブ計測点最小z座標
        %ドップラープローブ計測点配列を生成
        mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_flow
            %イオン温度、フローを計算
            [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow);
        elseif load_flow
            %保存済みイオン温度、フローを読み取り
            [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
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
        end
        %イオン温度、フローをプロット
        if plot_flow
            if not(isempty(V_i))
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,'ionflow')
                else
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionflow')
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
