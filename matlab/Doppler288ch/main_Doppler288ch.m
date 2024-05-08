%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%Doppler288chによるイオン温度をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all

%各PCのパスを定義
run define_path.m

%------【input】---------------------------------------------------
date = 230524;%【input】実験日
date2 = 210924;%(仮)
begin_cal = 15;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 15;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
IDS288ch.line = 'Ar';%【input】ドップラー発光ライン('Ar')

%-----------------------詳細設定【input】----------------------------
cal_Ti = false;%【input】イオン温度分布を計算(true,false)
save_Ti = false;%【input】イオン温度分布データを保存(true,false)
load_Ti = true;%【input】イオン温度分布データを読み込む(true,false)

cal_pcb = false;%【input】磁場を計算(true,false)
save_pcb = false;%【input】磁場データを保存(true,false)
load_pcb = true;%【input】磁場データを読み込む(true,false)

plot_Ti = true;%【input】イオン温度分布をプロット(true,false)
plot_psi = true;%【input】磁気面をプロット(true,false)
overlay_plot = false;%【input】イオン温度と磁気面を重ねる(true,false)

save_fig = true;%【input】イオン温度分布pngを保存(true,false)

dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
mesh_rz = 100;%【input】磁気プローブrz方向のメッシュ数(50)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
%------------------------------------------------------------------

%-----------------------解析オプション【input】----------------------------
cal_LineInt = false;%【input】線積分イオン温度、発光強度分布を計算
plot_2D_spectra = 'off';%【input】('off','all','good','bad')2次元スペクトル分布をプロット(cal_2D = trueが必要)

hw_lambda = 50;%【input】波長切り出し半幅
hw_ch = 5;%【input】CH方向切り出し半幅
hw_fit = 12;%【input】フィッティング波長切り出し半幅(< hw_lambda)
num_r = 30;%【input】r分割数(比例して計算時間が増える)
hw_lambdaA = 40;%【input】lambdaA半幅(< hw_lambda)
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
        a039shot = exp_log(i,index.a039);%a039ショット番号
        a039tfshot = exp_log(i,index.a039_TF);%a039TFショット番号
        expval.PF1 = exp_log(i,index.PF1);%PF1電圧(kV)
        expval.PF2 = exp_log(i,index.PF2);%PF2電圧(kV)
        expval.TF = exp_log(i,index.TF);%PF2電圧(kV)
        expval.EF = exp_log(i,index.EF);%EF電流(A)
        % IDS288ch.shot = exp_log(i,index.shot);%ショット番号
        % IDS288ch.exp_w = exp_log(i,index.IDS288ch_exp_w);%IDS288ch露光時間
        % IDS288ch.trg = exp_log(i,index.IDS288ch_trg);%IDS288chトリガ時間
        IDS288ch.shot = 1;%ショット番号
        IDS288ch.exp_w = 3;%IDS288ch露光時間
        IDS288ch.trg = 470;%IDS288chトリガ時間
        time = round(IDS288ch.trg+IDS288ch.exp_w/2);%磁気面プロット時間

        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_Ti
            %イオン温度分布を計算
            [z_IDS,r_IDS,Ti_IDS,Ti_max_IDS,Ti_min_IDS,Em_IDS] = cal_Doppler288ch(date,IDS288ch,pathname,cal_LineInt,plot_2D_spectra,save_Ti);
        elseif load_Ti
            %保存済みイオン温度、フローを読み取り
            [z_IDS,r_IDS,Ti_IDS,Ti_max_IDS,Ti_min_IDS,Em_IDS] = load_Doppler288ch(date2,IDS288ch,pathname);
        else
            z_IDS = char.empty;
            r_IDS = char.empty;
            Ti_IDS = char.empty;
            Ti_max_IDS = char.empty;
            Ti_min_IDS = char.empty;
            Em_IDS = char.empty;
        end

        %磁場を計算
        if cal_pcb
            [grid2D,data2D,ok_z,ok_r] = cal_pcb200ch(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
        elseif load_pcb
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
        else
            data2D = char.empty;
        end
        if not(isempty(data2D))
            %磁気面をプロット
            if plot_psi
                plot_psi200ch_at_t(time,trange,grid2D,data2D,ok_z,ok_r);
            end
        end

        %イオン温度をプロット
        if plot_Ti
            if not(isempty(Ti_IDS))
                if plot_psi
                    plot_iontemp(z_IDS,r_IDS,Ti_IDS,date2,expval,IDS288ch,pathname,overlay_plot,save_fig)
                else
                    plot_iontemp(z_IDS,r_IDS,Ti_IDS,date2,expval,IDS288ch,pathname,false,save_fig)
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
