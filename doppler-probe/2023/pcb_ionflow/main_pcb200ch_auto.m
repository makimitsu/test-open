%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を実験ログから
%自動取得して磁気面などをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%各PCのパスを定義
run define_path.m

%------【input】-------
date = 230526;%【input】実験日(yymmdd)
begin_cal = 7;%【input】計算始めshot番号(実験ログD列)
end_cal = 12;%【input】計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降同日の全shot計算)
n_plot = 1;%【input】磁気面subplot枚数(16以下の4の倍数or3以下)
t_start = 475;%【input】磁気面subplot開始時間(us)
dt = 2;%【input】磁気面subplot時間間隔(us)
cmap = false;%【input】磁気面カラーマップ('psi','Bz','Br','Bt','Jt','Et',false)
cbar = true;%【input】カラーバー(true,false)
%------詳細設定【input】-------
cal_pcb = true;%【input】磁場を計算(true,false)

save_pcb = true;%【input】磁場データを保存(true,false)
load_pcb = false;%【input】磁場データを読み込む(true,false)

tfoffset = true;%【input】tf真空磁場を含めるかどうか(true,false)

plot_psi = false;%【input】磁気面をプロット(true,false)

plot_Br = false;%【input】Brをプロット(true,false)
cut_z_Br = 2.1;%【input】Brプロットのz断面[cm]

plot_magpres = false;%【input】磁気圧をプロット(true,false)
cut_z_magpres = 1;%【input】磁気圧プロットのz断面[cm]

dtacq.num = 39;%【input】磁気プローブdtacq番号(ほぼ固定)
mesh_rz = 100;%【input】磁気プローブrz方向のメッシュ数(ほぼ固定)
trange = 430:590;%【input】磁気プローブ計算時間範囲(ほぼ固定)

%実験ログ読み取り
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------磁気面を計算------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        warning('end_cal must <= begin_cal.')
        return
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;
    else
        warning('end_cal must <= %d.', exp_log(end_row,4))
        return
    end
    for i = start_i:end_i
        a039shot = exp_log(i,index.a039);%a039ショット番号
        a039tfshot = exp_log(i,index.a039_TF);%a039TFショット番号
        expval.PF1 = exp_log(i,index.PF1);%PF1電圧(kv)
        expval.PF2 = exp_log(i,index.PF2);%PF2電圧(kv)
        expval.TF = exp_log(i,index.TF);%PF2電圧(kv)
        expval.EF = exp_log(i,index.EF);%EF電流
        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_pcb
            if tfoffset == false
                [grid2D,data2D,ok_z,ok_r] = cal_pcb200ch(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
            else
                [grid2D,data2D_tfoffset,ok_z,ok_r] = cal_pcb200ch_with_tfoffset(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
                data2D = char.empty;
            end
        elseif load_pcb
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
            [~,data2D_tfoffset,~,~] = load_pcb200ch_with_tfoffset(date,dtacq,pathname);
        end
        if not(isempty(data2D))
            if plot_psi
                plot_psi200ch_multi(n_plot,t_start,dt,trange,cmap,cbar,grid2D,data2D,ok_z,ok_r);
            end
            if plot_Br
                plot_Br200ch(cut_z_Br,n_plot,t_start,dt,trange,grid2D,data2D);
            end
            if plot_magpres
                plot_magpres200ch(cut_z_magpres,n_plot,t_start,dt,trange,grid2D,data2D,data2D_tfoffset);
            end
        end
    end
else
    warning('begin_cal must <= %d.', exp_log(end_row,4))
    return
end
