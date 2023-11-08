%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を手入力して
%磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%各PCのパスを定義
run define_path.m

%------【input】-------
date = 230816;%【input】実験日(yymmdd)
shotlist = 2178;%991:992;%【input】dtacqの保存番号
tfshotlist = 2177;%【input】dtacqのtfonly保存番号
EFlist = 150;%150;%【input】EF電流
n_plot = 16;%【input】磁気面subplot枚数(16以下の4の倍数or3以下)
t_start = 470;%【input】磁気面subplot開始時間(us)
dt = 2;%【input】磁気面subplot時間間隔(us)
cmap = false;%【input】磁気面カラーマップ('psi','Bz','Bt','Jt','Et',false)
cbar = true;%【input】カラーバー(true,false)

cal_pcb = true;%【input】磁場データを計算(true,false)
save_pcb = false;%【input】磁場データを保存(true,false)
load_pcb = false;%【input】磁場データを読み込む(true,false)
plot_psi = true;%【input】磁気面をプロット(true,false)

%------詳細設定【input】-------
dtacq.num = 39;%【input】磁気プローブdtacq番号(ほぼ固定)
mesh_rz = 50;%【input】磁気プローブrz方向のメッシュ数(ほぼ固定)
trange = 430:590;%【input】磁気プローブ計算時間範囲(ほぼ固定)

n_data=numel(shotlist);%計測データ数
for i=1:n_data
    dtacq.shot=shotlist(i);
    dtacq.tfshot=tfshotlist;
    expval.EF=EFlist;
    if cal_pcb
        [grid2D,data2D,ok_z,ok_r] = cal_pcb200ch(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
    elseif load_pcb
        [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
    end
    if not(isempty(data2D))
        if plot_psi
            plot_psi200ch_multi(n_plot,t_start,dt,trange,cmap,cbar,grid2D,data2D,ok_z,ok_r);
        end
    end
end
