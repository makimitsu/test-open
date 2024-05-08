%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%グループ化した散布図をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%各PCのパスを定義
run define_path.m

%------【input】---------------------------------------------------
date = 230315;%【input】実験日%230315
begin_cal = 4;%【input】磁気面&フロー計算始めshot番号(実験ログD列)%4
end_cal = 5;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)%35
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm](2.5)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm](4.2)
IDSP.line = "Ar";%【input】ドップラー発光ライン("Ar")
n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
t_range = [477 483];%【input】計測時刻範囲(us)
r_range = [17.5 22.5];%【input】計測点範囲(cm)
ng_shotlist = [6 9 10 26 33 34];%【input】磁場失敗ショット番号%[6 10 26 33 34]
x_type = "|Br|";%【input】x軸の変数("Vz","Vr","|Vz|","|Vr|","Bz","Br","|Bz|","|Br|","z","r","time")
y1_type = "|Vz|";%【input】y1軸の変数("Vz","Vr","|Vz|","|Vr|","Bz","Br","|Bz|","|Br|","z","r","time")
y2_type = "none";%【input】y2軸の変数("Vz","Vr","|Vz|","|Vr|","Bz","Br","z","r","time","none")
g_type = "time";%【input】グループ化基準("time","r")

%実験ログ読み取り
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

% figure("Position",[600 150 600 600])
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        error("end_cal must <= begin_cal.")
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_calからend_calまで計算
    else
        error("end_cal must <= %d.", exp_log(end_row,4))
    end

    %配列を作成
    n_scatter = end_i - start_i;%散布図プロット数
    sort = zeros(n_scatter,1);%範囲限定用
    V_z = zeros(n_scatter,1);%Vz
    V_r = zeros(n_scatter,1);%Vr
    T = zeros(n_scatter,1);%T
    B_z = zeros(n_scatter,1);%Bz
    B_r = zeros(n_scatter,1);%Br
    tp = zeros(n_scatter,1);%計測時刻
    zp = zeros(n_scatter,1);%計測点z座標
    rp = zeros(n_scatter,1);%計測点r座標
    for i = start_i:end_i
        step = i - start_i + 1;
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
        if time >= t_range(1) && time <= t_range(2) && isempty(intersect(IDSP.shot,ng_shotlist))
            %保存済みイオン温度、フローを読みこむ
            [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,IDSP,pathname);
            %磁場を読みこむ
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
            for k = 1:mpoints.n_r
                if mpoints.r(k,1) >= r_range(1) && mpoints.r(k,1) <= r_range(2)
                    IDX_z = knnsearch(grid2D.zq(1,:).',mpoints.z(k,1)*1e-2);%grid2D.zqの中で最も計測z座標に近いセル番号を取得
                    IDX_r = knnsearch(grid2D.rq(:,1),mpoints.r(k,1)*1e-2);%grid2D.rqの中で最も計測r座標に近いセル番号を取得
                    IDX_t = (time - trange(1) + 1);
                    sort(mpoints.n_r*(step - 1) + k,1) = 1;
                    V_z(mpoints.n_r*(step - 1) + k,1) = V_i(k,1);%Vz
                    V_r(mpoints.n_r*(step - 1) + k,1) = V_i(k,2);%Vr
                    T(mpoints.n_r*(step - 1) + k,1) = T_i(k,1);%イオン温度
                    B_z(mpoints.n_r*(step - 1) + k,1) = data2D.Bz(IDX_r,IDX_z,IDX_t)*1e3;%Bz[mT]
                    B_r(mpoints.n_r*(step - 1) + k,1) = data2D.Br(IDX_r,IDX_z,IDX_t)*1e3;%Br[mT]
                    tp(mpoints.n_r*(step - 1) + k,1) = trange(IDX_t);%計測時刻[us]
                    zp(mpoints.n_r*(step - 1) + k,1) = mpoints.z(k,1);%z座標[cm]
                    rp(mpoints.n_r*(step - 1) + k,1) = mpoints.r(k,1);%r座標[cm]
                end
            end
        end
    end
    indices = find(sort==0);
    V_z(indices) = [];
    V_r(indices) = [];
    T(indices) = [];
    B_z(indices) = [];
    B_r(indices) = [];
    tp(indices) = [];
    zp(indices) = [];
    rp(indices) = [];
    [n_data,~] = size(V_z);

    %グループ分けされた散布図を作成
    if y2_type == "none"
        ax_type = [x_type,y1_type];
    else
        ax_type = [x_type,y1_type,y2_type];
    end
    [~,n_ax] = size(ax_type);
    ax_val = zeros(n_data,n_ax);
    ax_name = strings(size(ax_type));
    for i = 1:n_ax
        switch ax_type(i)
            case "Vz"
                ax_val(:,i) = V_z;
                ax_name(i) = "V_z [km/s]";
            case "Vr"
                ax_val(:,i) = V_r;
                ax_name(i) = "V_r [km/s]";
            case "|Vz|"
                ax_val(:,i) = abs(V_z);
                ax_name(i) = "|V_z| [km/s]";
            case "|Vr|"
                ax_val(:,i) = abs(V_r);
                ax_name(i) = "|V_r| [km/s]";
            case "T"
                ax_val(:,i) = T;
                ax_name(i) = "T [eV]";
            case "Bz"
                ax_val(:,i) = B_z;
                ax_name(i) = "B_z [mT]";
            case "Br"
                ax_val(:,i) = B_r;
                ax_name(i) = "B_r [mT]";
            case "|Bz|"
                ax_val(:,i) = abs(B_z);
                ax_name(i) = "|B_z| [mT]";
            case "|Br|"
                ax_val(:,i) = abs(B_r);
                ax_name(i) = "|B_r| [mT]";
            case "time"
                ax_val(:,i) = tp;
                ax_name(i) = "time [us]";
            case "z"
                ax_val(:,i) = zp;
                ax_name(i) = "z [cm]";
            case "r"
                ax_val(:,i) = rp;
                ax_name(i) = "r [cm]";
            otherwise
                warning("Input error in axis type.")%ax_typeの入力エラー
                return;
        end
    end
    switch g_type
        case "time"
            g_val = tp;
            g_name = "time [us]";
        case "z"
            g_val = zp;
            g_name = "z [cm]";
        case "r"
            g_val = rp;
            g_name = "r [cm]";
        otherwise
            warning("Input error in group type.")%g_typeの入力エラー
            return;
    end
    n_group = length(unique(g_val, "rows"));
    % clr = hsv(n_group);
    clr = [0 90/255 1; 3/255 175/255 122/255; 1 165/255 0];
    gs1 = gscatter(ax_val(:,1),ax_val(:,2),g_val,clr,"o","filled");
    hold on
    if y2_type ~= "none"
        gs2 = gscatter(ax_val(:,1),ax_val(:,3),g_val,clr,"^");
        y_name = sprintf("%s  %s",ax_name(2),ax_name(3));
        title_name = sprintf("(%s, %s) VS %s grouped by %s",ax_name(2),ax_name(3),ax_name(1),g_name);
    else
        y_name = ax_name(2);
        title_name = sprintf("%s VS %s grouped by %s",ax_name(2),ax_name(1),g_name);
    end
    % title(title_name)
    % yline(0)%直線y = 0を挿入
    xlabel(ax_name(1))
    ylabel(y_name)
    pbaspect([1 2 1])
    lgd = legend("Location","northeastoutside");
    title(lgd,g_name)
    legend
    ax = gca;
    ax.FontSize = 16;
    hold off
else
    error("begin_cal must <= %d.", exp_log(end_row,4))
end
