close all

%各PCのパスを定義
run define_path.m

%------【input】-------
date = 230524;%【input】実験日230524,230526
shotlist_cal = 15:24;%【input】磁気面&フロー計算ショットリスト15:24,[2:5,7:12]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm](2.5)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm](4.2)
IDSP.line = 'Ar';%【input】ドップラー発光ライン('Ar')
n_CH = 20;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
start_r = 1;%【input】プロット始めCH
end_r = 5;%【input】プロット終わりCH
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
dtacq.num = 39;

mu0 = 4*pi*1e-7;%真空の透磁率
[~,n_shot] = size(shotlist_cal);%計算ショット数

%実験ログ読み取り
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%全ショットのデータを入れる配列を準備
V_i_all = zeros(n_CH/4,2*n_z,n_shot);
absV_all = zeros(n_CH/4,n_z,n_shot);
T_i_all = zeros(n_CH/4,n_z,n_shot);
magpres_t_all = zeros(100,1,n_shot);
magpres_tp_all = zeros(100,1,n_shot);
magpres_z_all = zeros(100,1,n_shot);
magpres_all = zeros(100,1,n_shot);
diff_magpres_all = zeros(100-1,1,n_shot);

for i = 1:n_shot
    i_log = shotlist_cal(1,i) + begin_row - 1;
    IDSP.shot = exp_log(i_log,index.shot);%ショット番号
    a039shot = exp_log(i_log,index.a039);%a039ショット番号
    a039tfshot = exp_log(i_log,index.a039_TF);%a039TFショット番号
    IDSP.trg = exp_log(i_log,index.IDSP_trg);%IDSPトリガ時間
    IDSP.exp_w = exp_log(i_log,index.IDSP_exp_w);%IDSP露光時間
    IDSP.gain = exp_log(i_log,index.IDSP_gain);%Andor gain
    time = round(IDSP.trg+IDSP.exp_w/2);%計測時刻
    min_r = exp_log(i_log,index.minR);%ドップラープローブ計測点最小r座標
    min_z = exp_log(i_log,index.minZ);%ドップラープローブ計測点最小z座標
    cut_z_magpres = min_z;%磁気圧計算z座標
    if dtacq.num == 39
        dtacq.shot = a039shot;
        dtacq.tfshot = a039tfshot;
    end
    %ドップラープローブ計測点配列を生成
    mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);
    %保存済みイオン温度、フローを読み取り
    % [V_i_all(:,:,i),absV_all(:,:,i),T_i_all(:,:,i)] = load_ionflow(date,IDSP,pathname);
    [V_i_all(:,:,i),absV_all(:,:,i),T_i_all(:,:,i),~,~,~,~,~,~,~,~] = load_ionvdist(date,IDSP,pathname);
    %保存済み磁場データを読み取り
    [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
    [~,data2D_tfoffset,~,~] = load_pcb200ch_with_tfoffset(date,dtacq,pathname);
    %磁気圧計算
    cut_z_magpres = cut_z_magpres*1e-2;%単位を[cm]から[m]に変換
    i_pcb = time-trange(1)+1; %磁場データの時刻index
    IDX = knnsearch(grid2D.zq(1,:).',cut_z_magpres);%grid2D.zqの中で最もcut_zに近いセル番号を取得
    z = grid2D.zq(1,IDX)*1e2;%z座標[cm]
    [nR,~,~] = size(data2D.Bz);
    magpres_z = zeros(nR,1);
    magpres_t = zeros(nR,1);
    magpres_tp = zeros(nR,1);
    magpres = zeros(nR,1);
    for j = 1:nR
        magpres_z(j,1) = data2D.Bz(j,IDX,i_pcb)^2/(2*mu0);
        magpres_t(j,1) = data2D_tfoffset.Bt(j,IDX,i_pcb)^2/(2*mu0);
        magpres_tp(j,1) = data2D.Bt(j,IDX,i_pcb)^2/(2*mu0);
        magpres(j,1) = magpres_z(j,1) + magpres_t(j,1);
    end
    magpres_t_all(:,1,i) = magpres_t;
    magpres_tp_all(:,1,i) = magpres_tp;
    magpres_z_all(:,1,i) = magpres_z;
    magpres_all(:,1,i) = magpres;

    for k = 1:99
        dr = grid2D.rq(2,1) - grid2D.rq(1,1);
        diff_magpres_all(k,1,i) = (magpres_all(k+1,1,i) - magpres_all(k,1,i))/(grid2D.rq(k+1,1) - grid2D.rq(k,1));
    end
end

if not(exist([pathname.mat,'/ionflow_magpres/',num2str(date)],'dir'))
    mkdir(sprintf("%s/ionflow_magpres", pathname.mat), sprintf("%s", num2str(date)));
end
save([pathname.mat,'/ionflow_magpres/',num2str(date),'/shot',num2str(shotlist_cal(1)),'-',num2str(shotlist_cal(n_shot)),'.mat'],'grid2D','mpoints','V_i_all','T_i_all','magpres_t_all','magpres_tp_all','magpres_z_all','magpres_all','diff_magpres_all')

% V_i_mean = mean(V_i_all,3);
% V_i_sigma = std(V_i_all,0,3);
% magpres_t_mean = mean(magpres_t_all,3);
% magpres_t_sigma = std(magpres_t_all,0,3);
% magpres_z_mean = mean(magpres_z_all,3);
% magpres_z_sigma = std(magpres_z_all,0,3);
% magpres_mean = mean(magpres_all,3);
% magpres_sigma = std(magpres_all,0,3);
% figure('Position',[600 800 600 600])
% %左縦軸
% yyaxis left
% errorbar(mpoints.r(start_r:end_r),-V_i_mean(start_r:end_r,2),V_i_sigma(start_r:end_r,2),'LineWidth',2)%Vr
% xlabel('R [cm]')
% ylabel('Ion Velocity [km/s]')
% left_y_upper = max(-V_i_mean(start_r:end_r,2) + V_i_sigma(start_r:end_r,2));
% left_y_lower = min(-V_i_mean(start_r:end_r,2) - V_i_sigma(start_r:end_r,2));
% left_ylim = round(max([left_y_upper -left_y_lower]))+1;
% ylim([-left_ylim left_ylim])
% yline(0,'--','LineWidth',3);
% %右縦軸
% yyaxis right
% unit = 1E2;
% errorbar(grid2D.rq(:,1)*unit,magpres_mean(:,1),magpres_sigma(:,1),'LineWidth',2)%磁気圧
% ylabel('Magneic Pressure [Pa]')
% right_y_upper = max(magpres_mean(:,1) + magpres_sigma(:,1));
% right_y_lower = min(magpres_mean(:,1) - magpres_sigma(:,1));
% right_ylim = round(max([right_y_upper -right_y_lower]))+1;
% ylim([-right_ylim right_ylim])
% xlim([8 20])
% % title([num2str(time),'us'],'Color','black','FontWeight','bold')
% % legend('V_z','V_r')
% ax = gca;
% ax.FontSize = 16;
