%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%------【input】-------
% % FC合体、X点R=0.2m、ExBアウトフロー小。IDSP->230828,230829(delay=480,484,488us)
% IDSP.date = 230828;%【input】IDSP実験日
% IDSPshotlist = [5 6 8:12 15:17 19:23 25:27 31:61 63];
% plot_time = 482;%【input】プロット時間[us]482,486

% SEP合体、X点R=0.26m、ExBアウトフロー大。IDSP->230830,230831(delay=468, 472, 476us)
IDSP.date = 230830;%【input】IDSP実験日
IDSPshotlist = [13 14 16 17 22:26 28 29 32 34 37 41:45 47 49 51 54 55 57 59 60];%[13 14 16 17 20 22 23 25 26 32 37 41:45 47 49 51 54 55 57];%【input】IDSPshot番号リスト
plot_time = 470;%【input】プロット時間[us]470,474

switch IDSP.date
    case 230828
        savename.ExB = [pathname.mat,'/ExB/','230828_shot5-63-a039_2301_', num2str(plot_time - 2), '_1_5.mat'];
        savename.curve = [pathname.mat,'/curve/','230828_a039_2301_', num2str(plot_time - 2), '_1_5.mat'];
        savename.nablaB = [pathname.mat,'/nablaB/','230828_a039_2301_', num2str(plot_time - 2), '_1_5.mat'];
    case 230830
        savename.ExB = [pathname.mat,'/ExB/','230830_shot11-60-a039_2437_' , num2str(plot_time - 2), '_1_5.mat'];
        savename.curve = [pathname.mat,'/curve/','230830_a039_2437_', num2str(plot_time - 2), '_1_5.mat'];
        savename.nablaB = [pathname.mat,'/nablaB/','230830_a039_2437_', num2str(plot_time - 2), '_1_5.mat'];
end

IDSP.n_CH = 28;%【input】IDSPファイバーCH数(28)
IDSP.n_r = 7;%【input】IDSPr方向データ数(7)
IDSP.n_z = 1;%【input】IDSPz方向データ数(1)
IDSP.int_r = 2.5;%【input】IDSPr方向計測点間隔[cm](2.5)
profileplot = 'Vr';%【input】プロット種類('Vr','Vz','T','offset')
cal_type = 'ionflow';%【input】IDSP計算法('ionflow','ionvdist')
ExB_mean = true;
curve_mean = true;
nablaB_mean = true;
magpresplot = 'Pz';%【input】プロット種類('Pr','Pz','Pt',Pall)
magpres_mean = true;
single_range = [1:4 6 7];%【input】プロットCH([1:4 7])
wav_range = [1:12 16:21];%【input】プロットr([1:21],[1:12 19:21])
plot_single = false;%【input】IDSPシングルショットをプロット
plot_VExB = false;%【input】ExBをプロット
plot_Vcurve = true;%【input】curve driftをプロット
plot_VnablaB = true;%【input】∇B driftをプロット
plot_WAV = false;%【input】IDSP重み付き平均をプロット
plot_magpres = false;%【input】磁気圧をプロット

%------IDSP詳細設定【input】------
plot_fit = false;%【input】ガウスフィッティングを表示(true,false)
save_fit = false;%【input】ガウスフィッティングpngを保存(true,false)
show_offset = false;%【input】分光offsetを表示(true,false)
plot_spectra = false;%【input】スペクトルをプロット(true,false)
plot_analisis = false;%【input】逆変換解析をプロット(true,false)
Ti_type = 'dispersion';%【input】イオン温度計算法('dispersion')
inversion_method = 5;%【input】速度分布逆変換手法(1~6)

%--------実験ログ取得---------
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,IDSP.date);
IDSPshotlist(ismember(IDSPshotlist,T.shot)==0) = [];
n_data=numel(IDSPshotlist);%計測shot数
PCBshotlist_a039 =T.a039(IDSPshotlist);
PCBshotlist_a040 = T.a040(IDSPshotlist);
PCBshotlist = [PCBshotlist_a039, PCBshotlist_a040];
PCBtfshotlist_a039 =T.a039_TF(IDSPshotlist);
PCBtfshotlist_a040 =T.a040_TF(IDSPshotlist);
PCBtfshotlist = [PCBtfshotlist_a039, PCBtfshotlist_a040];
PF1list=T.CB1_kV_(IDSPshotlist);
PF2list=T.CB2_kV_(IDSPshotlist);
TFlist=T.TF_kV_(IDSPshotlist);
EFlist=T.EF_A_(IDSPshotlist);
Gaslist=T.gas(IDSPshotlist);
IDSPminZlist=T.IDSPZ_cm_(IDSPshotlist);
IDSPminRlist=T.IDSPMinR_cm_(IDSPshotlist);
IDSPdelaylist=T.IDSPDelay_us_(IDSPshotlist);
IDSPwidthlist=T.IDSPWidth_us_(IDSPshotlist);
IDSPgainlist=T.IDSPGain(IDSPshotlist);
IDSPtimelist=round(IDSPdelaylist+IDSPwidthlist/2);

FIG.tate = 1;%プロット枚数(縦)
FIG.yoko = 1;%プロット枚数(横)
FIG.dt = 0;%プロット時間間隔[us]

buf_minR = sort(rmmissing(unique(IDSPminRlist)));
for i=1:size(buf_minR,1)
    buf_r = ((buf_minR(i):IDSP.int_r:buf_minR(i)+(IDSP.n_r-1)*IDSP.int_r)*1E-2)';
    if i == 1
        WAV_IDSPdata.r = buf_r;
    else
        WAV_IDSPdata.r = cat(1,WAV_IDSPdata.r,buf_r);%統合後IDSPdata R座標[m]
    end
end
WAV_IDSPdata.z = sort(rmmissing(unique(IDSPminZlist)))*1E-2;%統合後IDSPdata Z座標[m]
WAV_IDSPdata.time = sort(rmmissing(unique(IDSPtimelist)));%統合後IDSPdata 時間[us]
WAV_IDSPdata.V_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1)*2,size(WAV_IDSPdata.time,1));%V_i重み付き平均
WAV_IDSPdata.err_V_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1)*2,size(WAV_IDSPdata.time,1));%V_i重み付き平均の誤差分散
WAV_IDSPdata.T_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1),size(WAV_IDSPdata.time,1));%T_i重み付き平均
WAV_IDSPdata.err_T_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1),size(WAV_IDSPdata.time,1));%T_i重み付き平均の誤差分散
sum_weight_V_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1)*2,size(WAV_IDSPdata.time,1));
sum_weight_T_i = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1),size(WAV_IDSPdata.time,1));
cnt = zeros(size(WAV_IDSPdata.r,1),size(WAV_IDSPdata.z,1),size(WAV_IDSPdata.time,1));
legendStrings = "";
hs = char.empty;

figure('Position', [0 0 500 500],'visible','on');
if plot_magpres
    yyaxis left
end
for i=1:n_data
    IDSP.shot=IDSPshotlist(i);
    IDSP.line=Gaslist(i);
    IDSP.delay=IDSPdelaylist(i);
    IDSP.width=IDSPwidthlist(i);
    IDSP.gain=IDSPgainlist(i);
    IDSP.time=IDSPtimelist(i);
    IDSP.z = IDSPminZlist(i)*ones(7,1)*1E-2;%IDSP計測点Z座標[m]
    IDSP.r = ((IDSPminRlist(i):IDSP.int_r:IDSPminRlist(i)+(IDSP.n_r-1)*IDSP.int_r)*1E-2)';%IDSP計測点R座標[m]
    EXP.PF1=PF1list(i);
    EXP.PF2=PF2list(i);
    EXP.TF=TFlist(i);
    EXP.EF=EFlist(i);
    FIG.start = IDSP.time;%プロット開始時刻[us]
    if not(isnan(IDSP.delay)) && (IDSP.time == plot_time)
        %IDSP計算
        switch cal_type
            case 'ionflow'
                IDSPdata = cal_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit);
            case 'ionvdist'
                IDSPdata = cal_ionvdist(IDSP,pathname,show_offset,plot_spectra,inversion_method,plot_analisis,Ti_type);
        end
        if plot_single
            if IDSP.time == plot_time
                switch IDSP.time
                    case {470,482}
                        switch cal_type
                            case 'ionflow'
                                switch profileplot
                                    case 'Vz'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),IDSPdata.err_V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),IDSPdata.err_V_i(single_range,2),'k*-')
                                    case 'T'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),IDSPdata.err_T_i(single_range,1),'k*-')
                                    case 'offset'
                                        plot(single_range,IDSPdata.offset(single_range,1),'*-')
                                        % IDSPdata.offset(6,1) - IDSPdata.offset(3,1)
                                end
                            case 'ionvdist'
                                switch profileplot
                                    case 'Vz'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),'k*-')
                                    case 'T'
                                        plot(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),'k*-')
                                end
                        end
                    case {474,486}
                        switch cal_type
                            case 'ionflow'
                                switch profileplot
                                    case 'Vz'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),IDSPdata.err_V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),IDSPdata.err_V_i(single_range,2),'k*-')
                                    case 'T'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),IDSPdata.err_T_i(single_range,1),'k*-')
                                end
                            case 'ionvdist'
                                switch profileplot
                                    case 'Vz'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),'k*-')
                                    case 'T'
                                        plot(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),'k*-')
                                end
                        end
                    case {478,490}
                        switch cal_type
                            case 'ionflow'
                                switch profileplot
                                    case 'Vz'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),IDSPdata.err_V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),IDSPdata.err_V_i(single_range,2),'k*-')
                                    case 'T'
                                        errorbar(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),IDSPdata.err_T_i(single_range,1),'k*-')
                                end
                            case 'ionvdist'
                                switch profileplot
                                    case 'Vz'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,1),'k*-')
                                    case 'Vr'
                                        plot(IDSPdata.r(single_range),IDSPdata.V_i(single_range,2),'k*-')
                                    case 'T'
                                        plot(IDSPdata.r(single_range),IDSPdata.T_i(single_range,1),'k*-')
                                end
                        end
                end
            end
            hold on
        end
        switch cal_type
            case 'ionflow'
                %加重平均計算
                idx_r = find(WAV_IDSPdata.r==IDSPdata.r(1,1));
                idx_z = find(WAV_IDSPdata.z==IDSPdata.z(1,1));
                IDSPdata.time = round(IDSPdata.delay+IDSPdata.width/2);
                idx_time = find(WAV_IDSPdata.time==IDSPdata.time);
                weight_V_i = 1./IDSPdata.err_V_i.^2;%重み1/σ^2
                WAV_IDSPdata.V_i(idx_r:idx_r+IDSP.n_r-1,2*(idx_z-1)+1:2*(idx_z-1)+2,idx_time) = WAV_IDSPdata.V_i(idx_r:idx_r+IDSP.n_r-1,2*(idx_z-1)+1:2*(idx_z-1)+2,idx_time) + IDSPdata.V_i.*weight_V_i;
                sum_weight_V_i(idx_r:idx_r+IDSP.n_r-1,2*(idx_z-1)+1:2*(idx_z-1)+2,idx_time) = sum_weight_V_i(idx_r:idx_r+IDSP.n_r-1,2*(idx_z-1)+1:2*(idx_z-1)+2,idx_time) + weight_V_i;
                weight_T_i = 1./IDSPdata.err_T_i.^2;%重み1/σ^2
                WAV_IDSPdata.T_i(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) = WAV_IDSPdata.T_i(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) + IDSPdata.T_i.*weight_T_i;
                sum_weight_T_i(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) = sum_weight_T_i(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) + weight_T_i;
                cnt(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) = cnt(idx_r:idx_r+IDSP.n_r-1,idx_z,idx_time) +1;
            case 'ionvdist'
        end
    end
end
switch cal_type
    case 'ionflow'
        WAV_IDSPdata.V_i = WAV_IDSPdata.V_i./sum_weight_V_i;
        WAV_IDSPdata.err_V_i = 1./sqrt(sum_weight_V_i);
        WAV_IDSPdata.T_i = WAV_IDSPdata.T_i./sum_weight_T_i;
        WAV_IDSPdata.err_T_i = 1./sqrt(sum_weight_T_i);
    case 'ionvdist'
end

if plot_VExB
    if exist(savename.ExB,"file")
        load(savename.ExB,'ExBdata2D')
        idx_IDSP_z = knnsearch(ExBdata2D.zq(1,:)',IDSP.z(1));
        idx_IDSP_z_min = knnsearch(ExBdata2D.zq(1,:)',IDSP.z(1)-0.016);
        idx_IDSP_z_max = knnsearch(ExBdata2D.zq(1,:)',IDSP.z(1)+0.016);
        idx_ExB_time = knnsearch(ExBdata2D.time,plot_time);
        VExB_z_mean = mean(ExBdata2D.VExB_z(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        VExB_z_mean = mean(VExB_z_mean,2);
        VExB_r_mean = mean(ExBdata2D.VExB_r(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        VExB_r_mean = mean(VExB_r_mean,2);
        err_ExB_z_y = 3;
        err_ExB_r_y = 3;
        err_ExB_x = 0.01;
        switch plot_time
            case {470,482}
                if ExB_mean
                    switch profileplot
                        case 'Vz'
                            errorbar(ExBdata2D.rq(:,idx_IDSP_z),VExB_z_mean,err_ExB_z_y,err_ExB_z_y,err_ExB_x,err_ExB_x,'bo','LineWidth',3)
                        case 'Vr'
                            errorbar(ExBdata2D.rq(:,idx_IDSP_z),VExB_r_mean,err_ExB_r_y,err_ExB_r_y,err_ExB_x,err_ExB_x,'bo','LineWidth',3)
                    end
                else
                    switch profileplot
                        case 'Vz'
                            errorbar(ExBdata2D.rq(:,idx_IDSP_z),ExBdata2D.VExB_z(:,idx_IDSP_z,idx_ExB_time),err_ExB_z_y,err_ExB_z_y,err_ExB_x,err_ExB_x,'bo','LineWidth',3)
                        case 'Vr'
                            errorbar(ExBdata2D.rq(:,idx_IDSP_z),ExBdata2D.VExB_r(:,idx_IDSP_z,idx_ExB_time),err_ExB_r_y,err_ExB_r_y,err_ExB_x,err_ExB_x,'bo','LineWidth',3)
                    end
                end
        end
        hold on
        h = plot(nan,nan,'bo-','LineWidth',3);
        hold on
        if isempty(hs)
            hs = h;
        else
            hs = [hs h];
        end
        if legendStrings == ""
            legendStrings = "ExB";
        else
            legendStrings = [legendStrings "ExB"];
        end
    else
        warning([savename.ExB, 'does not exist.'])
    end
end

if plot_Vcurve
    if exist(savename.curve,"file")
        load(savename.curve,'Curvedata2D')
        idx_IDSP_z = knnsearch(Curvedata2D.zq(1,:)',IDSP.z(1));
        idx_IDSP_z_min = knnsearch(Curvedata2D.zq(1,:)',IDSP.z(1)-0.016);
        idx_IDSP_z_max = knnsearch(Curvedata2D.zq(1,:)',IDSP.z(1)+0.016);
        idx_curve_time = knnsearch(Curvedata2D.trange,plot_time);
        Vcurve_z_mean = mean(Curvedata2D.Vcurve_z(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        Vcurve_z_mean = mean(Vcurve_z_mean,2);
        Vcurve_r_mean = mean(Curvedata2D.Vcurve_r(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        Vcurve_r_mean = mean(Vcurve_r_mean,2);
        err_curve_z_y = 1E-2;
        err_curve_r_y = 1E-2;
        err_curve_x = 0.01;
        switch plot_time
            case {470,482}
                if curve_mean
                    switch profileplot
                        case 'Vz'
                            errorbar(Curvedata2D.rq(:,idx_IDSP_z),Vcurve_z_mean,err_curve_z_y,err_curve_z_y,err_curve_x,err_curve_x,'mo','LineWidth',3)
                        case 'Vr'
                            errorbar(Curvedata2D.rq(:,idx_IDSP_z),Vcurve_r_mean,err_curve_r_y,err_curve_r_y,err_curve_x,err_curve_x,'mo','LineWidth',3)
                    end
                else
                    switch profileplot
                        case 'Vz'
                            errorbar(Curvedata2D.rq(:,idx_IDSP_z),Curvedata2D.Vcurve_z(:,idx_IDSP_z,idx_curve_time),err_curve_z_y,err_curve_z_y,err_curve_x,err_curve_x,'mo','LineWidth',3)
                        case 'Vr'
                            errorbar(Curvedata2D.rq(:,idx_IDSP_z),Curvedata2D.Vcurve_r(:,idx_IDSP_z,idx_curve_time),err_curve_r_y,err_curve_r_y,err_curve_x,err_curve_x,'mo','LineWidth',3)
                    end
                end
        end
        hold on
        h = plot(nan,nan,'mo-','LineWidth',3);
        hold on
        if isempty(hs)
            hs = h;
        else
            hs = [hs h];
        end
        if legendStrings == ""
            legendStrings = "Curvature";
        else
            legendStrings = [legendStrings "Curvature"];
        end
    else
        warning([savename.curve, 'does not exist.'])
    end
end

if plot_VnablaB
    if exist(savename.nablaB,"file")
        load(savename.nablaB,'NablaBdata2D')
        idx_IDSP_z = knnsearch(NablaBdata2D.zq(1,:)',IDSP.z(1));
        idx_IDSP_z_min = knnsearch(NablaBdata2D.zq(1,:)',IDSP.z(1)-0.016);
        idx_IDSP_z_max = knnsearch(NablaBdata2D.zq(1,:)',IDSP.z(1)+0.016);
        idx_nablaB_time = knnsearch(NablaBdata2D.trange,plot_time);
        VnablaB_z_mean = mean(NablaBdata2D.VnablaB_z(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        VnablaB_z_mean = mean(VnablaB_z_mean,2);
        VnablaB_r_mean = mean(NablaBdata2D.VnablaB_r(:,idx_IDSP_z_min:idx_IDSP_z_max,:),3);
        VnablaB_r_mean = mean(VnablaB_r_mean,2);
        err_nablaB_z_y = 5E-2;
        err_nablaB_r_y = 5E-2;
        err_nablaB_x = 0.01;
        switch plot_time
            case {470,482}
                if nablaB_mean
                    switch profileplot
                        case 'Vz'
                            errorbar(NablaBdata2D.rq(:,idx_IDSP_z),VnablaB_z_mean,err_nablaB_z_y,err_nablaB_z_y,err_nablaB_x,err_nablaB_x,'co','LineWidth',3)
                        case 'Vr'
                            errorbar(NablaBdata2D.rq(:,idx_IDSP_z),VnablaB_r_mean,err_nablaB_r_y,err_nablaB_r_y,err_nablaB_x,err_nablaB_x,'co','LineWidth',3)
                    end
                else
                    switch profileplot
                        case 'Vz'
                            errorbar(NablaBdata2D.rq(:,idx_IDSP_z),NablaBdata2D.VnablaB_z(:,idx_IDSP_z,idx_nablaB_time),err_nablaB_z_y,err_nablaB_z_y,err_nablaB_x,err_nablaB_x,'co','LineWidth',3)
                        case 'Vr'
                            errorbar(NablaBdata2D.rq(:,idx_IDSP_z),NablaBdata2D.VnablaB_r(:,idx_IDSP_z,idx_nablaB_time),err_nablaB_r_y,err_nablaB_r_y,err_nablaB_x,err_nablaB_x,'co','LineWidth',3)
                    end
                end
        end
        hold on
        h = plot(nan,nan,'co-','LineWidth',3);
        hold on
        if isempty(hs)
            hs = h;
        else
            hs = [hs h];
        end
        if legendStrings == ""
            legendStrings = "GradB";
        else
            legendStrings = [legendStrings "GradB"];
        end
    else
        warning([savename.nablaB, 'does not exist.'])
    end
end

if plot_WAV
    idx_WAV_time = knnsearch(WAV_IDSPdata.time,plot_time);
    switch profileplot
        case 'Vz'
            buf = cat(2,WAV_IDSPdata.r,WAV_IDSPdata.V_i(:,1,idx_WAV_time));
            buf = sortrows(cat(2,buf,WAV_IDSPdata.err_V_i(:,1,idx_WAV_time)));
        case 'Vr'
            buf = cat(2,WAV_IDSPdata.r,WAV_IDSPdata.V_i(:,2,idx_WAV_time));
            buf = sortrows(cat(2,buf,WAV_IDSPdata.err_V_i(:,2,idx_WAV_time)));
        case 'T'
            buf = cat(2,WAV_IDSPdata.r,WAV_IDSPdata.T_i(:,1,idx_WAV_time));
            buf = sortrows(cat(2,buf,WAV_IDSPdata.err_T_i(:,1,idx_WAV_time)));
    end
    switch profileplot
        case 'Vz'
            err_WAV_y = buf(wav_range,3)*6;
        case 'Vr'
            err_WAV_y = buf(wav_range,3)*6;
        case 'T'
            err_WAV_y = buf(wav_range,3);
    end
    err_WAV_xpos = 0.01;
    err_WAV_xneg = 0.015;
    switch plot_time
        case {470,482}
            errorbar(buf(wav_range,1),buf(wav_range,2),err_WAV_y,err_WAV_y,err_WAV_xneg,err_WAV_xpos,'r+','LineWidth',3)
        case {474,486}
            errorbar(buf(wav_range,1),buf(wav_range,2),err_WAV_y,err_WAV_y,err_WAV_xneg,err_WAV_xpos,'r+','LineWidth',3)
        case {478,490}
            errorbar(buf(wav_range,1),buf(wav_range,2),err_WAV_y,err_WAV_y,err_WAV_xneg,err_WAV_xpos,'r+','LineWidth',3)
    end
    hold on
    h = plot(nan,nan,'ro-','LineWidth',3);
    hold on
    if isempty(hs)
        hs = h;
    else
        hs = [hs h];
    end
    if legendStrings == ""
        legendStrings = "Ion Flow";
    else
        legendStrings = [legendStrings "Ion Flow"];
    end
end
title([num2str(plot_time),'us'])
xlabel('R [m]')
switch profileplot
    case 'Vz'
        ylabel('V_Z [km/s]')
        yline(0,'--k','LineWidth',2)
    case 'Vr'
        ylabel('V_R [km/s]')
        yline(0,'--k','LineWidth',2)
    case 'T'
        ylabel('T_i [eV]')
end

if plot_magpres
    if exist(savename.ExB,"file")
        load(savename.ExB,'newPCBdata2D')
        mu0 = 4*pi*1E-7;
        idx_IDSP_z = knnsearch(newPCBdata2D.zq(1,:)',IDSP.z(1));
        idx_IDSP_z_min = knnsearch(newPCBdata2D.zq(1,:)',IDSP.z(1)-0.016);
        idx_IDSP_z_max = knnsearch(newPCBdata2D.zq(1,:)',IDSP.z(1)+0.016);
        idx_pcb_time = knnsearch(newPCBdata2D.trange',plot_time);
        magpres_z = newPCBdata2D.Bz(:,idx_IDSP_z,idx_pcb_time).^2/(2*mu0);
        magpres_r = newPCBdata2D.Br(:,idx_IDSP_z,idx_pcb_time).^2/(2*mu0);
        magpres_t = newPCBdata2D.Bt(:,idx_IDSP_z,idx_pcb_time).^2/(2*mu0);
        magpres_all = magpres_z + magpres_r + magpres_t;
        magpres_z_mean = mean(newPCBdata2D.Bz(:,idx_IDSP_z_min:idx_IDSP_z_max,idx_pcb_time-2:idx_pcb_time+2).^2/(2*mu0),3);
        magpres_z_mean = mean(magpres_z_mean,2);
        magpres_r_mean = mean(newPCBdata2D.Br(:,idx_IDSP_z_min:idx_IDSP_z_max,idx_pcb_time-2:idx_pcb_time+2).^2/(2*mu0),3);
        magpres_r_mean = mean(magpres_r_mean,2);
        magpres_t_mean = mean(newPCBdata2D.Bt_ext(:,idx_IDSP_z_min:idx_IDSP_z_max,idx_pcb_time-2:idx_pcb_time+2).^2/(2*mu0),3);
        magpres_t_mean = mean(magpres_t_mean,2);
        magpres_all_mean = magpres_z_mean + magpres_r_mean + magpres_t_mean;
        yyaxis right
        err_magpres_z_y = 2E1;
        err_magpres_r_y = 1E1;
        err_magpres_t_y = 3E3;
        err_magpres_x = 0.01;
        switch plot_time
            case {470,482}
                if magpres_mean
                    switch magpresplot
                        case 'Pz'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_z_mean,err_magpres_z_y,err_magpres_z_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pr'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_r_mean,err_magpres_r_y,err_magpres_r_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pt'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_t_mean,err_magpres_t_y,err_magpres_t_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pall'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_all_mean,err_magpres_t_y,err_magpres_t_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                    end
                else
                    switch magpresplot
                        case 'Pz'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_z,err_magpres_z_y,err_magpres_z_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pr'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_r,err_magpres_r_y,err_magpres_r_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pt'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_t,err_magpres_t_y,err_magpres_t_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                        case 'Pall'
                            errorbar(newPCBdata2D.rq(:,idx_IDSP_z),magpres_all,err_magpres_t_y,err_magpres_t_y,err_magpres_x,err_magpres_x,'g^','LineWidth',3)
                    end
                end
        end
        switch magpresplot
            case 'Pz'
                ylabel('B_Z^2/2\mu_0 [Pa]')
            case 'Pr'
                ylabel('B_R^2/2\mu_0 [Pa]')
            case 'Pt'
                ylabel('B_t^2/2\mu_0 [Pa]')
            case 'Pall'
                ylabel('B^2/2\mu_0 [Pa]')
        end
        hold on
        h = plot(nan,nan,'g^-','LineWidth',3);
        if isempty(hs)
            hs = h;
        else
            hs = [hs h];
        end
        if legendStrings == ""
            legendStrings = "Manetic Pressure";
        else
            legendStrings = [legendStrings "Manetic Pressure"];
        end
    else
        warning([savename.ExB, 'does not exist.'])
    end
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'g';
end

fontsize(20,"points")

xlim([0.07 0.27])

% xlabel('R [m]')
% switch profileplot
%     case 'Vz'
%         ylabel('V_Z [km/s]')
%     case 'Vr'
%         ylabel('V_R [km/s]')
% end
% legendStrings = string(WAV_IDSPdata.time(1:2)) +"us";
legend(hs,legendStrings,'Location','southwest')

% grid on


