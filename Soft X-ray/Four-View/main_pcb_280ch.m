clear
% close all
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%%%%%%%%%%%%%%%%%%%%
%280ch用新規pcbプローブのみでの磁気面（Bz）
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
% date=230714;
T=searchlog(T,node,date);
% IDXlist= 1; %[5:50 52:55 58:59];%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
n_data=numel(IDXlist);%計測データ数
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

trange=400:800;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数

doCheck = false;
doCheck = true;
if ~doCheck
    figure('Position', [0 0 1500 1500],'visible','on');
end

for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    if doCheck
        check_signal(date, shot, tfshot, pathname, n);
    else
        plot_psi200ch(date, shot, tfshot, pathname,n,i_EF,trange);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi200ch(date, shot, tfshot, pathname, n,i_EF,trange)
% filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% if exist(filename,"file")==0
%     disp('No rawdata file -- Start generating!')
%     rawdataPath = pathname.rawdata;
%     save_dtacq_data(dtacq_num, shot, tfshot,rawdataPath)
%     % return
% end
% load(filename,'rawdata');%1000×192
% 
% %正しくデータ取得できていない場合はreturn
% if numel(rawdata)< 500
%     return
% end

filename = strcat(pathname.pre_processed_directory,'/a039_',num2str(shot(1)),'.mat');
if exist(filename,'file') == 0
    doCalculation = true;
    disp('no processed data -- start calculation');
else
    doCalculation = false;
    disp('loading processed data');
end


if doCalculation
    %較正係数のバージョンを日付で判別
    sheets = sheetnames('coeff200ch.xlsx');
    sheets = str2double(sheets);
    sheet_date=max(sheets(sheets<=date));
    C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
    r_shift = 0.00;
    ok = logical(C(:,14));
    dtacq_num_list = C(:,1);
    dtaq_ch = C(:,2);
    polarity=C(:,13);
    coeff=C(:,12);
    zpos=C(:,9);
    rpos=C(:,10)+r_shift;
    ch=C(:,7);
    
    if ismember(39,dtacq_num_list)
        filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
        if exist(filename1,"file")==0
            disp('No rawdata file of a039 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
            % disp(['File:',filename1,' does not exit']);
            % return
        end
        a039_raw = importdata(filename1);
    end
    if ismember(40,dtacq_num_list)
        filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
        if exist(filename2,"file")==0
            disp('No rawdata file of a040 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
            % disp(['File:',filename2,' does not exit']);
            % return
        end
        a040_raw = importdata(filename2);
    end
    
    raw = zeros(1000,length(dtaq_ch));
    for i = 1:length(dtaq_ch)
        if dtacq_num_list(i) == 39
            raw(:,i) = a039_raw(:,dtaq_ch(i));
        elseif dtacq_num_list(i) == 40
            raw(:,i) = a040_raw(:,dtaq_ch(i));
        end
    end
    
    b=raw.*coeff';%較正係数RC/NS
    b=b.*polarity';%極性揃え
    
    %デジタイザchからプローブ通し番号順への変換
    bz=zeros(1000,100);
    bt=bz;
    ok_bz=false(100,1);
    ok_bt=ok_bz;
    zpos_bz=zeros(100,1);
    rpos_bz=zpos_bz;
    zpos_bt=zpos_bz;
    rpos_bt=zpos_bz;
    
    %digital filter
    windowSize = 8;
    bb = (1/windowSize)*ones(1,windowSize);
    aa = 1;
    
    for i=1:length(ch)
        b(:,i) = filter(bb,aa,b(:,i));
        b(:,i) = b(:,i) - mean(b(1:40,i));
        if rem(ch(i),2)==1
            bz(:,ceil(ch(i)/2))=b(:,i);
            ok_bz(ceil(ch(i)/2))=ok(i);
            zpos_bz(ceil(ch(i)/2))=zpos(i);
            rpos_bz(ceil(ch(i)/2))=rpos(i);
        elseif rem(ch(i),2)==0
            bt(:,ch(i)/2)=b(:,i);
            ok_bt(ceil(ch(i)/2))=ok(i);
            zpos_bt(ceil(ch(i)/2))=zpos(i);
            rpos_bt(ceil(ch(i)/2))=rpos(i);
        end
    end
    
    % zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
    zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
    rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
    rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
    [zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
    % [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
    [zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
    ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
    ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
    for i = 1:length(ok_bt)
        if rpos_bt(i) > (r_shift)
            index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
            ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
        end
        index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
        ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
    end
    
    grid2D=struct(...
        'zq',zq,...
        'rq',rq,...
        'zprobepcb',zprobepcb,...
        'rprobepcb',rprobepcb,...
        'rprobepcb_t',rprobepcb_t,...1
        'ok_bz_matrix',ok_bz_matrix,...
        'ok_bt_matrix',ok_bt_matrix);
    grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);
    
    clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix
    
    % probecheck_script;
    
    %data2Dcalc.m
    r_EF   = 0.5 ;
    n_EF   = 234. ;
    
    if date<221119
        z1_EF   = 0.68;
        z2_EF   = -0.68;
    else
        z1_EF   = 0.78;
        z2_EF   = -0.78;
    end
    [Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
    clear EF r_EF n_EF i_EF z_EF
    
    data2D=struct(...
        'psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
        'trange',trange);
    
    % ******************* no angle correction ********************
    for i=1:size(trange,2)
        t=trange(i);
    
        %Bzの二次元補間(線形fit)
        vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
        B_z = -Bz_EF+vq;
        B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
    
        % PSI計算
        data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
        % data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);
        % このままだと1/2πrが計算されてないので
        [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
        data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
        data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
        data2D.Bt(:,:,i)=B_t;
        data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
    end
else
    load(filename,'data2D','grid2D');
end
% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

% プロット部分
% figure('Position', [0 0 1500 1500],'visible','on');
start=60;
dt = 4;
%  t_start=470+start;
 for m=1:16 %図示する時間
     i=start+m.*dt; %end
     t=trange(i);
     subplot(4,4,m)
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
%     caxis([-0.01,0.01])%Bz
     % clim([-0.1,0.1])%Bt
    % clim([-5e-3,5e-3])%psi
%     caxis([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
     % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
    hold off
    title(string(t)+' us')
%     xlabel('z [m]')
%     ylabel('r [m]')
 end

 drawnow

if doCalculation
    clearvars -except data2D grid2D shot pathname;
    filename = strcat(pathname.pre_processed_directory,'/a039_',num2str(shot(1)),'.mat');
    save(filename)
end

end

function save_dtacq_data(dtacq_num,shot,tfshot,rawdataPath)

% rawdataPath=getenv('rawdata_path'); %保存先

% %%%%実験オペレーションの取得
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
% node='date';
% % pat=230707;
% pat=date;
% T=searchlog(T,node,pat);
% IDXlist=6;%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
% IDXlist = shot;
% n_data=numel(IDXlist);%計測データ数
% shotlist=T.a039(IDXlist);
% tfshotlist=T.a039_TF(IDXlist);
% EFlist=T.EF_A_(IDXlist);
% TFlist=T.TF_kV_(IDXlist);
% dtacqlist=39.*ones(n_data,1);

% dtacqlist=39;
% shotlist= 1819;%【input】dtacqの保存番号
% tfshotlist = 1817;

% dtacqlist=40;
% shotlist= 299;%【input】dtacqの保存番号
% tfshotlist = 297;

% date = 230706;%【input】計測日
% n=numel(shotlist);%計測データ数

%RC係数読み込み

% for i=1:n
%     dtacq_num=dtacqlist(i);
%     shot=shotlist(i);
%     tfshot=tfshotlist(i);
%     [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
%     save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
%     if tfshot>0
%         [rawdata]=getMDSdata(dtacq_num,shot,0);
%         % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
%         save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
%     end
% end

% dtacq_num=dtacqlist(i);
% shot=shotlist(i);
% tfshot=tfshotlist(i);

[rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
if tfshot==0
    [rawdata]=getMDSdata(dtacq_num,shot,0);
    % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
    save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
end

end

function check_signal(date, shot, tfshot, pathname,n)

sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
r_shift = 0.00;
ok = logical(C(:,14));
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
polarity=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10)+r_shift;
ch=C(:,7);

if ismember(39,dtacq_num_list)
    filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    if exist(filename1,"file")==0
        disp('No rawdata file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
    a039_raw = importdata(filename1);
end
if ismember(40,dtacq_num_list)
    filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    if exist(filename2,"file")==0
        disp('No rawdata file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
        % disp(['File:',filename2,' does not exit']);
        % return
    end
    a040_raw = importdata(filename2);
end

raw = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 39
        raw(:,i) = a039_raw(:,dtaq_ch(i));
    elseif dtacq_num_list(i) == 40
        raw(:,i) = a040_raw(:,dtaq_ch(i));
    end
end

b=raw.*coeff';%較正係数RC/NS
b=b.*polarity';%極性揃え
% b = smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

%digital filter
windowSize = 8;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:length(ch)
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
for i = 1:length(ok_bt)
    if rpos_bt(i) > (r_shift)
        index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
        ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
    end
    index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
    ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
end

grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'rprobepcb_t',rprobepcb_t,...1
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix);

figure_switch = ["on","on","on","on","off","off"];%bz1,bz2,bt1,bt2,bz_vs_z,bt_vs_z

r = 7;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.05;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.05;%縦軸プロット領域（b_z下限）
t_start=1;%横軸プロット領域（開始時間）
t_end=1000;%横軸プロット領域（終了時間）
t = 470;
% z_probe_pcb = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
z_probe_pcb = [-0.2975 -0.255 -0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17 0.255 0.2975];
n_z = length(z_probe_pcb);
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure(Visible=figure_switch(1));
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j));
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        ylim([y_lower_lim y_upper_lim]);
        % xticks([t_start t_end]);
        xticks([]);
        yticks([]);
    end
end
sgtitle('Bz signal probe1-5') 

f2=figure(Visible=figure_switch(2));
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        ylim([y_lower_lim y_upper_lim]);
        % xticks([t_start t_end]);
        xticks([]);
        yticks([]);
    end
end
sgtitle('Bz signal probe6-10')

f3=figure(Visible=figure_switch(3));
f3.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)));
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
        % xticks([t_start t_end]);
        xticks([]);
        yticks([]);
    end
end
sgtitle('Bt signal probe1-5')

f4=figure(Visible=figure_switch(4));
f4.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)));
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
        % xticks([t_start t_end]);
        xticks([]);
        yticks([]);
    end
end
sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

% 横軸z, 縦軸Bzのプロット
f5=figure(Visible=figure_switch(5));
f5.WindowState = 'maximized';
styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(z_probe_pcb,bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bz_interped = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bz_zline=Bz_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

% 横軸z, 縦軸Btのプロット
f6=figure(Visible=figure_switch(6));
f6.WindowState = 'maximized';
styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=bt(t,zline);
    bt_zline(ok_bt(zline)==false)=NaN;
    plot(z_probe_pcb,bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bt_interped = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=Bt_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

hidden = find(figure_switch == 'off');
figures = [f1,f2,f3,f4,f5,f6];
for i = hidden
    close(figures(i));
end

end