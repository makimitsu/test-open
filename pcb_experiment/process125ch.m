%%%%%%%%%%%%%%%%%%%%%%%%
% 125ch用 akimitsu pcbプローブのみでの磁気面（Bz）
% dtacqのshot番号を直接指定する場合
% 補間手順（scatteredInterpolant)
% 測定データ→死んだチャンネルを補間→プロット用meshgridに補間
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
clear all
pathname.ts3u=getenv('ts3u_path');% old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');% fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');% resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');% outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');% dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');% rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');% dtacqのrawdataの保管場所


%% 直接入力の場合
%{
shotlist=[11946]; %[10650:10692];%【input】dtacqの保存番号: shot_38
%tfshot=10646*ones(size(shotlist));
%tfshot= [10646];
tfshot=zeros(size(shotlist));%【input】dtacqのTFのみ番号: tfshot_38
date = 230914;%【input】計測日
i_EF = 120;%【input】EF電流
probecheck_mode = 1; % 【input】TF only の時は必ずtrue(1)にして生信号をcheck
interp_method = 1; % 0: 'scatteredInterpolant', 1: 'bz_rbfinterp', 2: 'spline'
trange=400:800;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数
r_shift = 0.00; %【input】プローブの差し込み具合を変更した場合は記入
%{
​probecheck_mode=1;%【input】TF only の時は必ず1にして生信号をcheck
interp_method = 1; %【input】0: 'scatteredInterpolant', 1: 'bz_rbfinterp', 2: 'spline'
%}
for i = 1:length(shotlist)
    shot = shotlist(i)
    process_psi125ch(date,shot,tfshot,pathname,n,i_EF,trange,r_shift,probecheck_mode,interp_method);
    disp('pcb:1/1')
end
%}


%%　実験オペレーションの取得　自動入力の場合
prompt = {'Date:','Shot number:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2double(cell2mat(answer(2)));
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);%計測データ数
shotlist_a038 =T.a038(IDXlist);
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a038, shotlist_a039, shotlist_a040];
tfshotlist_a038 =T.a038_TF(IDXlist);
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a038, tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=38.*ones(n_data,1);

probecheck_mode = 1; %【input】TF only の時は必ずtrue(1)にして生信号をcheck
interp_method = 1;   % 【input】0: 'scatteredInterpolant', 1: 'bz_rbfinterp', 2: 'spline'
trange=400:800;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数
r_shift = 0; % 【input】プローブの差し込み具合を変更した場合は記入​
for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    process_psi125ch(date,shot,tfshot,pathname,n,i_EF,trange,r_shift,probecheck_mode,interp_method);
    disp('pcb:1/1')
    plot_psi125ch(shot,trange);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function process_psi125ch(date, shot, tfshot, pathname, n,i_EF,trange, r_shift,probecheck_mode,interp_method)

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff125ch.xlsx');
sheets = str2double(sheets);
sheet_date = max(sheets(sheets <= date)); % 計測日以前で最新バージョンの較正係数を使用
C = readmatrix('coeff125ch.xlsx', 'Sheet', num2str(sheet_date)); % 較正係数の読み込み
ok = logical(C(:,14)); % chが生きていれば1，死んでいれば0
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
probe_num_list = C(:,5);
polarity=C(:,13); % 極性
coeff=C(:,12); % 較正係数 RC/NS
zpos=C(:,9); % z位置[m]
rpos=C(:,10)+r_shift; % r位置[m]
ch=C(:,7); % デジタイザch番号
%disp(shot(1))
p_ch= readmatrix('coeff125ch.xlsx','Sheet','p_ch');

if shot(1) == tfshot(1)
        tfshot(1) = 0;
end

if ismember(38,dtacq_num_list)
    filename38 = strcat(pathname.rawdata,'\rawdata_dtacq',num2str(38),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    filename38_extf = strcat(pathname.rawdata,'rawdata_dtacq',num2str(38),'_shot',num2str(tfshot(1)),'_tfshot0.mat');
    if exist(filename38,"file")==0
        disp(['File: shot',num2str(shot(1)),' does not exist: start save_dtacqdata'])
       % res = save_dtacqdata_func(shot,tfshot);
    end
    if tfshot(1) ~= 0 && exist(filename38_extf,"file")==0
        disp(['File:',filename38_extf,' does not exit']);
        return
    end
    a038_raw = importdata(filename38);
    if tfshot(1) == 0
        a038_raw_extf = zeros(size(a038_raw));
    else
        a038_raw_extf = importdata(filename38_extf);
    end
end


raw = zeros(1000,length(dtaq_ch));
raw_extf = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 38
        raw(:,i) = a038_raw(:,dtaq_ch(i));
        raw_extf(:,i) = a038_raw_extf(:,dtaq_ch(i));
    end
end

%% 今３本目と４本目のプローブを入れ替えている2023/08/31~

if date > 230831
    [real_probe3_d_ch,Locb3] = find(probe_num_list==3);
    [real_probe4_d_ch,Locb4] = find(probe_num_list==4);
    % Lia3にはprobe3の要素にd_chが入る
    tmp = raw;
    raw(:,real_probe3_d_ch) = tmp(:,real_probe4_d_ch);
    raw(:,real_probe4_d_ch) = tmp(:,real_probe3_d_ch);
end


b=raw.*coeff';%較正係数RC/NS
b=b.*polarity';%極性揃え
b_extf = raw_extf.*coeff';
b_extf = b_extf.*polarity';

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,125);
bz_ex = bz;
ok_bz=false(125,1);
zpos_bz=zeros(125,1);
rpos_bz=zpos_bz;


%digital filter 移動平均フィルター（ノイズを含む信号の平滑化）
windowSize = 3;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

% ch: pcb coilのch number
% dtacq_ch: digitizerのch,つまり保存される磁場信号bzの番号
nanlist = [find(isnan(b(1,:))),63]; % ch 63, 127, 128

for i=1:length(ch)
    
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    b_extf(:,i) = filter(bb,aa,b_extf(:,i));
    b_extf(:,i) = b_extf(:,i) - mean(b_extf(1:40,i));
    
    
    if ~ismember(i, nanlist)
        bz(:, ch(i))=b(:,i);
        bz_ex(:,ch(i)) = b_extf(:,i);
        ok_bz(ch(i))=ok(i);
        zpos_bz(ch(i))=zpos(i);
        rpos_bz(ch(i))=rpos(i);
    end
 
end
bz(:,63) = [];%1列詰められる
bz_ex(:,63) = [];
ok_bz(63) = [];
zpos_bz(63) = [];
rpos_bz(63) = [];

zprobepcb    = [-0.0525 -0.021 0 0.021 0.0525]; % ch1-25, ch26-50, ch51-76, ch77-101, ch102-126
rprobepcb    = [0.0600 0.0800 0.1000 0.1200 0.1400 0.1500...
                0.1550 0.1600 0.1650 0.1700 0.1750 0.1800...
                0.1850 0.1900 0.1950 0.2000 0.2050 0.2100... 
                0.2150 0.2200 0.2250 0.2300 0.2350 0.2400...
                0.2450]...
                +r_shift;

[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
[zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);

ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
for i = 1:length(ok_bz)
    index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
    ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
end

grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'ok_bz_matrix',ok_bz_matrix);
grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb);

clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb ok_bz_matrix

%% probecheck_script

if probecheck_mode
   %probecheck_script125ch;
   probecheck_script_test;
end

%% calculate data

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
    'Bt_ex',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Babs',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'trange',trange);

%{
%% **************** angle correction **************** 
% 死んだchの内挿 → プローブ1本ごとに角度補正 → プロット用gridへの内挿
% 内挿は Delaunay 三角形分割による
sheets_angle = sheetnames('angle.xlsx');sheets_angle = str2double(sheets_angle);
sheet_angle_date=max(sheets_angle(sheets_angle<=date));
angle_file = readmatrix('angle.xlsx','Sheet',num2str(sheet_angle_date));
angle_zpos = angle_file(1,2:end);
angle = angle_file(2,2:end);
sin_lis = sin(angle);
cos_lis = cos(angle);
sin_matrix = repmat(sin(angle),10,1);%sin_matrix(1:10,:) = zeros(10,10);
cos_matrix = repmat(cos(angle),10,1);%cos_matrix(1:10,:) = ones(10,10);

% B_z_calibrated = zeros(length(grid2D.rprobepcb),length(grid2D.zprobepcb),length(trange));
% B_t_calibrated = zeros(length(grid2D.rprobepcb_t),length(grid2D.zprobepcb),length(trange));
% B_t_ex_calibrated = zeros(length(grid2D.rprobepcb_t),length(grid2D.zprobepcb),length(trange));

bz_ok = bz(:,ok_bz);
zpos_bz_ok = zpos_bz(ok_bz);
rpos_bz_ok = rpos_bz(ok_bz);
bz_ex_ok = bz_ex(:,ok_bz);

B_z_noncalib = zeros(length(trange),size(bz,2));
B_z_ex_noncalib = B_z_noncalib;

B_z_calibrated_restored = bz;


Fz = scatteredInterpolant(zpos_bz_ok, rpos_bz_ok, bz_ok(1,:)');
Fz_ex = scatteredInterpolant(zpos_bz_ok, rpos_bz_ok, bz_ex_ok(1,:)');

for i=1:length(trange)
    t=trange(i);
    Fz.Values = bz_ok(t,:)';
    Fz_ex.Values = bz_ex_ok(t,:)';
    B_z_noncalib(i,:) = Fz(zpos_bz, rpos_bz)';
    B_z_ex_noncalib(i,:) = Fz_ex(zpos_bz, rpos_bz)';
end


% 
% for i=1:length(angle) % プローブ1本(10ch)ごとに角度較正
%     idx = 10*(i-1)+1:10*i; % プローブのidxリスト
%     cos_val = cos_lis(angle_zpos == zpos_bz(10*i));
%     sin_val = sin_lis(angle_zpos == zpos_bz(10*i));
%     B_z_calibrated_restored(trange,idx) = B_z_noncalib(:,idx).*cos_val - B_t_noncalib(:,idx).*sin_val;
%     B_t_calibrated_restored(trange,idx) = B_z_noncalib(:,idx).*sin_val + B_t_noncalib(:,idx).*cos_val;
%     B_t_ex_calibrated_restored(trange,idx) = B_z_ex_noncalib(:,idx).*sin_val + B_t_ex_noncalib(:,idx).*cos_val;
% end



B_z_splined = zeros(length(grid2D.rprobepcb),length(grid2D.zprobepcb),length(trange));

for ch = 1:size(B_z_noncalib,2)
    B_z_splined(grid2D.rprobepcb==rpos_bz(ch),grid2D.zprobepcb==zpos_bz(ch),:) = B_z_calibrated_restored(trange,ch);
end
%}


%% **************** interpolation and grid 2D calculation **************** 
% 3 methods:scatterInterpolant or bz_rbfinterp or spline
% bz_rbfinterp: linear, multiquadric...



 %Fz_grid = scatteredInterpolant(zpos_bz, rpos_bz, B_z_calibrated_restored(1,:)');
 

 for i=1:length(trange)
    t=trange(i);
    
    if interp_method == 0 %'scatteredInterpolant'
            Fz_grid.Values = B_z_calibrated_restored(t,:)';
          

            grid_z_lis = reshape(grid2D.zq,[],1);
            grid_r_lis = reshape(grid2D.rq,[],1);
            vq = reshape(Fz_grid(grid_z_lis,grid_r_lis),[n,n]);
            B_z = -Bz_EF + vq;
  
        

    elseif interp_method == 1 %'bz_rbfinterp'
            %Bzの二次元補間(線形fit)
            vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
            B_z = -Bz_EF+vq;
            

    elseif interp_method == 2 %'spline'
         vq = interp2(grid2D_probe.zq,grid2D_probe.rq,B_z_splined(:,:,i),grid2D.zq,grid2D.rq, 'spline');
         B_z = -Bz_EF+vq;
        
    end

            

            % PSI計算
            data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
            %data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);

            %このままだと1/2πrが計算されてないので
            [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
            data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
            data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
            %data2D.Bt(:,:,i)=B_t;
            %data2D.Bt_ex(:,:,i)=B_t_ex;
            data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
            %data2D.Lambda(:,:,i) = (2*pi*grid2D.rq.*data2D.Bt(:,:,i))./(data2D.psi(:,:,i));
            %data2D.Babs(:,:,i) = sqrt(data2D.Br(:,:,i).^2 + data2D.Bz(:,:,i).^2 + (data2D.Bt(:,:,i)+data2D.Bt_ex(:,:,i)).^2);
 
 end
% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

clearvars -except data2D grid2D shot date;
filename = strcat("G:\My Drive\X-ray\Data\A038FILE\shot",num2str(shot(1)),'.mat');%保存ファイル名
save(filename)
end

function plot_psi125ch(shot,trange)

load(strcat("G:\My Drive\X-ray\Data\A038FILE\shot\",num2str(shot),".mat"));

figure('Position', [0 0 1500 1500],'visible','on');
start=15;
% t_start=470+start;
 for m=1:10 % 図示する時間
     i=start+m.*2; %end
     t=trange(i);
     subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    caxis([-2*1e+6,1.0*1e+6]) %カラーバーの軸の範囲
%     caxis([-0.1,0.1])
    colorbar('Location','eastoutside')
    % カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
    plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    title(string(t)+' us')
    xlabel('z [m]')
    ylabel('r [m]')
 end

end
