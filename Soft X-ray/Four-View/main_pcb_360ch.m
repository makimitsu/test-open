clearvars -except date IDXlist start dt
addpath 'C:/Users/yuleo/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath 'C:/Users/yuleo/Documents/GitHub/test-open'
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.processed_data = "G:\My Drive\X-ray\Data\PSIOUT\MATDATA";%processed data の保管場所
pathname.fig_psi = "G:\My Drive\X-ray\Data\PSIOUT\FIG"; % plot画像の保存先

% 実験オペレーション取得
prompt = {'Date:','Shot number:','start:','dt'};
definput = {'','','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('start','var')
    definput{3} = num2str(start);
end
if exist('dt','var')
    definput{4} = num2str(dt);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
start = str2num(cell2mat(answer(3)));
dt = str2num(cell2mat(answer(4)));

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';% スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);% 計測データ数
shotlist_a038 =T.a038(IDXlist);
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040, shotlist_a038];
tfshotlist_a038 =T.a038_TF(IDXlist);
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040,tfshotlist_a038];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
probecheck_mode = 0; % TFonly の時は必ずtrue(1)にして生信号をcheck
interp_method = 1;   % 0: 'scatteredInterpolant', 1: 'bz_rbfinterp', 2: 'spline'
calc_xpoint_mode = 1;
sheet_plot_mode = 0;
trange=400:800;% 計算時間範囲
n=70; % rz方向のメッシュ数
r_shift38 = 0.005; % プローブの差し込み具合を変更した場合は記入
r_shift39 = 0.000;
doCheck = false;

for processnum = 1:n_data

shot = shotlist(processnum,:);
tfshot=tfshotlist(processnum,:);

% MDSplusからの保存
filename3 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(38),'_shot',num2str(shot(3)),'_tfshot',num2str(tfshot(3)),'.mat');
if exist(filename3,"file")==0
    disp('No rawdata file of a038 -- Start generating!')
    rawdataPath = pathname.rawdata;
    save_dtacq_data(38, shot(3), tfshot(3), rawdataPath)
end
a038_raw = importdata(filename3);

filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
if exist(filename1,"file")==0
    disp('No rawdata file of a039 -- Start generating!')
    rawdataPath = pathname.rawdata;
    save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
end
a039_raw = importdata(filename1);

filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
if exist(filename2,"file")==0
    disp('No rawdata file of a040 -- Start generating!')
    rawdataPath = pathname.rawdata;
    save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
end
a040_raw = importdata(filename2);

% 3839.matの作成
filename3839 = strcat(pathname.processed_data,'\a03839_',num2str(shot(1)),'.mat');
if exist(filename3839,'file') == 0
    disp(strcat('File:',filename3839,' not found -> Initiate processing'));
    i_EF=EFlist(processnum);
    process_psi365ch(date, shot, tfshot, pathname, n,i_EF,trange, r_shift38,r_shift39,probecheck_mode,interp_method,calc_xpoint_mode)
end

load(filename3839);

if isstruct(grid2D) == 0 || isstruct(data2D) == 0
    disp('data incomplete/corrupted');
    return
end

% Psi図面の作成
figure('Position',[100 100 640 480],'visible','on');
set(gcf,'Name','shot' + string(shot(2)),'NumberTitle','off');

total = 6;

if sheet_plot_mode
    for m = 1:total
        grid2D = calc_sheet_width(grid2D, data2D,start,dt,total,m);
    end
end

tile = tiledlayout(2,3); % Figureを4行4列に分割。左上から右向きm番目の位置に図を描画
for m=1:total
    i = start + (m-1).*dt + 1 - 400;
    t = data2D.trange(i); % trange = 400:600
    nexttile

    colorplot = 'Jt';

    switch colorplot
        case 'Psi'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none'); clim([-8e-3,8e-3]); % psi
            label = 'Psi [Wb]';
        case 'Bz'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none'); clim([-0.1, 0.1]); % Bz
            label = 'Bz [T]';
        case 'Bt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-50e-3:0.2e-3:50e-3,'LineStyle','none'); clim([-0.06, 0.06]); % Bt
            label = 'Bt [T]';
        case 'Bt_ex'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt_ex(:,:,i),0.05:0.01:0.3,'LineStyle','none'); clim([0.05, 0.3]); % Bt_ex
            label = 'Bt_{ex} [T]';
        case 'Babs'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Babs(:,:,i),0:0.01:0.3,'LineStyle', 'none'); clim([0, 0.3]); % Babs
            label = '|B| [T]';
        case 'mg_P'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.mg_P(:,:,i),1e2:0.1e2:1.5e4,'LineStyle', 'none'); clim([1e2, 1.5e4]); % Babs
                label = 'B^2/(2mu)';
        case 'Jt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),30,'LineStyle','none'); clim([-0.7e+6, 0.7*1e+6]); % Jt
            label = 'Jt [A/m^2]';
            hold on
            % Jt等高線の表示
            contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),[2e5:2e5:10e5],'white','LineWidth',1)
            hold off
        case 'Jz'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jz(:,:,i),30,'LineStyle','none'); clim([-0.8*1e+6, 0.8*1e+6]); % Jt
            label = 'Jz [A/m^2]';
             % Jt等高線の表示
            contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),3.0e5:1.5e5:6e5,'white','LineWidth',1)
        case 'Jr'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jr(:,:,i),40,'LineStyle','none'); %clim([-8e-3,8e-3]); 
            label = 'Psi [Wb]';
        case 'Fabs'
            %contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Fabs(:,:,i),40,'LineStyle','none'); %clim([-8e-3,8e-3]); 
            label = '|j×B| [N]';
            hold on
            q = quiver(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.Fz(:,:,i)),squeeze(data2D.Fr(:,:,i)),1.5);
            q.Color = "k";
            q.LineWidth = 2;
            hold off
        case 'Et'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none'); clim([-500, 400]); % Et
            label = 'Et [V/m]';
        case 'eta'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.eta(:,:,i),20,'LineStyle','none'); clim([-10, 200]); % eta
            label = 'eta[Ω m]';
        case 'q'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.q(:,:,i),10,'LineStyle','none'); clim([-2, 2]); % q
            label = 'q';
        case 'p'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.P(:,:,i),20,'LineStyle','none'); clim([0, 3000]); % P
            label = 'P [Pa]';
        case 'dpdr'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.dPdr(:,:,i),10,'LineStyle','none'); %clim([-3000, 3000]); % dPdr
            label = 'dP/dr [Pa/m]';
        case 'Ip_area'
            zpos = grid2D.zq(1,:);
            rpos = grid2D.rq(:,1);
            % セパラトリクスとなるpsiを求める
            max_psi = zeros(size(zpos));
            for j = 1:size(zpos,2)
                max_psi(j) = max(data2D.psi(:,j,i));
            end
            max_psi_zng = max_psi(1:size(zpos,2)/2);
            max_psi_zpos = max_psi(size(zpos,2)/2+1:size(zpos,2));
            min_psi_zng = min(max_psi_zng(max_psi_zng>=0));
            min_psi_zpos = min(max_psi_zpos(max_psi_zpos>=0));
            for j = 1:size(zpos,2)
                if max_psi(j) == min_psi_zng
                    z_sep_idx_ng = j;
                elseif max_psi(j) == min_psi_zpos
                    z_sep_idx_pos = j;
                end
            end
            psi_sep = max(min_psi_zng,min_psi_zpos);
            for k = 1:size(rpos,1)
                r = rpos(k);
                for j = 1:size(zpos,2)
                    z = zpos(j);
                    if data2D.psi(k,j,i) >= psi_sep && j>=z_sep_idx_ng && j <= z_sep_idx_pos
                        plot(z,r,'r.','MarkerSize',6);
                        hold on
                    end
                end
            end
            % hold off
            label = 'None';
        otherwise
            disp('error: incorrect input for colorplot')
            return
    end
    
    % colormap(whitejet)
    colormap(turbo)
    axis image % 各軸のデータ単位を合わせる
    axis tight manual % 軸の範囲をデータ範囲と合わせて固定

    hold on
    % psi等高線の表示
    contour(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),-20e-3:0.9e-3:40e-3,'black','LineWidth',1.5) % psi
    %yline(grid2D.xpos(2,i))

    if sheet_plot_mode == 1
        xline(grid2D.lower_half_width(i),'b','LineWidth',1)
        xline(grid2D.upper_half_width(i),'b','LineWidth',1)
        yline(grid2D.lower_half_height(i),'b','LineWidth',1)
        yline(grid2D.upper_half_height(i),'b','LineWidth',1)
        
    end
    
    xlim([-0.0525 0.0525])
    xticks([-0.05 0 0.05])
    xticklabels({})
    yticklabels({})
    xtickangle(0)
    daspect([0.8 1 1]);
    hold off
    title(string(t) + 'us','FontSize',10)
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = label;
tile.TileSpacing = 'compact'; % 各tileの間隔を縮める
tile.Padding = 'compact'; % 各title内の余白を縮める

fontsize(gcf,25,'pixels')
% save plot image

foldername_psi = strcat(pathname.fig_psi, '/',num2str(date));
if exist(foldername_psi, 'dir') == 0
   mkdir(foldername_psi);
end
savename_psi = strcat(foldername_psi, '/shot', num2str(shot(2), '%04i'), '_', colorplot, '_unified_with_axis_',num2str(start),'to',num2str(t),'.png');
savename_psi_fig = strcat(foldername_psi, '/shot', num2str(shot(2), '%04i'), '_', colorplot, '_unified_with_axis_',num2str(start),'to',num2str(t),'.fig');
exportgraphics(gcf,savename_psi, 'Resolution',300);
saveas(gcf,savename_psi_fig);

end

function process_psi365ch(date, shot, tfshot, pathname, n,i_EF,trange, r_shift38,r_shift39,probecheck_mode,interp_method,calc_xpoint_mode)
%% read doi probe file
% 較正係数のバージョンを日付で判別
sheets_a39 = sheetnames('coeff200ch.xlsx');
sheets_a39 = str2double(sheets_a39);
sheet_date_a39 = max(sheets_a39(sheets_a39 <= date)); % 計測日以前で最新バージョンの較正係数を使用
C_a39 = readmatrix('coeff200ch.xlsx', 'Sheet', num2str(sheet_date_a39)); % 較正係数の読み込み
ok_a39 = logical(C_a39(:,14)); % chが生きていれば1，死んでいれば0
dtacq_num_list_a39 = C_a39(:,1);
dtacq_ch_a39= C_a39(:,2);
probe_num_list_a39 = C_a39(:,5);
probe_ch_a39 = C_a39(:,6);
polarity_a39=C_a39(:,13); % 極性
coeff_a39=C_a39(:,12); % 較正係数 RC/NS
zpos_a39=C_a39(:,9); % z位置[m]
rpos_a39=C_a39(:,10)+r_shift39; % r位置[m]
ch_a39=C_a39(:,7); % デジタイザch番号

%% read akimitsu calibration file

%較正係数のバージョンを日付で判別
sheets_a38 = sheetnames('coeff125ch.xlsx');
sheets_a38 = str2double(sheets_a38);
sheet_date_a38 = max(sheets_a38(sheets_a38 <= date)); % 計測日以前で最新バージョンの較正係数を使用
C_a38 = readmatrix('coeff125ch.xlsx', 'Sheet', num2str(sheet_date_a38)); % 較正係数の読み込み
ok_a38 = logical(C_a38(:,14)); % chが生きていれば1，死んでいれば0
dtacq_num_list_a38 = C_a38(:,1);
dtacq_ch_a38 = C_a38(:,2);
probe_num_list_a38 = C_a38(:,5);
probe_ch_a38 = C_a38(:,6);
probe_num_list_a38 = C_a38(:,5);
polarity_a38=C_a38(:,13); % 極性
coeff_a38=C_a38(:,12); % 較正係数 RC/NS
zpos_a38=C_a38(:,9); % z位置[m]
rpos_a38=C_a38(:,10)+r_shift38; % r位置[m]
ch_a38=C_a38(:,7); % デジタイザch番号
%disp(shot(1))
p_ch= readmatrix('coeff125ch.xlsx','Sheet','p_ch');



%% read rawdata file

if ismember(38,dtacq_num_list_a38)
    filename38 = strcat(pathname.rawdata,'rawdata_dtacq',num2str(38),'_shot',num2str(shot(3)),'_tfshot',num2str(tfshot(3)),'.mat');
    filename38_extf = strcat(pathname.rawdata,'rawdata_dtacq',num2str(38),'_shot',num2str(tfshot(3)),'_tfshot0.mat');
    if exist(filename38,"file")==0
        disp(['File: shot',num2str(shot(3)),' does not exist: start save_dtacqdata'])
    end
    if tfshot(3) ~= 0 && exist(filename38_extf,"file")==0
        disp(['File:',filename38_extf,' does not exit']);
        return
    end
    a038_raw = importdata(filename38);
    if tfshot(3) == 0
        a038_raw_extf = zeros(size(a038_raw));
    else
        a038_raw_extf = importdata(filename38_extf);
    end
end

if ismember(39,dtacq_num_list_a39)
    filename39 = strcat(pathname.rawdata,'rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    filename39_extf = strcat(pathname.rawdata,'rawdata_dtacq',num2str(39),'_shot',num2str(tfshot(1)),'_tfshot0.mat');
    if exist(filename39,"file")==0
        disp(['File: shot',num2str(shot(2)),' does not exist: start save_dtacqdata'])
    end
    if tfshot(1) ~= 0 && exist(filename39_extf,"file")==0
        disp(['File:',filename39_extf,' does not exit']);
        return
    end
    a039_raw = importdata(filename39);
    if tfshot(1) == 0
        a039_raw_extf = zeros(size(a039_raw));
    else
        a039_raw_extf = importdata(filename39_extf);
    end
end

if ismember(40,dtacq_num_list_a39)
    filename40 = strcat(pathname.rawdata,'rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    filename40_extf = strcat(pathname.rawdata,'rawdata_dtacq',num2str(40),'_shot',num2str(tfshot(2)),'_tfshot0.mat');
    if exist(filename40,"file")==0
        disp(['File:',filename40,' does not exist']);
        return
    end
    if tfshot(2)~= 0 && exist(filename40_extf,"file")==0
        disp(['File:',filename40_extf,' does not exist']);
        return
    end
    a40_raw = importdata(filename40);
    if tfshot(2) == 0
        a40_raw_extf = zeros(size(a40_raw));
    else
        a40_raw_extf = importdata(filename40_extf);
    end
end

%% read doi probe rawdata
raw_a39 = zeros(1000,length(dtacq_ch_a39));
raw_extf_a39 = zeros(1000,length(dtacq_ch_a39));
raw_a38 = zeros(1000,length(dtacq_ch_a38));
raw_extf_a38 = zeros(1000,length(dtacq_ch_a38));

%read doi probe rawdata
for i = 1:length(dtacq_ch_a39)
    if dtacq_num_list_a39(i) == 39
        raw_a39(:,i) = a039_raw(:,dtacq_ch_a39(i));
        raw_extf_a39(:,i) = a039_raw_extf(:,dtacq_ch_a39(i));
    elseif dtacq_num_list_a39(i) == 40
        raw_a39(:,i) = a40_raw(:,dtacq_ch_a39(i));
        raw_extf_a39(:,i) = a40_raw_extf(:,dtacq_ch_a39(i));
    end
end

% read akimitsu probe rawdata
for i = 1:length(dtacq_ch_a38)
    if dtacq_num_list_a38(i) == 38
        raw_a38(:,i) = a038_raw(:,dtacq_ch_a38(i));
        raw_extf_a38(:,i) = a038_raw_extf(:,dtacq_ch_a38(i));
    end
end

%% 今３本目と４本目のプローブを入れ替えている2023/08/31~

if date > 230831
    [real_probe3_d_ch,Locb3] = find(probe_num_list_a38==3);
    [real_probe4_d_ch,Locb4] = find(probe_num_list_a38==4);
    [real_probe5_d_ch,Locb5] = find(probe_num_list_a38==5);
    probe3_ch = probe_ch_a38(real_probe3_d_ch);
    probe4_ch = probe_ch_a38(real_probe4_d_ch);
    tmp = raw_a38;
    for i = 1:length(probe3_ch)
        %probe3のch1~ch25をそれぞれprobe4のch1~ch25に対応させてrawdataを入れ替える
        ind3 = find(probe3_ch == probe4_ch(i));
        raw_a38(:,real_probe3_d_ch(ind3)) = tmp(:,real_probe4_d_ch(i));
        ind4 = find(probe4_ch == probe3_ch(i));
        raw_a38(:,real_probe4_d_ch(ind4)) = tmp(:,real_probe3_d_ch(i));

        % real_probe3にはprobe3のd_chが入る
        clear ind3 ind4
    end
end


%% read b data a038 & a039
b_a38=raw_a38.*coeff_a38';%較正係数RC/NS
b_a38=b_a38.*polarity_a38';%極性揃え
b_extf_a38 = raw_extf_a38.*coeff_a38';
b_extf_a38 = b_extf_a38.*polarity_a38';

b_a39=raw_a39.*coeff_a39';%較正係数RC/NS
b_a39=b_a39.*polarity_a39';%極性揃え
b_extf_a39 = raw_extf_a39.*coeff_a39';
b_extf_a39 = b_extf_a39.*polarity_a39';

%% デジタイザchからプローブ通し番号順への変換 a038 & a039
bz_a39=zeros(1000,140);
bt=bz_a39;
bz_ex_a39 = bz_a39;
bt_ex = bz_a39;
ok_bz_a39=false(140,1);
ok_bt=ok_bz_a39;
zpos_bz_a39=zeros(140,1);
rpos_bz_a39=zpos_bz_a39;
zpos_bt=zpos_bz_a39;
rpos_bt=zpos_bz_a39;

%デジタイザchからプローブ通し番号順への変換 a038
bz_a38=zeros(1000,125);
bz_ex_a38 = bz_a38;
ok_bz_a38=false(125,1);
zpos_bz_a38=zeros(125,1);
rpos_bz_a38=zpos_bz_a38;

%% digital filter 移動平均フィルター（ノイズを含む信号の平滑化）
windowSize = 3;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

%% read bz bt data a038 & a039
for i=1:length(ch_a39)
    b_a39(:,i) = filter(bb,aa,b_a39(:,i));
    b_a39(:,i) = b_a39(:,i) - mean(b_a39(1:40,i));
    b_extf_a39(:,i) = filter(bb,aa,b_extf_a39(:,i));
    b_extf_a39(:,i) = b_extf_a39(:,i) - mean(b_extf_a39(1:40,i));

        if rem(ch_a39(i),2)==1
            bz_a39(:,ceil(ch_a39(i)/2))=b_a39(:,i);
            bz_ex_a39(:,ceil(ch_a39(i)/2)) = b_extf_a39(:,i);
            ok_bz_a39(ceil(ch_a39(i)/2))=ok_a39(i);
            zpos_bz_a39(ceil(ch_a39(i)/2))=zpos_a39(i);
            rpos_bz_a39(ceil(ch_a39(i)/2))=rpos_a39(i);
        elseif rem(ch_a39(i),2)==0
            bt(:,ch_a39(i)/2)=b_a39(:,i);
            bt_ex(:,ch_a39(i)/2) = b_extf_a39(:,i);
            ok_bt(ceil(ch_a39(i)/2))=ok_a39(i);
            zpos_bt(ceil(ch_a39(i)/2))=zpos_a39(i);
            rpos_bt(ceil(ch_a39(i)/2))=rpos_a39(i);
        end
end

nanlist = [find(isnan(b_a38(1,:))),63]; % ch 63, 127, 128
for i=1:length(ch_a38)
    b_a38(:,i) = filter(bb,aa,b_a38(:,i));
    b_a38(:,i) = b_a38(:,i) - mean(b_a38(1:40,i));
    b_extf_a38(:,i) = filter(bb,aa,b_extf_a38(:,i));
    b_extf_a38(:,i) = b_extf_a38(:,i) - mean(b_extf_a38(1:40,i));  
    if ~ismember(i, nanlist)
        bz_a38(:, ch_a38(i))=b_a38(:,i);
        bz_ex(:,ch_a38(i)) = b_extf_a38(:,i);
        ok_bz_a38(ch_a38(i))=ok_a38(i);
        zpos_bz_a38(ch_a38(i))=zpos_a38(i);
        rpos_bz_a38(ch_a38(i))=rpos_a38(i);
    end
end
bz_a38(:,63) = [];%1列詰められる
bz_ex_a38(:,63) = [];
ok_bz_a38(63) = [];
zpos_bz_a38(63) = [];
rpos_bz_a38(63) = [];


zprobepcb_a38    = [-0.0525 -0.021 0 0.021 0.0525]; % ch1-25, ch26-50, ch51-76, ch77-101, ch102-126
rprobepcb_a38    = [0.0600 0.0800 0.1000 0.1200 0.1400 0.1500...
                0.1550 0.1600 0.1650 0.1700 0.1750 0.1800...
                0.1850 0.1900 0.1950 0.2000 0.2050 0.2100... 
                0.2150 0.2200 0.2250 0.2300 0.2350 0.2400...
                0.2450]...
                +r_shift38;

zprobepcb_a39    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
rprobepcb_a39    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift39;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift39;

zprobepcb= [zprobepcb_a38,zprobepcb_a39];
rprobepcb = [rprobepcb_a38,rprobepcb_a39];
zpos_bz = vertcat(zpos_bz_a38,zpos_bz_a39);
rpos_bz = vertcat(rpos_bz_a38,rpos_bz_a39);
ok_bz = vertcat(ok_bz_a38,ok_bz_a39);
ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
ok_bt_matrix = false(length(rprobepcb_t),length(zprobepcb));
bz = horzcat(bz_a38,bz_a39);
bz_ex = horzcat(bz_ex_a38,bz_ex_a39);


[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));

%% probeの位置はakimitsu probeとdoi probeで分けてプロットする
[zq_probepcb_a38,rq_probepcb_a38]=meshgrid(zprobepcb_a38,rprobepcb_a38);
ok_bt_matrix_a38 = false(length(rprobepcb_a38),length(zprobepcb_a38));
ok_bz_matrix_a38 = false(length(rprobepcb_a38),length(zprobepcb_a38));

[zq_probepcb_a39,rq_probepcb_a39]=meshgrid(zprobepcb_a39,rprobepcb_a39);
ok_bt_matrix = false(length(rprobepcb_t),length(zprobepcb_a39));
ok_bz_matrix_a39 = false(length(rprobepcb_a39),length(zprobepcb_a39));


for i = 1:length(ok_bz_a38)
   
        index_r_a38 = (abs(rpos_bz_a38(i)-rprobepcb_a38)<0.001);index_z_a38 = (zpos_bz_a38(i)==zprobepcb_a38);
        ok_bz_matrix_a38 = ok_bz_matrix_a38 + rot90(index_r_a38,-1)*index_z_a38*ok_bz_a38(i);
     
end

for i = 1:length(ok_bz_a39)
   
        index_r_a39 = (abs(rpos_bz_a39(i)-rprobepcb_a39)<0.001);index_z_a39 = (zpos_bz_a39(i)==zprobepcb_a39);
        ok_bz_matrix_a39 = ok_bz_matrix_a39 + rot90(index_r_a39,-1)*index_z_a39*ok_bz_a39(i);
     
end

for i = 1:length(ok_bt)
    if rpos_bt(i) > (r_shift39)
        index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb_a39);
        ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ここまで編集済み
total = 9;
grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'rprobepcb_t',rprobepcb_t,...
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix,...
    'zprobepcb_a38',zprobepcb_a38,...
    'rprobepcb_a38',rprobepcb_a38,...
    'ok_bz_matrix_a38',ok_bz_matrix_a38,...
    'ok_bt_matrix_a38',ok_bt_matrix_a38,...
    'zprobepcb_a39',zprobepcb_a39,...
    'rprobepcb_a39',rprobepcb_a39,...
    'ok_bz_matrix_a39',ok_bz_matrix_a39,...
    'ok_bt_matrix_a39',ok_bt_matrix,...
    'lower_half_width',zeros(1,size(trange,2)),...
    'upper_half_width',zeros(1,size(trange,2)),...
    'lower_half_height',zeros(1,size(trange,2)),...
    'upper_half_height',zeros(1,size(trange,2)),...
    'width',zeros(n,total),...
    'height',zeros(n,total)...
    );
 
clear zq rq rq_t zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix 
clear sorted_zprobepcb sorted_rprobepcb sorted_zq_probepcb sorted_rq_probepcb sorted_ok_bz_matrix 
clear zprobepcb_a38 rprobepcb_a38 zq_probepcb_a38 rq_probepcb_a38 ok_bz_matrix_a38 
clear zprobepcb_a39 rprobepcb_a39 zq_probepcb_a39 rq_probepcb_a39 ok_bz_matrix_a39 


%% probecheck_script

if probecheck_mode
    probecheck_script;
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
    'mg_P',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'mg_Pz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Ft',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Fr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Fz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Fabs',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'eta',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'trange',trange,...
    'psi_pr',zeros(3,size(trange,2)),...
    'fitrate',zeros(1,size(trange,2)),...
    'xJt',zeros(1,size(trange,2)),...
    'xEt',zeros(1,size(trange,2)),...
    'xVr',zeros(1,size(trange,2)),...
    'xVz',zeros(1,size(trange,2)),...
    'mg_Px',zeros(1,size(trange,2)),... 
    'mg_Pz_x',zeros(1,size(trange,2)),...
    'xpos',zeros(2,size(trange,2)),... %X点座標(z,r)
    'pileup',zeros(1,total)...
    );


B_z_calibrated_restored = bz;
B_t_calibrated_restored = bt;
B_t_ex_calibrated_restored = bt_ex;

 for i=1:length(trange)
    t=trange(i);

    if interp_method == 1 %'bz_rbfinterp'
            %Bzの二次元補間(線形fit)
            vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
            B_z = -Bz_EF+vq;
            B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
            B_t_ex = bz_rbfinterp(rpos_bt-0.01, zpos_bt, grid2D, B_t_ex_calibrated_restored, ok_bt, t);

    elseif interp_method == 2 %'spline'
         vq = interp2(grid2D_probe.zq,grid2D_probe.rq,B_z_splined(:,:,i),grid2D.zq,grid2D.rq, 'spline');
         B_z = -Bz_EF+vq;
         B_t = interp2(grid2D_probe.zq,grid2D_probe.rq_t,B_t_splined(:,:,i),grid2D.zq,grid2D.rq, 'spline');
         B_t_ex = interp2(grid2D_probe.zq,grid2D_probe.rq_t,B_t_ex_splined(:,:,i),grid2D.zq,grid2D.rq, 'spline');   
    end

            % PSI計算
            % data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
            data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);

            %このままだと1/2πrが計算されてないので
            [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
            data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
            data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
            data2D.Bt(:,:,i)=B_t;
            data2D.Bt_ex(:,:,i)=B_t_ex;
            data2D.Babs(:,:,i) = sqrt(data2D.Br(:,:,i).^2 + data2D.Bz(:,:,i).^2 + (data2D.Bt(:,:,i)+data2D.Bt_ex(:,:,i)).^2); 

            mu=4*pi*1e-7;
            data2D.mg_P(:,:,i) = (data2D.Babs(:,:,i).^2)./(2*mu);
                        
            data2D.Bt_total(:,:,i) = (data2D.Bt(:,:,i)+data2D.Bt_ex(:,:,i));
            %Jr, jzをrotBから計算
            [data2D.Jr(:,:,i),data2D.Jz(:,:,i)]=gradient(data2D.Bt_total(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1));
            data2D.Jr(:,:,i)=-data2D.Jr(:,:,i)./mu;
            data2D.Jz(:,:,i)=data2D.Jz(:,:,i)./mu;
            data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
            data2D.Fr(:,:,i) = data2D.Jt(:,:,i).*data2D.Bz(:,:,i) - data2D.Bt_total(:,:,i).*data2D.Jz(:,:,i);
            data2D.Ft(:,:,i) = data2D.Jz(:,:,i).*data2D.Br(:,:,i) - data2D.Jr(:,:,i).*data2D.Bz(:,:,i);
            data2D.Fz(:,:,i) = data2D.Jr(:,:,i).*data2D.Bt_total(:,:,i) - data2D.Jt(:,:,i).*data2D.Br(:,:,i);
            data2D.Fabs(:,:,i) = sqrt(data2D.Fr(:,:,i).^2 + data2D.Fz(:,:,i).^2 + data2D.Ft(:,:,i).^2); 
            data2D.Lambda(:,:,i) = (2*pi*grid2D.rq.*data2D.Bt(:,:,i))./(data2D.psi(:,:,i));
            %data2D.dPdr(:,:,i) = (data2D.Jt*data2D.Bz - data2D.Jz*data2D.Bt);
          
 end

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

clearvars -except data2D grid2D shot date;
filename = strcat('G:\My Drive\X-ray\Data\PSIOUT\MATDATA\a03839_',num2str(shot(1)),'.mat');%保存ファイル名
save(filename)
end

function save_dtacq_data(dtacq_num,shot,tfshot,rawdataPath)

[rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
if tfshot==0
    [rawdata]=getMDSdata(dtacq_num,shot,0);
    save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
end

end