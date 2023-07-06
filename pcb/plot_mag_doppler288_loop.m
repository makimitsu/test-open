clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%
%  288chDopplerとmrdの磁気面を重ねてプロットして保存するコード
%　288chDopplerはIDLで作ったsavのデータをあらかじめfourier\md0\makimitsu\date(yyddmm)に保存してあるものを読み込む
%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath


setenv('a038_path','I:\a038'); %fourier/md0/a038までのパスをa038_pathとして環境変数設定して、内部ネットで実行できる
yourname = 'C:\Users\Moe Akimitsu\';
f = fullfile(yourname,'Documents','GitHub','test-open');
addpath(genpath(f));
%addpath(fullfile(yourname,'Documents','GitHub','test-open','pcb'));
%f = fullfile(yourname,'Documents','GitHub','SXR_test');
%addpath(f);

%%%適宜変更


DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得
 
%shotlist=[2946:2950];%IDX（ログの通し番号を入れる）
shotlist=[2946];

%subT=searchlog(T,'date',211224);
%shotlist=subT.number;
subT=T(shotlist,:);
IDXlist=subT{:,1};
IDXlist=IDXlist(isfinite(subT.DopplerDelay)&isfinite(subT.d_tacq))';
%IDX=IDXlist(1,88);

for IDX=IDXlist 

date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end
t=T.DopplerDelay(IDX);
n=50;


%%%DopplerはIDLのファイル

%共有フォルダ以下から目的ショットのファイルを探す
%filepath.rgw=strcat(pathname.ts3u, '\', string(date),'\' ...
%    ,string(date),num2str(shot,'%03i'),'.rgw');
%filepath.D288=dir(strcat(pathname.fourier,'\Doppler\288CH\20',string(date),'\*shot',num2str(shot),'*.asc'));
if shot<10
    filepath.D288=dir(strcat(pathname.fourier,'\makimitsu\',string(date),'\doppler2D_shot',num2str(shot),'_.sav'));
else
    filepath.D288=dir(strcat(pathname.fourier,'\makimitsu\',string(date),'\doppler2D_shot',num2str(shot),'.sav'));
end    
    %filepath.Dhighspeed=dir(strcat(pathname.NIFS,'\Doppler\Photron\',string(date),'\**\*shot',num2str(shot),'*.tif'));
    filepath.SXR=strcat(pathname.NIFS,'\X-ray\',string(date),'\shots\',string(date),num2str(shot,'%03i'),'.tif');
if numel(filepath.D288)==0
    continue
end
 f = fullfile(pathname.ts3u,num2str(date),strcat(num2str(date),num2str(shot,'%03i'),'.mrd'));
if isfile(f)==0
    continue
end

if isfile( fullfile(filepath.D288.folder,filepath.D288.name))
    [B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num] = get_B_z(date,TF_shot,shot,true,i_EF,pathname.ts3u);

    B_z = B_z([2,3,4,6,7,8],2:end,:);
    data = data([2,3,4,6,7,8],2:end,:);
    z_probe = z_probe(2:end);
    ch_dist = ch_dist([2,3,4,6,7,8],2:end);
    r_probe = r_probe([2,3,4,6,7,8]);

    psi = get_psi(B_z,r_probe,t);
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
    [probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    %restore_idl( fullfile(filepath.D288.folder,filepath.D288.name),'lowercase','create'); %.savファイルを読み込む
    filenameasc=dir(strcat('I:\Doppler\288CH\20',num2str(date),'\shot',num2str(shot),'_*.asc'));
    if numel(filenameasc)==0
        continue
    end    
    doppler=doppler288ch(fullfile(filenameasc.folder,filenameasc.name),date);


    f = figure;
    % f.WindowState = 'maximized';
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.6];
    pos1 = [0.07,0.2,0.35,0.6];
    pos2 = [0.58,0.2,0.35,0.6];

    subplot('Position',pos1);
    contourf(doppler.z,doppler.yy,doppler.emission,[0:0.1:1]*3e5,'LineStyle','none')
    colormap(jet)
    colorbar('Location','eastoutside')
    hold on
    contour(psi_mesh_z,psi_mesh_r,psi,30,'-b');
    hold off
    axis image
    axis tight manual
    title(string(t)+'us,emiision')
    xlabel('z')
    ylabel('r')
    ylim([0.075 0.25])
    caxis([0,3e5])
   % doppler.ti_2d(doppler.ti_2d>100)=NaN;
    subplot('Position',pos2);
    contourf(doppler.z,doppler.yy,doppler.ti_2d,[0:0.05:1]*100,'LineStyle','none')
    colormap(jet)
    colorbar('Location','eastoutside')
    hold on
    contour(psi_mesh_z,psi_mesh_r,psi,30,'-b');
    hold off
    axis image
    axis tight manual
    caxis([0,100])
     ylim([0.075 0.25])
    title(string(t)+'us,ti')
    xlabel('z')
    ylabel('r')
    filename = strcat('I:\makimitsu\',num2str(date),'\Doppler_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
% plot_psi_SXR_at_t(B_z,r_probe,z_probe,date,shot,473,false,false,T.Period_StartTime_(IDX),5,false,filepath.SXR);
     saveas(gcf,strcat(filename,'_all','.png'))
    close
end


% %[B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] 
% [B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num]= get_B_z(date,TF_shot,shot,offset_TF,T.EF_A_(IDX),pathname.ts3u);
% B_z = B_z([2,3,4,6,7,8],2:end,:);
% data = data([2,3,4,6,7,8],2:end,:);
% z_probe = z_probe(2:end);
% ch_dist = ch_dist([2,3,4,6,7,8],2:end);
% r_probe = r_probe([2,3,4,6,7,8]);



%plot_B_z_in_time(B_z,ch_dist,350,600);
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,layer,area,start,exposure,SXRfilename)
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,true,true,T.Period_StartTime_(IDX),2,filepath.SXR)
end
%

