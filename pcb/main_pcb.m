%%%%%%%%%%%%%%%%%%%%%%%%
%pcbプローブ(bz)の解析
%%%%%%%%%%%%%%%%%%%%%%%%

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
pathname.ts3u='ts3u_path';%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier='fourier_path';%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS='NIFS_path';%resultsまでのpath（ドップラー、SXR）
pathname.save='C:\Users\uswk0\OneDrive\デスクトップ\data\a038_out\'; %保存先

%%%%(2)ログから解析したいデータを検索
%Github/test-open/searchlog.mを使用

node='date';  % 【input】検索する列の名前. T.Properties.VariableNamesで一覧表示できる
 pat=211223;   % 【input】検索パターン（数値なら一致検索、文字なら含む検索）　

searchlog(T,node,pat); % ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%%%%(3)指定したshotの解析
IDXlist=2870:2921; %【input】テーブルから解析したいshot番号を抽出して入力
%for IDX=IDXlist(1,37)
%plot_psi(T, pathname,IDX); %通常の時系列プロット
%plot_position(T, pathname, IDX); %計測位置、各位置での生信号も含めた確認用プロット
%end
IDX = 2911;
plot_psi(T, pathname,IDX);
%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, plot_psi, plot_posision)
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%getinput:実験ログ（T）からshot番号（IDX）におけるオペレーションの値を出力
%出力は構造体に変更しても良いかも
function [date, shot, TF_shot, offset_TF, i_EF, start, Doppler_t, d_tacq, d_tacqTF, trange, t, n] = getinput(T,IDX)
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
offset_TF=isfinite(TF_shot);

if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end

start=T.Period_StartTime_(IDX);
Doppler_t=T.DopplerDelay(IDX);

d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);

trange=460:490;
t=T.DopplerDelay(IDX);
n=50; %rz方向のメッシュ数
end


%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット
%Github/test-open/pcbdata.mを使用->[チャンネルごとの生信号のプロット]の項目をコメントアウトしなければ、別figureで確認用にプロットされる
function plot_psi(T, pathname,IDX)
[date, shot, TF_shot, offset_TF, i_EF, start, Doppler_t, d_tacq, d_tacqTF,trange, t, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力

[grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
    maxrange=max(abs(data2D.Jt),[],'all');

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,p]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%磁気面時間発展プロット
f=figure;
f.WindowState = 'maximized';
 start=13; %460+?
 t_start=460+start;
 for m=1:10 %図示する時間
     i=start+m*2; %end
     t=trange(i);
     subplot(2,5,m)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    %     xlim([-0.02 0.02])
    %     ylim([0.12 0.27])
    caxis([-3*1e+6,3*1e+6]) %カラーバーの軸の範囲
    %caxis([-maxrange,maxrange])
    colorbar('Location','eastoutside')
    %zlim([-1 1])
    %colormap(bone)
    %%カラーバーのラベル付け
    %c = colorbar;
    %c.Label.String = 'Jt';
    hold on
    plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
    %plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
    hold off
    title(string(t)+'us','Fontsize', 12)
    xlabel('z [m]','FontSize',12)
    ylabel('r [m]','FontSize',12)
 end
 
 %sgtitle(strcat('IDX=',num2str(IDX),': shot=',num2str(date),num2str(shot,'%03i'),': dtacq=',num2str(T.d_tacq(IDX))))

%figureの保存
saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_time',num2str(t_start),'.png'))
%saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.png'))
%saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'.png'))
%save(strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.mat'))
%save(strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'.mat'))
%close
end

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット+測定プローブ位置の表示
%Github/test-open/pcbdata.mを使用->[チャンネルごとの生信号のプロット]の項目をコメントアウトしなければ、別figureで確認用にプロットされる
%Github/test-open/getpcbbz.mを使用
function plot_position(T, ~, IDX)
[date, shot, ~, ~, i_EF, ~, ~, d_tacq, d_tacqTF,trange, ~, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力

load('rc_coeff2020.mat'); 
[rawdata]=getvalue(d_tacq,d_tacqTF); % rawdata
[ok, bz, rpos, zpos, ~] = getpcbbz(rawdata, coeff, date);
x = zpos(ok); %z方向の生きているチャンネル
y = rpos(ok); %r方向の生きているチャンネル

[grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
    maxrange=max(abs(data2D.Jt),[],'all');

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,p]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%磁気面時間発展プロット+計測位置プロット
figure
 start=9; %460+?
 for m=1:10 %図示する時間
     i=start+(m*2+1); %end
     t=trange(i);
     subplot(2,5,m)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    %     xlim([-0.02 0.02])
    %     ylim([0.12 0.27])
    caxis([-3*1e+6,3*1e+6])
    %caxis([-maxrange,maxrange])
    colorbar('Location','eastoutside')
    %zlim([-1 1])
    %colormap(bone)
    hold on
    plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black','LineWidth',0.6)
    plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
    plot(x,y,"k.",'MarkerSize', 10)
    hold off
    title(string(t)+'us')
    xlabel('z [m]','FontSize',12)
    ylabel('r [m]','FontSize',12)
end

sgtitle(strcat('IDX=',num2str(IDX),': shot=',num2str(date),num2str(shot,'%03i'),': dtacq=',num2str(T.d_tacq(IDX))))

saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.png'))
%save(strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.mat')
%close
end