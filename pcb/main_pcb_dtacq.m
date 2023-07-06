%%%%%%%%%%%%%%%%%%%%%%%%
%pcbプローブ(bz)の解析
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
%pathname.ts3u='/Users/mgar/koala_home/pub/mnt/old-koala/experiment/results/ts-3u/';%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier='/Users/mgar/data/';%fourierのmd0（データックのショットが入ってる）までのpath
%pathname.NIFS='';%resultsまでのpath（ドップラー、SXR）
pathname.save='/Users/mgar/pcb_save'; %保存先

%%%%(3)指定したshotの解析
dtacqlist=3122:3149; %【input】テーブルから解析したいshot番号を抽出して入力

date = 211126;
d_tacqTF = 3123;
i_EF = 150;
trange=460:490;
n=50; %rz方向のメッシュ数

for d_tacq=dtacqlist(1,7)
plot_psi(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname); %通常の時系列プロット
%plot_position(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname); %計測位置、各位置での生信号も含めた確認用プロット
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, plot_psi, plot_posision)
%%%%%%%%%%%%%%%%%%%%%%%%

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット
%Github/test-open/pcbdata.mを使用->[チャンネルごとの生信号のプロット]の項目をコメントアウトしなければ、別figureで確認用にプロットされる
function plot_psi(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname)

[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF,trange, [], n, i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
    maxrange=max(abs(data2D.Jt),[],'all');

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2); %各r,tでのpsiの最小値,時間
[opoint,p]=islocalmin(psimid,1); %全rでのpsiの極小値
[xpoint,~]=islocalmax(psimid,1); %全rでのpsiの極大値
[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%磁気面時間発展プロット
f=figure;
f.WindowState = 'maximized';
 start=0; %460+?
 t_start=460+start;
 for m=1:10 %図示する時間
     i=start+m; %end
     t=trange(i);
     subplot(2,5,m)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet) %jet/parula
    axis image
    axis tight manual
    %     xlim([-0.02 0.02])
    %     ylim([0.12 0.27])
    caxis([-10*1e+6,10*1e+6]) %カラーバーの軸の範囲
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
    plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
    hold off
    title(string(t)+'us')
    xlabel('z')
    ylabel('r')
end

sgtitle(strcat('date=',num2str(date),': dtacq=',num2str(d_tacq)));

%figureの保存
%saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_time',num2str(t_start),'.png'))
%save(strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'.mat'))

%close
end

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット+測定プローブ位置の表示
%Github/test-open/pcbdata.mを使用->[チャンネルごとの生信号のプロット]の項目をコメントアウトしなければ、別figureで確認用にプロットされる
%Github/test-open/getpcbbz.mを使用
function plot_position(date, d_tacq, d_tacqTF, trange, n, i_EF, ~)

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
 start=11; %460+?
 for m=1:10 %図示する時間
     i=start+m; %end
     t=trange(i);
     subplot(2,5,m)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    %     xlim([-0.02 0.02])
    %     ylim([0.12 0.27])
    %caxis([-1.8*1e+6,1.8*1e+6])
    caxis([-maxrange,maxrange])
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
    xlabel('z')
    ylabel('r')
end

sgtitle(strcat('date=',num2str(date),': dtacq=',num2str(d_tacq)))

%saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.png'))
%save(strcat(pathname.save,'\IDX',num2str(IDX),'_shot',num2str(date),num2str(shot,'%03i'),'_dtacq',num2str(T.d_tacq(IDX)),'_cr',num2str(cr_time),'us.mat')
%close
end