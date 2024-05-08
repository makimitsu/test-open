%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%プロット枚数、プロット開始時間、プロット時間間隔を
%指定して指定したz断面(z=cut_z[cm])のBrをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_Br200ch(cut_z_Br,n_plot,t_start,dt,trange,grid2D,data2D)
cut_z_Br = cut_z_Br*1e-2;%単位を[cm]から[m]に変換
IDX = knnsearch(grid2D.zq(1,:).',cut_z_Br);%grid2D.zqの中で最もcut_zに近いセル番号を取得
z = grid2D.zq(1,IDX)*1e2;%z座標[cm]

fig = figure('Position', [0 0 1500 1500],'visible','on');
for m=1:n_plot %図示する時間
    i=(t_start-trange(1)+1)+(m-1)*dt; %end
    t=trange(i);
    if n_plot == 1
    elseif n_plot == 4
        subplot(2,2,m);
    elseif mod(n_plot,4) == 0
        subplot(n_plot/4,4,m)
    elseif n_plot < 4
        subplot(1,n_plot,m)
    else
        close(fig)
        warning('Input error in n_plot.')%n_plotの入力エラー
        return;
    end
    plot(grid2D.rq(:,1),data2D.Br(:,IDX,i)*1e3)

    title([num2str(t),' us, z = ',num2str(z),'cm'])
    xlabel('r [m]')
    ylabel('Br [mT]')
end