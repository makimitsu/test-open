%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%プロット枚数、プロット開始時間、プロット時間間隔を
%指定して指定したz断面(z=cut_z[cm])の磁気圧をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_magpres200ch(cut_z_magpres,n_plot,t_start,dt,trange,grid2D,data2D,data2D_tfoffset)
cut_z_magpres = cut_z_magpres*1e-2;%単位を[cm]から[m]に変換
IDX = knnsearch(grid2D.zq(1,:).',cut_z_magpres);%grid2D.zqの中で最もcut_zに近いセル番号を取得
z = grid2D.zq(1,IDX)*1e2;%z座標[cm]
mu0 = 4*pi*1e-7;%真空の透磁率
[nR,~,~] = size(data2D.Bz);
magpres_z = zeros(nR,1);
magpres_t = zeros(nR,1);
magpres = zeros(nR,1);

% fig = figure('Position', [0 0 1500 1500],'visible','on');
fig = figure('Position', [0 0 400 500],'visible','on');
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
    for j = 1:nR
        magpres_z(j,1) = data2D.Bz(j,IDX,i)^2/(2*mu0);
        magpres_t(j,1) = data2D_tfoffset.Bt(j,IDX,i)^2/(2*mu0);
        magpres(j,1) = magpres_z(j,1) + magpres_t(j,1);
    end
    plot(grid2D.rq(:,1),magpres_z(:,1),'LineWidth',3)
    hold on
    plot(grid2D.rq(:,1),magpres_t(:,1),'LineWidth',3)
    hold on
    plot(grid2D.rq(:,1),magpres(:,1),'LineWidth',3)

    % legend('$\frac{B_z^2}{2\mu_0}$',Interpreter="latex")
    lgd =legend('$\frac{B_z^2}{2\mu_0}$','$\frac{B_t^2}{2\mu_0}$','$\frac{B_z^2+B_t^2}{2\mu_0}$',Interpreter="latex");
    lgd.FontSize = 20;
    % title(sprintf('%d us, z = %.1f cm',t,z))
    title(sprintf('%d [us]',t))
    xlabel('r [m]')
    ylabel('Magnetic Pressure [Pa]')
    xlim([0.15 0.25])
    ax = gca;
    ax.FontSize = 16;
end