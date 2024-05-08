NofCH = 32;%チャンネル数
NT = 1;%計測トリガ時刻数
NData = 1;%同一計測トリガでの計測ショット数
nz = 1;
factor = 0.3;
r_measured = zeros(NofCH/8,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/8,nz);%ベクトルプロットz座標1列目data1、2列目data2
a_measured = zeros(NofCH/8,nz);%ベクトルプロットa座標1列目data1、2列目data2
r_measured(:,1) = [10 15 20 25];%data1のr座標
z_measured(:,1) = 0;%data1のz座標
a_measured(:,1) = 0;%data1のa座標

Va_err = zeros(NofCH/8,nz);
Vz_err = zeros(NofCH/8,nz);
Vr_err = zeros(NofCH/8,nz);

figure('Position',[100 150 1500 400])
for t = 1:NT
    Ti_multi = zeros(NofCH/4,NData);
    for ndata = 1:NData
        load(['magflow/mat/3d_',num2str(462+4*t),'us_',num2str(ndata),'.mat'],'Va','Vz','Vr','Ti')
        Ti_multi(:,ndata) = Ti(1:NofCH/4,1);
    end
    A = [Ti_multi(1:2:8,:),Ti_multi(2:2:8,:)];
    B = [Vr(1:2:8,:),Vr(2:2:8,:)];
    if NData <= 2
        Ti_mean = trimmean(A,0,2);
        Va_mean = trimmean(Va,0,2);    
        Vz_mean = trimmean(Vz,0,2);
        Vr_mean = trimmean(B,0,2);
    else
        Ti_mean = trimmean(A,60,2);    
        Va_mean = trimmean(Va,60,2);    
        Vz_mean = trimmean(Vz,60,2);   
        Vr_mean = trimmean(B,60,2);
    end
    Ti_err = sqrt(var(A,0,2))/sqrt(NT*NData);%標準誤差
    Va_err = sqrt(var(Va,0,2))/sqrt(NT*NData);%標準誤差
    Vz_err = sqrt(var(Vz,0,2))/sqrt(NT*NData);%標準誤差
    Vr_err = sqrt(var(Vr,0,2))/sqrt(NT*NData);%標準誤差
    absV_mean = zeros(NofCH/8,1);
    if t == 1
        Ti_offset = Ti_mean;
    end
    for i = 1:NofCH/8
        absV_mean(i,1) = sqrt(Va_mean(i,1)^2 + Vz_mean(i,1)^2 + Vr_mean(i,1)^2);
    end
    
%     %平均流速、平均温度をプロット
%     subplot(1,NT,t);
%     %     Ticon = repmat(Ti_mean,1,2);%等高線図用平均イオン温度
%     Ticon = repmat(Ti_mean - Ti_offset,1,2);%等高線図用平均イオン温度上昇
%     zcon = z_measured(1,1);
%     s = pcolor([zcon-1.8 zcon+0.6],r_measured,Ticon);
%     s.FaceColor = 'interp';
%     s.FaceAlpha = 0.6;
%     s.EdgeAlpha = 0;
%     colormap('jet')
%         caxis([55 110])
%     caxis([0 50])
%     hold on


    subplot(1,NT,t)
    plot3(a_measured,z_measured,r_measured,'xr','MarkerSize',5,'LineWidth',2);
    hold on
    if nz == 1
        q = quiver3(a_measured,z_measured,r_measured,Va_mean*factor,Vz_mean*factor,Vr_mean*factor);
    end
    if nz == 2
        q = quiver3(a_measured,z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
        %V(:,1)などになってるのがよくわからないので要修正
    end
    
    q.LineWidth = 2;
    q.MaxHeadSize = 30;
    q.AutoScale = 'off';
    q.Color = 'black';
    xlim([min(a_measured,[],'all')-1 max(a_measured,[],'all')+1])
    ylim([min(z_measured,[],'all')-1 max(z_measured,[],'all')+1])
    zlim([min(r_measured,[],'all')-1 max(r_measured,[],'all')+1])
    %     title({[num2str(463+2*t),'us'];'Ion Flow [km/s]'},'Color','black','FontWeight','bold')
%     title([num2str(466+4*t),'us'],'Color','black','FontWeight','bold')
    absV_mean = round(absV_mean,2);
    for j = 1:nz
        for i = 1:NofCH/8
%             txt = text(a_measured(i,j)-1.0,z_measured(i,j)+1.0,r_measured(i,j)+0.5,num2str(absV_mean(i,j)));
            T_temp = round(Ti_mean(i,j));
            str1 = {num2str(T_temp)};
            str2 = {'[eV]'};
            txt = text(a_measured(i,j)-1.0,z_measured(i,j)+1.0,r_measured(i,j)-1.5,append(str1,str2));
            txt.FontSize = 14;
            txt.Color = 'k';
            txt.FontWeight = 'bold';
        end
    end
    xlabel('\theta [cm]')
    ylabel('Z [cm]')
    zlabel('R [cm]')
    ax = gca;
    ax.FontSize = 12;
    grid on
    daspect([0.5 0.5 1])
    
    %カラーマップをつける
    [X,Z] = meshgrid([-0.3,0.3],linspace(10,25,5));
    Y = ones(size(X));
    Ticon = [Ti_mean;Ti_mean(NofCH/8)];
    C = ([Ticon,Ticon]);
    [xq,zq] = meshgrid([-0.3,0.3],linspace(10,25,50));
    D = griddata(X,Z,C,xq,zq);
    yq = ones(size(xq));
    surf(xq,yq,zq,D,'EdgeColor','none')
    colormap('jet')
    if t == NT
        colorbar
        
    end
    
    
    for j = 1:nz
        for i = 1:NofCH/8
            str = {"V\theta = "+num2str(round(Va_mean(i,j),1))+'±'+num2str(round(Va_err(i),2)),"Vz = "+num2str(round(Vz_mean(i,j),1))+'±'+num2str(round(Vz_err(i),2)),"Vr = "+num2str(round(Vr_mean(i,j),1))+'±'+num2str(round(Vr_err(i),2)),"|V| = "+num2str(round(absV_mean(i,j)))+'[km/s]'};
            txt = text(a_measured(i,j)+0.3,z_measured(i,j)-1.2,r_measured(i,j)-0.2,str);
            txt.FontSize = 14;
            txt.Color = 'k';
            txt.FontWeight = 'bold';
        end
    end
    
end


% hp = get(subplot(1,NT,NT),'Position');
% c = colorbar('Position', [hp(1)+hp(3)+0.01  hp(2)+0.21  0.02  hp(3)*5.5]);
% % c.Label.String = 'Ion Temperature [eV]';
% c.Label.String = 'Increase of Ion Temperature [eV]';
% c.FontSize = 12;
