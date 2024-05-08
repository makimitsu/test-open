NofCH = 32;%チャンネル数
NT = 4;%計測トリガ時刻数
NData = 4;%同一計測トリガでの計測ショット数
nz = 1;
factor = 0.7;
theta_ratio = 2;
theta_base = 2;
r_measured = zeros(NofCH/8,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/8,nz);%ベクトルプロットz座標1列目data1、2列目data2
a_measured = zeros(NofCH/8,nz);%ベクトルプロットa座標1列目data1、2列目data2
r_measured(:,1) = [10 15 20 25];%data1のr座標
% r_measured(:,1) = [12.5 17.5 22.5 27.5];%data1のr座標
% r_measured(:,1) = [17.5 22.5 27.5 32.5];%data1のr座標
% z_measured(:,1) = 0;%data1のz座標
% z_measured(:,1) = [-0.3 -0.3 -0.3 -0.3 0 0 0 0 0.3 0.3 0.3 0.3];%data1のz座標
% z_measured(:,1) = [13 13 13 13 0 0 0 0 -13 -13 -13 -13];%data1のz座標
% r_measured(:,1) = [10 15 20 25 10 15 20 25 10 15 20 25];%data1のr座標
a_measured(:,1) = 0;%data1のa座標

Va_err = zeros(NofCH/8,nz);
Vz_err = zeros(NofCH/8,nz);
Vr_err = zeros(NofCH/8,nz);

figure('Position',[400 150 1200 400])
for t = 1:NT
    Ti_multi = zeros(NofCH/4,NData);
    for ndata = 1:NData
        load(['magflow/mat/3d_',num2str(462+4*t),'us_',num2str(ndata),'.mat'],'Va','Vz','Vr','Ti')
        Ti_multi(:,ndata) = Ti(1:NofCH/4,1);
    end
    A = [Ti_multi(1:2:NofCH/4,:),Ti_multi(2:2:NofCH/4,:)];
    B = [Vr(1:2:NofCH/4,:),Vr(2:2:NofCH/4,:)]
    if NData <= 2
        Ti_mean = trimmean(A,0,2);
        Va_mean = trimmean(Va,50,2);    
        Vz_mean = trimmean(Vz,50,2);
        Vr_mean = trimmean(B,50,2);
    else
        Ti_mean = trimmean(A,60,2);    
        Va_mean = trimmean(Va,60,2);    
        Vz_mean = trimmean(Vz,60,2);   
        Vr_mean = trimmean(B,60,2);
    end
    Ti_err = sqrt(var(A,0,2))/sqrt(NT*NData);%標準誤差
    Va_err = sqrt(var(Va,0,2))/sqrt(NT*NData);%標準誤差
    Vz_err = sqrt(var(Vz,0,2))/sqrt(NT*NData);%標準誤差
    Vr_err = sqrt(var(B,0,2))/sqrt(NT*NData);%標準誤差
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
%     plot(z_measured,r_measured,'xr','MarkerSize',2,'LineWidth',2);

    hold on
    if nz == 1
        q = quiver(z_measured,r_measured,Vz_mean*factor,Vr_mean*factor);
    end
    if nz == 2
        q = quiver(z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
        %V(:,1)などになってるのがよくわからないので要修正
    end
    
    q.LineWidth = 2;
    q.MaxHeadSize = 30;
    q.AutoScale = 'off';
    q.Color = 'red';
    xlim([min(z_measured,[],'all')-5 max(z_measured,[],'all')+5])
    ylim([min(r_measured,[],'all')-7.5 max(r_measured,[],'all')+7])
%         title({[num2str(463+2*t),'us'];'Ion Flow [km/s]'},'Color','black','FontWeight','bold')
        title({'Ion Flow [km/s]'},'Color','black','FontWeight','bold')
%     title([num2str(462+4*t),'us'],'Color','black','FontWeight','bold')
    absV_mean = round(absV_mean,2);
%     for j = 1:nz
%         for i = 1:NofCH/8
% %             txt = text(a_measured(i,j)-1.0,z_measured(i,j)+1.0,r_measured(i,j)+0.5,num2str(absV_mean(i,j)));
%             T_temp = round(Ti_mean(i,j));
%             str1 = {num2str(T_temp)};
%             str2 = {'[eV]'};
%             txt = text(z_measured(i,j)-5.0,r_measured(i,j)-0.5,append(str1,str2));
%             txt.FontSize = 12;
%             txt.Color = 'k';
%             txt.FontWeight = 'bold';
%         end
%     end
    xlabel('Z [cm]')
    ylabel('R [cm]')
    ax = gca;
    ax.FontSize = 12;
    grid off
    daspect([0.5 0.5 1])
    
    
    
    for j = 1:nz
        for i = 1:NofCH/8
            str = {"Vr = "+num2str(round(Vr_mean(i,j),1)),"Vz = "+num2str(round(Vz_mean(i,j),1)),"V\theta = "+num2str(round(Va_mean(i,j),1))};
%             str = {"Vz = "+num2str(round(Vz_mean(i,j),1)),"Vr = "+num2str(round(Vr_mean(i,j),1)),"|V| = "+num2str(round(absV_mean(i,j)))+'[km/s]'};
%             str = {"Vz = "+num2str(round(Vz_mean(i,j),1))+'±'+num2str(round(Vz_err(i),2)),"Vr = "+num2str(round(Vr_mean(i,j),1))+'±'+num2str(round(Vr_err(i),2)),"|V| = "+num2str(round(absV_mean(i,j)))+'[km/s]'};
            txt = text(z_measured(i,j)+2,r_measured(i,j)-0.0,str);
            txt.FontSize = 10;
            txt.Color = 'k';
%             txt.FontWeight = 'bold';
            if round(Va_mean(i,j),1) > 0
                a = round(Va_mean(i,j)*theta_ratio)+theta_base;
                plot(z_measured(i,j),r_measured(i,j),'ob','MarkerSize',a,'LineWidth',1);
            elseif round(Va_mean(i,j),1) < 0
                a = (round(-1*Va_mean(i,j)*theta_ratio)+theta_base);
                plot(z_measured(i,j),r_measured(i,j),'xb','MarkerSize',a,'LineWidth',1);
                plot(z_measured(i,j),r_measured(i,j),'ob','MarkerSize',a,'LineWidth',1);
            end
        end
    end
    
%     % *** 軸などの削除は任意 ***
%     ax = gca;
%     ax.XColor = 'none';
%     ax.YColor = 'none';
%     ax.Color = 'none';
%     ax.Box = 'off';
%     ax.AmbientLightColor = 'none';
%     % *************************


    
end
% SaveTransparentFig(fig);

