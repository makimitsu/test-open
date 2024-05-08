NofCH = 32;%チャンネル数
NT = 4;%計測トリガ時刻数
NData = 4;%同一計測トリガでの計測ショット数
nz = 1;
factor = 0.6;
overlap = 1;  %1:重ねない,2:重ねる，3:横軸Rでプロット
PT = 2; %1:PF,2:TF
r_measured = zeros(NofCH/8,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/8,nz);%ベクトルプロットz座標1列目data1、2列目data2
a_measured = zeros(NofCH/8,nz);%ベクトルプロットa座標1列目data1、2列目data2
r_measured(:,1) = [10 15 20 25];%data1のr座標
z_measured(:,1) = 0;%data1のz座標
a_measured(:,1) = 0;%data1のa座標

Va_err = zeros(NofCH/8,nz);
Vz_err = zeros(NofCH/8,nz);
Vr_err = zeros(NofCH/8,nz);

%それぞれのグラフの幅
Tiyl = -30;
Tiyh = 250;
a = -8;
b = 8;
Vryl = a;
Vryh = b;
Vzyl = a;
Vzyh = b;
Vtyl = a;
Vtyh = b;

figure('Position',[100 150 500 700])
% figure('Position',[100 150 500 700])
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

    R = linspace(10,25,4);
    
    if overlap == 1
%         subplot(4,NT,t)
%         errorbar(R,Ti_mean,Ti_err)
%         title('Ti','Color','black','FontWeight','bold')
%         xlabel('R [cm]')
%         ylabel('Ti [eV]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([8 27]);
%         ylim([Tiyl Tiyh]);
        
        subplot(3,NT,t)
        errorbar(R,Vr_mean,Vr_err)
        title('Vr','Color','black','FontWeight','bold')
        xlabel('R [cm]')
        ylabel('Vr [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([8 27]);
        ylim([Vryl Vryh]);
        
        
        subplot(3,NT,t+1*NT)
        errorbar(R,Vz_mean,Vz_err)
        title('Vz','Color','black','FontWeight','bold')
        xlabel('R [cm]')
        ylabel('Vz [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([8 27]);
        ylim([Vzyl Vzyh]);
        
        subplot(3,NT,t+2*NT)
        errorbar(R,Va_mean,Va_err)
        title('V\theta','Color','black','FontWeight','bold')
        xlabel('R [cm]')
        ylabel('V\theta [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([8 27]);
        ylim([Vtyl Vtyh]);
        
    elseif overlap == 2
%         ax1 = subplot(4,1,1);
%         errorbar(R,Ti_mean,Ti_err)
%         set(ax1,'NextPlot','add');
%         title('Ti','Color','black','FontWeight','bold')
%         xlabel('R [cm]')
%         ylabel('Ti [eV]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([8 27]);
%         ylim([Tiyl Tiyh]);
        
        ax2 = subplot(3,1,1);
        errorbar(R,Vr_mean,Vr_err,"LineWidth",1.5)
        set(ax2,'NextPlot','add');
        title('Vr','Color','black','FontWeight','bold')
        if t == NT
            if PT == 1
            legend('PF = 20kV','PF = 25kV','PF = 30kV','PF = 35kV','PF = 39kV','Location','best')
            else
            legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
            end
        end
        xlabel('R [cm]')
        ylabel('Vr [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([9 26]);
        ylim([Vryl Vryh]);
        
        
        ax3 = subplot(3,1,2);
        errorbar(R,Vz_mean,Vz_err,"LineWidth",1.5)
        set(ax3,'NextPlot','add');
        title('Vz','Color','black','FontWeight','bold')
        if t == NT
            if PT == 1
            legend('PF = 39kV','PF = 35kV','PF = 30kV','PF = 25kV','PF = 20kV','Location','best')
            else
            legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
            end
        end
        xlabel('R [cm]')
        ylabel('Vz [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([9 26]);
        ylim([Vzyl Vzyh]);
        
        ax4 = subplot(3,1,3);
        errorbar(R,Va_mean,Va_err,"LineWidth",1.5)
        set(ax4,'NextPlot','add');
        title('V\theta','Color','black','FontWeight','bold')
        if t == NT
            if PT == 1
            legend('PF = 39kV','PF = 35kV','PF = 30kV','PF = 25kV','PF = 20kV','Location','best')
            else
            legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
            end
        end
        xlabel('R [cm]')
        ylabel('V\theta [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([9 26]);
        ylim([Vtyl Vtyh]);
        
    elseif overlap == 3
        
        subplot(3,NT,t)
        errorbar(Vr_mean,R,Vr_err,'horizontal',"-o")
        title('Vr','Color','black','FontWeight','bold')
        ylabel('R [cm]')
        xlabel('Vr [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([Vryl Vryh]);
        ylim([8 27]);
        
        
        subplot(3,NT,t+1*NT)
        errorbar(Vz_mean,R,Vz_err,'horizontal',"-o")
        title('Vz','Color','black','FontWeight','bold')
        ylabel('R [cm]')
        xlabel('Vz [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([Vzyl Vzyh]);
        ylim([8 27]);
        
        subplot(3,NT,t+2*NT)
        errorbar(Va_mean,R,Va_err,'horizontal',"-o")
        title('V\theta','Color','black','FontWeight','bold')
        ylabel('R [cm]')
        xlabel('V\theta [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        xlim([Vtyl Vtyh]);
        ylim([8 27]);
    end
end


