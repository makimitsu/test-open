NofCH = 32;%チャンネル数
NT = 4;%計測トリガ時刻数
NData = 4;%同一計測トリガでの計測ショット数
nz = 1;
factor = 0.6;
flag = 2;%1:PF,2:TF
r_measured = zeros(NofCH/8,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/8,nz);%ベクトルプロットz座標1列目data1、2列目data2
a_measured = zeros(NofCH/8,nz);%ベクトルプロットa座標1列目data1、2列目data2
r_measured(:,1) = [10 15 20 25];%data1のr座標
z_measured(:,1) = 0;%data1のz座標
a_measured(:,1) = 0;%data1のa座標

Va_err = zeros(NofCH/8,nz);%Vtの誤差
Vz_err = zeros(NofCH/8,nz);%Vzの誤差
Vr_err = zeros(NofCH/8,nz);%Vrの誤差

%時間ごとに全ての情報を入れていく
Va_time = zeros(NofCH/8,NT);
Vz_time = zeros(NofCH/8,NT);
Vr_time = zeros(NofCH/8,NT);
Va_err_time = zeros(NofCH/8,NT);%時間ごとのVtの誤差
Vz_err_time = zeros(NofCH/8,NT);%時間ごとのVzの誤差
Vr_err_time = zeros(NofCH/8,NT);%時間ごとのVrの誤差

%それぞれのグラフの幅
Tiyl = -30;
Tiyh = 250;
Vryl = -20;
Vryh = 11;
Vzyl = -8;
Vzyh = 8;
Vtyl = -5;
Vtyh = 9;

figure('Position',[100 150 800 300])
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
    
    Vr_time(:,t) = Vr_mean;
    Vz_time(:,t) = Vz_mean;
    Va_time(:,t) = Va_mean;
    Vr_err_time(:,t) = Vr_err;
    Vz_err_time(:,t) = Vz_err;
    Va_err_time(:,t) = Va_err;
    
    
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
    PF = linspace(39,33,3);
    TF = [0.2095029 0.2078232 0.2005398 0.1896198];
    xneg = [0.000260837 0.002252045 0.002200151 0.000283717];
    xpos = [0.000260837 0.002252045 0.002200151 0.000471567];
    
%     %         ax1 = subplot(4,1,1);
%     %         errorbar(R,Ti_mean,Ti_err)
%     %         set(ax1,'NextPlot','add');
%     %         title('Ti','Color','black','FontWeight','bold')
%     %         xlabel('R [cm]')
%     %         ylabel('Ti [eV]')
%     %         ax = gca;
%     %         ax.FontSize = 12;
%     %         grid on
%     %
%     %         xlim([8 27]);
%     %         ylim([Tiyl Tiyh]);
%     
%     Vr_time(t,:)
%     
%     ax2 = subplot(3,1,1);
%     errorbar(R,Vr_mean,Vr_err,"LineWidth",1.5)
%     set(ax2,'NextPlot','add');
%     title('Vr','Color','black','FontWeight','bold')
%     if t == NT
%         legend('PF = 39kV','PF = 36kV','PF = 33kV','Location','best')
%         %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%     end
%     xlabel('R [cm]')
%     ylabel('Vr [km/s]')
%     ax = gca;
%     ax.FontSize = 12;
%     grid on
%     
%     xlim([9 26]);
%     ylim([Vryl Vryh]);
%     
%     
%     ax3 = subplot(3,1,2);
%     errorbar(R,Vz_mean,Vz_err,"LineWidth",1.5)
%     set(ax3,'NextPlot','add');
%     title('Vz','Color','black','FontWeight','bold')
%     if t == NT
%         legend('PF = 39kV','PF = 36kV','PF = 33kV','Location','best')
%         %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%     end
%     xlabel('R [cm]')
%     ylabel('Vz [km/s]')
%     ax = gca;
%     ax.FontSize = 12;
%     grid on
%     
%     xlim([9 26]);
%     ylim([Vzyl Vzyh]);
%     
%     ax4 = subplot(3,1,3);
%     errorbar(R,Va_mean,Va_err,"LineWidth",1.5)
%     set(ax4,'NextPlot','add');
%     title('V\theta','Color','black','FontWeight','bold')
%     if t == NT
%         legend('PF = 39kV','PF = 36kV','PF = 33kV','Location','best')
%         %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%     end
%     xlabel('R [cm]')
%     ylabel('V\theta [km/s]')
%     ax = gca;
%     ax.FontSize = 12;
%     grid on
%     
%     xlim([9 26]);
%     ylim([Vtyl Vtyh]);
end

if flag == 1%PFの場合
    for t = 1:NofCH/8%-1はRが大きい点を無視してる（光量が足りてなさそうで信頼性に欠けるため）
        
        ax1 = subplot(1,1,1);
        errorbar(PF,abs(Vr_time(t,:)),Vr_err_time(t,:),"-o","LineWidth",1.5)
%         errorbar(PF,Vr_time(t,:),Vr_err_time(t,:),"LineWidth",1.5)
        set(ax1,'NextPlot','add');
%         title('Vr','Color','black','FontWeight','bold')
        if t == NofCH/8
            legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
            %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
        end
        xlabel('PF充電電圧 [kV]')
        ylabel('|Vr| [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
            
        xlim([32.5 39.5]);
        ylim([-2 11]);
        
%         ax2 = subplot(3,1,2);
%         %     errorbar(PF,abs(Vz_time(t,:)),Vz_err_time(t,:),"LineWidth",1.5)
%         errorbar(PF,Vz_time(t,:),Vz_err_time(t,:),"LineWidth",1.5)
%         set(ax2,'NextPlot','add');
%         title('Vz','Color','black','FontWeight','bold')
%         if t == NofCH/8-1
%             legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
%             %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%         end
%         xlabel('PF充電電圧 [kV]')
%         ylabel('Vz[km/s]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([32.5 39.5]);
%         %     ylim([Vryl Vryh]);
%         
%         ax3 = subplot(3,1,3);
%         %     errorbar(PF,abs(Vr_time(t,:)),Vr_err_time(t,:),"LineWidth",1.5)
%         errorbar(PF,Va_time(t,:),Va_err_time(t,:),"LineWidth",1.5)
%         set(ax3,'NextPlot','add');
%         title('Vt','Color','black','FontWeight','bold')
%         if t == NofCH/8-1
%             legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
%             %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%         end
%         xlabel('PF充電電圧 [kV]')
%         ylabel('V\theta [km/s]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([32.5 39.5]);
%         ylim([Vtyl Vtyh]);
    end
elseif flag == 2%TFの場合
    for t = 1:NofCH/8%-1はRが大きい点を無視してる（光量が足りてなさそうで信頼性に欠けるため）
        
        ax1 = subplot(1,1,1);
%         errorbar(TF,abs(Vr_time(t,:)),Vr_err_time(t,:),"LineWidth",1.5)
        errorbar(TF,abs(Vr_time(t,:)),Vr_err_time(t,:),Vr_err_time(t,:),xneg,xpos,"-o","LineWidth",1.5)
%         errorbar(TF,Vr_time(t,:),Vr_err_time(t,:),"LineWidth",1.5)
        set(ax1,'NextPlot','add');
%         title('Vr','Color','black','FontWeight','bold')
        if t == NofCH/8
            legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
            %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
        end
%         xlabel('TF充電電圧 [kV]')
        xlabel('Bt [T]')
        ylabel('|Vr| [km/s]')
        ax = gca;
        ax.FontSize = 12;
        grid on
        
        R = corrcoef(TF,abs(Vr_time(t,:)))
        
        xlim([0.1875 0.212]);
        ylim([-2 11]);
        
%         ax2 = subplot(3,1,2);
%         %     errorbar(PF,abs(Vz_time(t,:)),Vz_err_time(t,:),"LineWidth",1.5)
%         errorbar(TF,Vz_time(t,:),Vz_err_time(t,:),"LineWidth",1.5)
%         set(ax2,'NextPlot','add');
%         title('Vz','Color','black','FontWeight','bold')
%         if t == NofCH/8-1
%             legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
%             %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%         end
%         xlabel('TF充電電圧 [kV]')
%         ylabel('Vz [km/s]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([3.35 4.05]);
%         %     ylim([Vryl Vryh]);
%         
%         ax3 = subplot(3,1,3);
%         %     errorbar(PF,abs(Vr_time(t,:)),Vr_err_time(t,:),"LineWidth",1.5)
%         e = errorbar(TF,Va_time(t,:),Va_err_time(t,:),"LineWidth",1.5);
%         e.LineStyle = '--';
% %         e.LineSpec = "--";
%         set(ax3,'NextPlot','add');
%         title('Vt','Color','black','FontWeight','bold')
%         if t == NofCH/8-1
%             legend('R = 10cm','R = 15cm','R = 20cm','R = 25cm','Location','best')
%             %             legend('TF = 4.0kV','TF = 3.8kV','TF = 3.6kV','TF = 3.4kV','Location','best')
%         end
%         xlabel('TF充電電圧 [kV]')
%         ylabel('V\theta [km/s]')
%         ax = gca;
%         ax.FontSize = 12;
%         grid on
%         
%         xlim([3.35 4.05]);
%         ylim([Vtyl Vtyh]);
    end
end