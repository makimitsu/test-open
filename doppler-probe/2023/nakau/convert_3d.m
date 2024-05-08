NofCH = 32;%チャンネル数
NT = 4;%計測トリガ時刻数
NData = 4;%同一計測トリガでの計測ショット数
nz = 1;
r_measured = zeros(NofCH/8,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/8,nz);%ベクトルプロットz座標1列目data1、2列目data2
% z_measured(:,1) = 0;%data1のz座標
% r_measured(:,1) = [10, 15, 20, 25];%data1のr座標
% z_measured(:,1) = [-0.3 -0.3 -0.3 -0.3 0 0 0 0 0.3 0.3 0.3 0.3];%data1のz座標
% r_measured(:,1) = [10 15 20 25 10 15 20 25 10 15 20 25];%data1のr座標

for t = 1:NT
    Vk = zeros(NofCH/8,NData);
    Vl = zeros(NofCH/8,NData);
    Vr = zeros(NofCH/4,NData);
    Va = zeros(NofCH/8,NData);
    Vz = zeros(NofCH/8,NData);
    for ndata = 1:NData
        load(['magflow/mat/save_',num2str(462+4*t),'us_',num2str(ndata),'.mat'],'V','Ti')
        Vk(:,ndata) = V(1:2:NofCH/4,1);%座標をk,lに変換
        Vl(:,ndata) = -V(2:2:NofCH/4,1);
%         Vr(:,ndata) = (V(1:2:NofCH/4,2) + V(2:2:NofCH/4,2))/2 %Vrは同じ値になるはずだが平均をとっておく
        Vr(:,ndata) = V(1:NofCH/4,2);

        Va(:,ndata) = Vk(:,ndata) / sqrt(2) - Vl(:,ndata) / sqrt(2);%廊下側からみて下がθ+
        Vz(:,ndata) = -(Vk(:,ndata) / sqrt(2) + Vl(:,ndata) / sqrt(2));%シールドルームから奥側が+になる
        save(['magflow/mat/3d_',num2str(462+4*t),'us_',num2str(ndata),'.mat'],'Va','Vz','Vr','Ti')
    end
    
    a = 0 + t
end

% %再度座標をa,z,rに変換して3d直交座標で表示できるようにする
% for t = 1:NT
% 
%     for ndata = 1:NData
%     end
% end







