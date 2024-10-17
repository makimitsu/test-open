function plot_psi280ch(PCB,pathname)
shot = PCB.shot;
trange = PCB.trange;
start = PCB.start;
% filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% if exist(filename,"file")==0
%     disp('No rawdata file -- Start generating!')
%     rawdataPath = pathname.rawdata;
%     save_dtacq_data(dtacq_num, shot, tfshot,rawdataPath)
%     % return
% end
% load(filename,'rawdata');%1000×192
% 
% %正しくデータ取得できていない場合はreturn
% if numel(rawdata)< 500
%     return
% end

% filename = strcat(pathname.pre_processed_directory,'/a039_',num2str(shot(1)),'.mat');
% if exist(filename,'file') == 0
%     doCalculation = true;
%     disp('no processed data -- start calculation');
% else
%     doCalculation = false;
%     disp('loading processed data');
% end
% 
% 
% if doCalculation
%     %較正係数のバージョンを日付で判別
%     sheets = sheetnames('coeff200ch.xlsx');
%     sheets = str2double(sheets);
%     sheet_date=max(sheets(sheets<=date));
%     C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
%     r_shift = 0.00;
%     ok = logical(C(:,14));
%     dtacq_num_list = C(:,1);
%     dtaq_ch = C(:,2);
%     polarity=C(:,13);
%     coeff=C(:,12);
%     zpos=C(:,9);
%     rpos=C(:,10)+r_shift;
%     ch=C(:,7);
% 
%     if ismember(39,dtacq_num_list)
%         filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
%         if exist(filename1,"file")==0
%             disp('No rawdata file of a039 -- Start generating!')
%             rawdataPath = pathname.rawdata;
%             save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
%             % disp(['File:',filename1,' does not exit']);
%             % return
%         end
%         a039_raw = importdata(filename1);
%     end
%     if ismember(40,dtacq_num_list)
%         filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
%         if exist(filename2,"file")==0
%             disp('No rawdata file of a040 -- Start generating!')
%             rawdataPath = pathname.rawdata;
%             save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
%             % disp(['File:',filename2,' does not exit']);
%             % return
%         end
%         a040_raw = importdata(filename2);
%     end
% 
%     raw = zeros(1000,length(dtaq_ch));
%     for i = 1:length(dtaq_ch)
%         if dtacq_num_list(i) == 39
%             raw(:,i) = a039_raw(:,dtaq_ch(i));
%         elseif dtacq_num_list(i) == 40
%             raw(:,i) = a040_raw(:,dtaq_ch(i));
%         end
%     end
% 
%     b=raw.*coeff';%較正係数RC/NS
%     b=b.*polarity';%極性揃え
% 
%     %デジタイザchからプローブ通し番号順への変換
%     bz=zeros(1000,100);
%     bt=bz;
%     ok_bz=false(100,1);
%     ok_bt=ok_bz;
%     zpos_bz=zeros(100,1);
%     rpos_bz=zpos_bz;
%     zpos_bt=zpos_bz;
%     rpos_bt=zpos_bz;
% 
%     %digital filter
%     windowSize = 8;
%     bb = (1/windowSize)*ones(1,windowSize);
%     aa = 1;
% 
%     for i=1:length(ch)
%         b(:,i) = filter(bb,aa,b(:,i));
%         b(:,i) = b(:,i) - mean(b(1:40,i));
%         if rem(ch(i),2)==1
%             bz(:,ceil(ch(i)/2))=b(:,i);
%             ok_bz(ceil(ch(i)/2))=ok(i);
%             zpos_bz(ceil(ch(i)/2))=zpos(i);
%             rpos_bz(ceil(ch(i)/2))=rpos(i);
%         elseif rem(ch(i),2)==0
%             bt(:,ch(i)/2)=b(:,i);
%             ok_bt(ceil(ch(i)/2))=ok(i);
%             zpos_bt(ceil(ch(i)/2))=zpos(i);
%             rpos_bt(ceil(ch(i)/2))=rpos(i);
%         end
%     end
% 
%     % zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
%     zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
%     rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
%     rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
%     [zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
%     % [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
%     [zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
%     ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
%     ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
%     for i = 1:length(ok_bt)
%         if rpos_bt(i) > (r_shift)
%             index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
%             ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
%         end
%         index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
%         ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
%     end
% 
%     grid2D=struct(...
%         'zq',zq,...
%         'rq',rq,...
%         'zprobepcb',zprobepcb,...
%         'rprobepcb',rprobepcb,...
%         'rprobepcb_t',rprobepcb_t,...1
%         'ok_bz_matrix',ok_bz_matrix,...
%         'ok_bt_matrix',ok_bt_matrix);
%     grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);
% 
%     clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix
% 
%     % probecheck_script;
% 
%     %data2Dcalc.m
%     r_EF   = 0.5 ;
%     n_EF   = 234. ;
% 
%     if date<221119
%         z1_EF   = 0.68;
%         z2_EF   = -0.68;
%     else
%         z1_EF   = 0.78;
%         z2_EF   = -0.78;
%     end
%     [Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
%     clear EF r_EF n_EF i_EF z_EF
% 
%     data2D=struct(...
%         'psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%         'trange',trange);
% 
%     % ******************* no angle correction ********************
%     for i=1:size(trange,2)
%         t=trange(i);
% 
%         %Bzの二次元補間(線形fit)
%         vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
%         B_z = -Bz_EF+vq;
%         B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
% 
%         % PSI計算
%         data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
%         % data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);
%         % このままだと1/2πrが計算されてないので
%         [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
%         data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
%         data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
%         data2D.Bt(:,:,i)=B_t;
%         data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
%     end
% else
%     load(filename,'data2D','grid2D');
% end

[grid2D,data2D] = process_PCBdata_280ch(PCB,pathname);
% [grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);

% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D); %時間ごとの磁気軸、X点を検索

% プロット部分
figure('Position', [0 0 1500 1500],'visible','on');

% figure('Position', [0 0 1500 1500],'visible','on');
% start=30;
dt = PCB.dt;
%  t_start=470+start;
Et_t = zeros(1,16);
times = start+dt:dt:start+dt*16;
 for m=1:16 %図示する時間
     i=start+m.*dt; %end
     t=trange(i);
     subplot(4,4,m)
     switch(PCB.type)
         case 1
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');clim([-5e-3,5e-3])%psi
         case 2
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Bz
         case 3
             % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),40,'LineStyle','none');clim([0,0.3])%Bt
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
         case 4
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Br
         case 5
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none');clim([-500,400])%Et
             % 時刻iにおけるx点の[r,z]座標（インデックス）を取得
             idxR = knnsearch(grid2D.rq(:,1),xPointList.r(i));
             idxZ = knnsearch(grid2D.zq(1,:).',xPointList.z(i));
             % ±3程度の範囲でEtの平均を計算
             Et_t(m) = mean(data2D.Et(max(1,idxR-3):min(idxR+3,numel(grid2D.rq(:,1))),max(1,idxZ-3):min(numel(grid2D.zq(1,:)),idxZ+3),i),'all');
         case 6
             contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');clim([-0.8*1e+6,0]);%clim([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
     end

     % ポロイダル磁場ベースでのプロット
     % Bp = sqrt((data2D.Br(:,:,i)).^2+(data2D.Bz(:,:,i)).^2);
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),Bp,30,'LineStyle','none');clim([0,0.02]);
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),abs(data2D.Bt(:,:,i)./data2D.Bz(:,:,i)),30,'LineStyle','none');clim([-1,1]);
     % sepa = zeros(size(data2D.psi(:,:,i)));
     % x_psi = xPointList.psi(i);
     % sepa(((data2D.psi(:,:,i)<=x_psi+5e-4&data2D.psi(:,:,i)>=x_psi-5e-5)|(data2D.psi(:,:,i)<=x_psi*1.1&data2D.psi(:,:,i)>=x_psi*0.9)))=1;
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),sepa,'LineStyle','none');
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');clim([-5e-3,5e-3])%psi
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Bz
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Br
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),40,'LineStyle','none');clim([0,0.3])%Bt
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');clim([-0.8*1e+6,0]);%clim([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none');clim([-500,400])%Et
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Bt_th(:,:,i),20,'LineStyle','none');clim([0,0.3])%Bt
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
    % clim([-0.1,0.1])%Bz
     % clim([0,0.3])%Bt
    % clim([-5e-3,5e-3])%psi
    % clim([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
     % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置

    % timeIndex = find(trange==t);
    % [magaxis,xpoint] = get_axis_x(grid2D,data2D,t);
    % plot(magaxis.z,magaxis.r,'ro');
    % plot(xpoint.z,xpoint.r,'rx');
    plot(magAxisList.z(:,i),magAxisList.r(:,i),'ko');
    plot(xPointList.z(i),xPointList.r(i),'kx');

    hold off
    title(string(t)+' us')
%     xlabel('z [m]')
%     ylabel('r [m]')
 end

 sgtitle(strcat('shot',num2str(shot)));

 if PCB.type==5
     figure;plot(times+400,Et_t,'LineWidth',3);
     ax = gca;
     ax.FontSize = 18;
     xlabel('time [us]');ylabel('Electric field [V/m]');
 end
 % drawnow

end