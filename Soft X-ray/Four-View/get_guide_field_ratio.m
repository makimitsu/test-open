function [B_r,B_t,b] = get_guide_field_ratio(PCB,pathname)

% magDataDir = pathname.MAGDATA;
% magDataFile = strcat(magDataDir,'/',num2str(PCB.data),'.mat');
% if exist(magDataFile,'file')
%     load(magDataFile,'idxList','BrList','BtList','bList');
% end

trange = PCB.trange;
% newTimeRange = find(trange>=440&trange<=500);
newTimeRange = 1:numel(trange);
[grid2D,data2D] = process_PCBdata_280ch(PCB,pathname);
Br = data2D.Br;
rq = grid2D.rq;
zq = grid2D.zq;

% PCB_ref = PCB;
% PCB_ref.shot = PCB.tfshot;
% PCB_ref.tfshot = [0,0];
% [~,data2D_ref] = process_PCBdata_280ch(PCB_ref,pathname);
% Bt_ref = data2D_ref.Bt;

trange = trange(newTimeRange);
Br = Br(:,:,newTimeRange);

B_reconnection = zeros(1,numel(trange));
B_guide = zeros(1,numel(trange));
mergingRatio = zeros(1,numel(trange));

[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D);
[I_TF,x,aquisition_rate] = get_TF_current(PCB,pathname);
% timing = x/aquisition_rate==time;
m0 = 4*pi*10^(-7);
% r = xpoint.r;
% Bt = m0*I_TF(timing)*1e3*12/(2*pi()*r);

% 20240108 できるだけループを使わずに一括で処理したい
% 特にget_axis_xを呼び出しすぎてる
for i = 1:numel(trange)
    time = trange(i);
    % [magaxis,xpoint] = get_axis_x(grid2D,data2D,time);
    magaxis.r = magAxisList.r(:,i);
    magaxis.z = magAxisList.z(:,i);
    % if numel(magaxis.r) == 2
    if magaxis.z(1)~=magaxis.z(2) && ~isnan(magaxis.r(1))
        range = rq>=min(magaxis.r)&rq<=max(magaxis.r)&zq>=min(magaxis.z)&zq<=max(magaxis.z);
        % 値の取り方は考えた方がいい、Btは理論値でもよさそう
        Br_t = Br(:,:,i);
        B_reconnection(1,i) = max(abs(Br_t(range)),[],'all');
        % diffusionRegion = (abs(rq-xpoint.r)+abs(zq-xpoint.z))<=0.02;
        % Bt_t = Bt_ref(:,:,i);
        % B_guide(1,i) = mean(Bt_t(diffusionRegion));
        % B_guide(1,i) = get_B_troidal(PCB,grid2D,data2D,pathname,time);
        timing = x/aquisition_rate==time;
        B_guide(1,i) = m0*I_TF(timing)*1e3*12/(2*pi()*xPointList.r(1,i));
        mergingRatio(1,i) = xPointList.psi(i)/mean(magAxisList.psi(:,i));
    else
        B_reconnection(1,i) = NaN;
        B_guide(1,i) = NaN;
        mergingRatio(1,i) = NaN;
    end
end

nanMask = isnan(B_reconnection);
nanMask(2:end-1) = nanMask(1:end-2)&nanMask(3:end);
B_reconnection(nanMask) = NaN;
B_guide(nanMask) = NaN;
mergingRatio(nanMask) = NaN;

% figure;
% subplot(1,3,1);plot(trange,mergingRatio);xlabel('time');ylabel('merging ratio');
% subplot(1,3,2);plot(trange,B_reconnection);xlabel('time');ylabel('B_r');
% subplot(1,3,3);plot(trange,B_guide);xlabel('time');ylabel('B_g');


timing = knnsearch(mergingRatio.',0.1);
% timing = find(mergingRatio==0,1,'last');
B_r = B_reconnection(1,timing);
B_t = B_guide(1,timing);
b = B_t/B_r;

% if exist(magDataFile,'file')
%     idxList(end+1) = PCB.idx;
%     BrList(end+1) = B_r;
%     BtList(end+1) = B_t;
%     bList(end+1) = b;
% else
%     idxList = PCB.idx;
%     BrList = B_r;
%     BtList = B_t;
%     bList = b;
% end
% 
% save(magDataFile,'idxList','BrList','BtList','bList');

end