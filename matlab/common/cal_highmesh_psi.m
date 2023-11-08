function [highmesh_PCBgrid2D,highmesh_PCBdata2D] = cal_highmesh_psi(PCB,pathname)

filename = strcat(pathname.processeddata,'/mesh',num2str(PCB.mesh),'_a039_',num2str(PCB.shot(1)),'.mat');
if exist(filename,'file') == 0
    doCalculation = true;
else
    doCalculation = false;
end

if doCalculation
    %較正係数のバージョンを日付で判別
    sheets = sheetnames('coeff200ch.xlsx');
    sheets = str2double(sheets);
    sheet_date=max(sheets(sheets<=PCB.date));
    C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
    r_shift = 0.00;
    ok = logical(C(:,14));
    dtacq_num_list = C(:,1);
    dtaq_ch = C(:,2);
    polarity=C(:,13);
    coeff=C(:,12);
    zpos=C(:,9);
    rpos=C(:,10)+r_shift;
    ch=C(:,7);

    if ismember(39,dtacq_num_list)
        filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(PCB.shot(1)),'_tfshot',num2str(PCB.tfshot(1)),'.mat');
        filename3 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(PCB.shot(1)),'_tfshot0.mat');
        if exist(filename1,"file")==0
            disp('No rawdata file of a039 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(39, PCB.shot(1), PCB.tfshot(1),rawdataPath)
            % disp(['File:',filename1,' does not exit']);
            % return
        end
        if exist(filename3,"file")==0
            disp('No rawdata0 file of a039 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(39, PCB.shot(1), 0,rawdataPath)
            % disp(['File:',filename1,' does not exit']);
            % return
        end
        a039_raw = importdata(filename1);
        a039_raw0 = importdata(filename3);
    end
    if ismember(40,dtacq_num_list)
        filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(PCB.shot(2)),'_tfshot',num2str(PCB.tfshot(2)),'.mat');
        filename4 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(PCB.shot(2)),'_tfshot0.mat');
        if exist(filename2,"file")==0
            disp('No rawdata file of a040 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(40, PCB.shot(2), PCB.tfshot(2),rawdataPath)
            % disp(['File:',filename2,' does not exit']);
            % return
        end
        if exist(filename4,"file")==0
            disp('No rawdata0 file of a040 -- Start generating!')
            rawdataPath = pathname.rawdata;
            save_dtacq_data(40, PCB.shot(2), 0,rawdataPath)
            % disp(['File:',filename1,' does not exit']);
            % return
        end
        a040_raw = importdata(filename2);
        a040_raw0 = importdata(filename4);
    end

    raw = zeros(1000,length(dtaq_ch));
    raw0 = zeros(1000,length(dtaq_ch));
    for i = 1:length(dtaq_ch)
        if dtacq_num_list(i) == 39
            raw(:,i) = a039_raw(:,dtaq_ch(i));
            raw0(:,i) = a039_raw0(:,dtaq_ch(i));
        elseif dtacq_num_list(i) == 40
            raw(:,i) = a040_raw(:,dtaq_ch(i));
            raw0(:,i) = a040_raw0(:,dtaq_ch(i));
        end
    end

    b=raw.*coeff';%較正係数RC/NS
    b=b.*polarity';%極性揃え
    b0=raw0.*coeff';%較正係数RC/NS
    b0=b0.*polarity';%極性揃え

    %デジタイザchからプローブ通し番号順への変換
    bz=zeros(1000,100);
    bt=bz;
    ok_bz=false(100,1);
    ok_bt=ok_bz;
    zpos_bz=zeros(100,1);
    rpos_bz=zpos_bz;
    zpos_bt=zpos_bz;
    rpos_bt=zpos_bz;

    %digital filter
    windowSize = 8;
    bb = (1/windowSize)*ones(1,windowSize);
    aa = 1;

    for i=1:length(ch)
        b(:,i) = filter(bb,aa,b(:,i));
        b(:,i) = b(:,i) - mean(b(1:40,i));
        b0(:,i) = filter(bb,aa,b0(:,i));
        b0(:,i) = b0(:,i) - mean(b0(1:40,i));
        if rem(ch(i),2)==1
            bz(:,ceil(ch(i)/2))=b(:,i);
            ok_bz(ceil(ch(i)/2))=ok(i);
            zpos_bz(ceil(ch(i)/2))=zpos(i);
            rpos_bz(ceil(ch(i)/2))=rpos(i);
        elseif rem(ch(i),2)==0
            bt(:,ch(i)/2)=b(:,i);
            bt0(:,ch(i)/2)=b0(:,i);
            ok_bt(ceil(ch(i)/2))=ok(i);
            zpos_bt(ceil(ch(i)/2))=zpos(i);
            rpos_bt(ceil(ch(i)/2))=rpos(i);
        end
    end

    % zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
    zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
    rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
    rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
    [zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),PCB.mesh),linspace(min(rpos_bz),max(rpos_bz),PCB.mesh));
    % [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
    [zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
    ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
    ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
    for i = 1:length(ok_bt)
        if rpos_bt(i) > (r_shift)
            index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
            ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
        end
        index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
        ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
    end

    highmesh_PCBgrid2D=struct(...
        'zq',zq,...
        'rq',rq,...
        'zprobepcb',zprobepcb,...
        'rprobepcb',rprobepcb,...
        'rprobepcb_t',rprobepcb_t,...
        'ok_bz_matrix',ok_bz_matrix,...
        'ok_bt_matrix',ok_bt_matrix);
    grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);

    clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix

    % probecheck_script;

    %data2Dcalc.m
    r_EF   = 0.5 ;
    n_EF   = 234. ;

    if PCB.date<221119
        z1_EF   = 0.68;
        z2_EF   = -0.68;
    else
        z1_EF   = 0.78;
        z2_EF   = -0.78;
    end
    [Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,PCB.i_EF,n_EF,highmesh_PCBgrid2D.rq,highmesh_PCBgrid2D.zq,false);
    clear EF r_EF n_EF PCB.i_EF z_EF

    highmesh_PCBdata2D=struct(...
        'psi',zeros(size(highmesh_PCBgrid2D.rq,1),size(highmesh_PCBgrid2D.rq,2),size(PCB.trange,2)),...
        'Bt',zeros(size(highmesh_PCBgrid2D.rq,1),size(highmesh_PCBgrid2D.rq,2),size(PCB.trange,2)),...
        'Bt_ext',zeros(size(highmesh_PCBgrid2D.rq,1),size(highmesh_PCBgrid2D.rq,2),size(PCB.trange,2)),...
        'trange',PCB.trange);

    % ******************* no angle correction ********************
    for i=1:size(PCB.trange,2)
        t=PCB.trange(i);
        %Bzの二次元補間(線形fit)
        vq = bz_rbfinterp(rpos_bz, zpos_bz, highmesh_PCBgrid2D, bz, ok_bz, t);
        B_z = -Bz_EF+vq;
        B_t = bz_rbfinterp(rpos_bt, zpos_bt, highmesh_PCBgrid2D, bt0, ok_bt, t);
        % PSI計算
        highmesh_PCBdata2D.psi(:,:,i) = cumtrapz(highmesh_PCBgrid2D.rq(:,1),2*pi*B_z.*highmesh_PCBgrid2D.rq(:,1),1);
        highmesh_PCBdata2D.Bt(:,:,i)=B_t;
    end
    highmesh_PCBdata2D.Bt_ext = cal_Bt_theory(PCB,highmesh_PCBgrid2D);
else
    load(filename,'highmesh_PCBdata2D','highmesh_PCBgrid2D');
end
% ***********************************************

if isstruct(highmesh_PCBgrid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

if doCalculation
    clearvars -except highmesh_PCBdata2D highmesh_PCBgrid2D PCB pathname filename;
    save(filename)
end

end

