%%%%%%%%%%%%%%%%%
% pcbプローブの角度計算
% TF only のshotを使って計算する
%%%%%%%%%%%%%%%%%

clear all
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所


shot = [2417;896]; % TF only のshot番号
tfshot = [0;0];
date = 230830;%【input】計測日
figure_switch = [1 1 0 0]; %bz1,bz2,bt1,bt2

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date = max(sheets(sheets <= date)); % 計測日以前で最新バージョンの較正係数を使用
C = readmatrix('coeff200ch.xlsx', 'Sheet', num2str(sheet_date)); % 較正係数の読み込み
r_shift = 0.00; % 【input】プローブの差し込み具合を変更した場合は記入
ok = logical(C(:,14)); % chが生きていれば1，死んでいれば0
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
polarity=C(:,13); % 極性
coeff=C(:,12); % 較正係数 RC/NS
zpos=C(:,9); % z位置[m]
rpos=C(:,10)+r_shift; % r位置[m]
ch=C(:,7); % デジタイザch番号


% if ismember(38,dtacq_num_list)
%     filename38 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(38),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
%     if exist(filename38,"file")==0
%         disp(['File:',filename38,' does not exit']);
%         return
%     end
%     a039_raw = importdata(filename39);
% end
if ismember(39,dtacq_num_list)
    filename39 = strcat(pathname.rawdata,'rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    if exist(filename39,"file")==0
        disp(['File:',filename39,' does not exit']);
        return
    end
    a039_raw = importdata(filename39);
end
if ismember(40,dtacq_num_list)
    filename40 = strcat(pathname.rawdata,'rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    if exist(filename40,"file")==0
        disp(['File:',filename40,' does not exit']);
        return
    end
    a040_raw = importdata(filename40);
end

raw = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 39
        raw(:,i) = a039_raw(:,dtaq_ch(i));
    elseif dtacq_num_list(i) == 40
        raw(:,i) = a040_raw(:,dtaq_ch(i));
    end
end

b=raw.*coeff';%較正係数RC/NS
b=b.*polarity';%極性揃え

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

%digital filter 移動平均フィルター（ノイズを含む信号の平滑化）
windowSize = 3;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:length(ch)
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

rad_lis = zeros(1,14);
zpos_lis = rad_lis;
idx_lis = [3:12,1,2,13,14];
for i=1:14
    idx_st = (i-1)*10 + 1;
    idx_ed = i*10;
    bz_max = maxabs(bz(300:700,idx_st:idx_ed));
    bt_max = maxabs(bt(300:700,idx_st:idx_ed));
    ok_bz_p = ok_bz(idx_st:idx_ed).';
    ok_bt_p = ok_bt(idx_st:idx_ed).';
    ok_p = logical(ok_bz_p.*ok_bt_p);
    % disp(bz_max)
    % disp(bt_max)
    % disp(ok_bz_p)
    % disp(ok_bt_p)
    % disp(ok_p)
    tan = bz_max./bt_max.*ok_p;
    rad = atan(tan);
    deg = rad*180/pi;
    rad_mean = mean(rad(ok_p));
    % disp(rad)
    % disp(rad_mean)
    p = idx_lis(i);
    rad_lis(p) = rad_mean;
    zpos_lis(p) = zpos_bt(idx_st);
end







%生信号描画用パラメータ
r = 7;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
bz_lim = [-0.1 0.1]; % bzのylim
bt_lim = [-0.4 0.4]; %[-0.1 0.4]; % btのylim
t_start = 100; % 横軸プロット領域（開始時間）
t_end =1000; % 横軸プロット領域（終了時間）


if figure_switch(1)
    f1=figure;
    f1.WindowState = 'maximized';
    for i=1:r
        for j=1:col
            subplot(r,col,(i-1)*col+j)
            if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i-1)+j)-1));
            xticks([t_start t_end]);
            xlim([t_start t_end]);
            ylim(bz_lim);
        end
    end
    sgtitle('Bz signal probe1-7')
end

if figure_switch(2)
    f2=figure;
    f2.WindowState = 'maximized';
    for i=1:r
        for j=1:col
            subplot(r,col,(i-1)*col+j)
            if ok_bz(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i+r-1)+j)-1));
            xticks([t_start t_end]);
            xlim([t_start t_end]);
            ylim(bz_lim);
        end
    end
    sgtitle('Bz signal probe8-14')
end

if figure_switch(3)
    f3=figure;
    f3.WindowState = 'maximized';
    for i=1:r
        for j=1:col
            subplot(r,col,(i-1)*col+j)
            if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i-1)+j)));
            xticks([t_start t_end]);
            %ylim([-0.2 0.2]);
            xlim([t_start t_end]);
            ylim(bt_lim);
        end
    end
    sgtitle('Bt signal probe1-7')
end

if figure_switch(4)
    f4=figure;
    f4.WindowState = 'maximized';
    for i=1:r
        for j=1:col
            subplot(r,col,(i-1)*col+j)
            if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i+r-1)+j)));
            xticks([t_start t_end]);
            %ylim([-0.2 0.2]);
            xlim([t_start t_end]);
            ylim(bt_lim);
        end
    end
    sgtitle('Bt signal probe8-14')
end

%disp(zpos_lis)
disp(rad_lis)
