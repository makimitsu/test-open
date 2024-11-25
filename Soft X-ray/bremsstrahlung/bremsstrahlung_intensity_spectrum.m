clearvars -except date IDXlist doSave doFilter doNLR
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath '/Users/shohgookazaki/Documents/GitHub/test-open'/'Soft X-ray'/Four-View; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
% ~/Documents/MATLAB にてstartup.mを作って、その中でsetenv('パス名','アドレス')していくと自動になる。
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory_path=getenv('pre_processed_directory_path');


%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','Plot:'};
definput = {'','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('plot','var')
    definput{3} = num2str(plot);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
doplot = logical(str2num(cell2mat(answer(3))));


%-----------スプレッドシートからデータ抜き取り--------------------%
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

T=searchlog(T,'date',date);
if isnan(T.shot(1))
    T(1, :) = [];
end
n_data=numel(IDXlist);%計測データ
shotlist = [T.a039(IDXlist), T.a040(IDXlist)];
tfshotlist = [T.a039_TF(IDXlist), T.a040_TF(IDXlist)];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1); % 39が計測データ数だけ縦に並ぶ。
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);


% Determine the total number of time steps across all shots
n_times_per_shot = 8; % since times = start:interval:(start+interval*7)
total_rows = n_data * n_times_per_shot;
transmission_n = 10;

transmissionfile = sprintf('transmission_data%d.mat',transmission_n);
if ~isfile(transmissionfile)
    transmission(transmission_n);
end
load(transmissionfile, 'transmission_matrix','means', 'U','s','v','M','K');

% Preallocate space for I_data and x_data
I_data = zeros(total_rows, transmission_n); % Preallocate for I values
x_data = zeros(total_rows, transmission_n);             % Preallocate for time values
t_data = zeros(total_rows, 1);
shot_data = zeros(total_rows, 1);          % For storing shotlist values


row_idx = 1;

PCB.trange=400:800;%【input】計算時間範囲
PCB.n=50; %【input】rz方向のメッシュ数
PCB.restart = 0;

for i = 1:n_data
    shot = IDXlist(i);
    start = startlist(i);
    interval = intervallist(i);

    PCB.idx = IDXlist(i);
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.shot == PCB.tfshot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    PCB.date = date;
    TF=TFlist(i);
    [grid2D,data2D] = process_PCBdata_200ch(PCB,pathname); %process_PCBdata_200ch.mに行く
    [magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D); %時間ごとの磁気軸、X点を検索

    options = 'LF_TP';
    
    times = start:interval:(start+interval*7);

    dirPath = getenv('SXR_MATRIX_DIR');
    matrixFolder = strcat(dirPath,'/',options,'/',num2str(date),'/shot',num2str(shot));


    for t = times
        disp(strcat(num2str(shot),num2str(t)));
        number = (t-start)/interval+1;
        matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
        
        load(matrixPath,'EE1','EE2','EE3','EE4');
        % EE  = cat(3,EE1,EE2,EE3,EE4);
        if t == times(1)
            
            n = size(EE1, 1);

            %-----------下流領域よう-------------------
            % 長方形の開始位置とサイズ
            % rect_height = floor(n / 5);       % 長方形の高さ
            % rect_width = floor(4 * n / 5);    % 長方形の幅
            % 
            % % 長方形の位置 (行方向は下部、列方向は中央)
            % start_row = n - rect_height + 1;         % 行の開始位置
            % end_row = n;                             % 行の終了位置
            % start_col = floor((n - rect_width) / 2) + 1;  % 列の開始位置
            % end_col = start_col + rect_width - 1;    % 列の終了位置
            % 
            % % 長方形領域の抽出
            % rect_region = EE1(start_row:end_row, start_col:end_col);
            % disp(size(rect_region));
            % disp(n)
            % % 平均値の計算
            % rect_mean = mean(rect_region(:));

            %--------------------X点用----------------------------------------------
            %5分割した時の中心がX点だと仮定した場合の計算。これを初期Epにする
            % 中心部分の開始位置を計算
            startIndex = round(n/2 - n/6) + 1; % n/10を引き、切り上げて中心を確保
            endIndex = startIndex + round(n/3) - 1;
    
            % 中心部分を抽出
            centerBlock1 = EE1(startIndex:endIndex, startIndex:endIndex);
            centerBlock2 = EE2(startIndex:endIndex, startIndex:endIndex);
            centerBlock3 = EE3(startIndex:endIndex, startIndex:endIndex);
            centerBlock4 = EE4(startIndex:endIndex, startIndex:endIndex);
            % 平均を計算
            E1 = mean(centerBlock1, 'all');
            E2 = mean(centerBlock2, 'all');
            E3 = mean(centerBlock3, 'all');
            E4 = mean(centerBlock4, 'all');
            


            Ep_previous = cat(1,E1,E2,E3,E4);
            % Ep_previous = cat(1,E1,E2,E3);
        end 
        Ep = xpointarea(EE1, EE2, EE3, EE4, xPointList,t,data2D, Ep_previous);
        Ep_previous = Ep;
        % I = transmission_matrix\Ep;
        [I,conv] = get_distribution(M,K,U,s,v,Ep, transmission_matrix, means);
        if conv
            % Store I, x, and t values in preallocated arrays
            I_data(row_idx, :) = I'; % Store I as a row
            x_data(row_idx, :) = means;     % store means
            t_data(row_idx) = t;     % Store t in t_data as well
            shot_data(row_idx) = shot; % Store shot in shot_data
            row_idx = row_idx + 1;   % Increment row index
        end    
        
    end
end

figure; % 新しい図を作成

% データセットが複数ある場合、例えばIとmeansがそれぞれの列にデータを持つと仮定
for i = 1:size(I_data, 1)
    plot(x_data(i, :).', I_data(i, :), 'o-'); % 各データセットをプロット
    hold on;
end
hold off;
xlabel('Energy[eV]'); % x軸のラベル
ylabel('Intensity[a.u.]'); % y軸のラベル
title('Bremsstrahlung Intensity'); % タイトル
grid on; % グリッドを表示
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% xlim([10 300]);
% ylim([1e-4 1e2]);

% Concatenate x_data, t_data, and I_data for export
output_data = [shot_data, t_data, x_data, I_data];
writematrix(output_data, 'xpoint_data.xlsx', 'Sheet', 1, 'Range', 'A1');