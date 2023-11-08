%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%実験日(date)を含む実験ログを保存して
%実験日のログ位置(行番号)を取得
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exp_log,index,begin_row,end_row] = load_log(date)
%時短用(固定値で良い)
load_s = 5300;%実験ログ読み始め行番号(230309~)
load_f = 10000;%実験ログ読み終わり行番号(固定)
%実験ログがcdになければ最新の実験ログをwebから取得

%%実験中など毎回ログを更新したい時用%%%%%%%%%%%%%%%%%%%%
% save_log();
% disp('最新の実験ログ(exp_log.xlsx)を保存しました。');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) >= date
        disp(append(string(FileDate), '更新の実験ログを使用します。'));
    else
        save_log();
        disp('最新の実験ログ(exp_log.xlsx)を保存しました。');
    end
else
    save_log();
    disp('最新の実験ログ(exp_log.xlsx)を保存しました。');
end
%実験ログ中の実験日に対応する範囲を特定
exp_index = readcell('exp_log.xlsx','Sheet','log','Range','A1:BC1');

index.date = find(strncmp('date', exp_index,4));%列番号を取得
index.shot = find(strncmp('shot', exp_index,4));%列番号を取得
index.a039 = find(strncmp('a039', exp_index,4));%列番号を取得
index.a039 = index.a039(1,1);%列番号を取得
index.a039_TF = find(strncmp('a039_TF', exp_index,7));%列番号を取得
index.PF1 = find(strncmp('CB1[kV]', exp_index,7));%列番号を取得
index.PF2 = find(strncmp('CB2[kV]', exp_index,7));%列番号を取得
index.TF = find(strncmp('TF[kV]', exp_index,6));%列番号を取得
index.EF = find(strncmp('EF[A]', exp_index,5));%列番号を取得
index.IDS288ch_trg = find(strncmp('DopplerDelay[us]', exp_index,16));%列番号を取得
index.IDS288ch_exp_w = find(strncmp('DopplerWidth[us]', exp_index,16));%列番号を取得
index.IDSP_trg = find(strncmp('IDSPDelay[us]', exp_index,13));%列番号を取得
index.IDSP_exp_w = find(strncmp('IDSPWidth[us]', exp_index,13));%列番号を取得
index.IDSP_gain = find(strncmp('IDSPGain', exp_index,8));%列番号を取得
index.IDSP_minR = find(strncmp('IDSPminR[cm]', exp_index,9));%列番号を取得
index.IDSP_minZ = find(strncmp('IDSPminZ[cm]', exp_index,9));%列番号を取得

exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':BC' num2str(load_f)]);
[n_row,~] = size(exp_log);

begin_row = find(exp_log(:,3) == date,1,"first");%実験日の最初のshotの行番号を取得
if isempty(begin_row)
    error('実験日が実験ログ中に存在しません。')
end
end_row = begin_row;
while (end_row<n_row && not(isnan(exp_log(end_row+1,4))) && (isnan(exp_log(end_row+1,3)) || (exp_log(end_row+1,3) == date)))%日付がNaN or date&&shot番号が記入済=実験日のshot
    end_row = end_row+1;
end%実験日の最後のshotの行番号を取得
end
