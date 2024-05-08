function [center] = load_calibration(date,ICCD)

switch ICCD.line
    case 'Ar'%アルゴンの時
        if date < 230303
            warning('Input error in date.')%ICCD.lineの入力エラー
        elseif date < 230521
            center_file = '230303_Xe4834_calibation.txt';%中心データファイル名
        elseif 230521 <= date
            center_file = '230521_Xe4834_calibation.txt';%中心データファイル名
        else
            warning('Input error in date.')%ICCD.lineの入力エラー
        end
    case 'H'%水素の時
        % center_file = 'Hbeta_calibration.txt';%中心データファイル名
        warning('Sorry, not ready for H experiment.')%ICCD.lineの入力エラー
        return;
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        return;
end
center = importdata(center_file);%校正中心座標を取得