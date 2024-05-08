function [center] = load_IDSPcalib(IDSP)

switch char(IDSP.line)
    case 'Ar'%アルゴンの時
        sheets = sheetnames('coeffIDSP_Ar.xlsx');
        sheets = str2double(sheets);
        sheet_date=max(sheets(sheets<=IDSP.date));
        center = readmatrix('coeffIDSP_Ar.xlsx','Sheet',num2str(sheet_date));
    case 'H'%水素の時
        % center_file = 'Hbeta_calibration.txt';%中心データファイル名
        warning('Sorry, not ready for H experiment.')%ICCD.lineの入力エラー
        return;
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        return;
end
