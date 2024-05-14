%%ST40データ取り込み練習

%mdsplusのパスを通す
setenv('MDSPLUS_DIR','/usr/local/mdsplus');
addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
%ST40に接続
mdsconnect('87.224.94.202:8000');
%shotにアクセス
mdsopen('ST40',11215);
%データ取得
d=mdsvalue("_d="+'.MAG.BEST.ROG.MCT:I');
% %データに対応する時刻を取得
% t = mdsvalue("dim_of(_d)");
% plot(t,d)

% mdsclose;
% mdsdisconnect;
