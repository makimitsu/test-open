%%logから日付や語句検索でショットを探す関数
% node='テーブルのヘッダー'：　検索するもの　'operator','date',など 
% pat=数値 or str :　抽出するパターン
function output=searchlog(inputtable,node,pat)
%%名前,commentで検索(~を含むの検索)
if ischar(pat)
rows = contains( inputtable.(node),pat,'IgnoreCase',true); 
%%日付,shotなどで検索(一致検索)
else
 rows = inputtable.(node) == pat ;
end
output=inputtable(rows,:);
end