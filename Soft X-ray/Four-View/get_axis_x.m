function [magaxis,xpoint] = get_axis_x(grid2D,data2D,timeIndex)

psi = data2D.psi;
psi = psi(:,:,timeIndex);
rq = grid2D.rq;
zq = grid2D.zq;

% 各z座標について、r方向の最大値とそのインデックスを取得
% [M,I] = max(psi,[],1);
[M,I] = max(psi,[],1,'linear');
% 取得されたr方向最大値ベクトルからの極大値（軸位置候補）を2つ取得
TF = islocalmax(M,'MaxNumExtrema',2);

% ここから条件分岐
% 軸位置を2つ取得できた場合
if numel(I(TF)) == 2
    % 軸位置のrz座標を取得
    magaxis.r = rq(I(TF));
    magaxis.z = zq(I(TF));
    magaxis.psi = psi(I(TF));
    % 2つの軸の間でx点を検索
    M_x = M(find(TF,1):find(TF,1,'last'));
    I_x = I(find(TF,1):find(TF,1,'last'));
    TF_x = islocalmin(M_x,'MaxNumExtrema',1);
    if ~isempty(TF_x)
        xpoint.r = rq(I_x(TF_x));
        xpoint.z = zq(I_x(TF_x));
        xpoint.psi = psi(I_x(TF_x));
    else
        xpoint.r = NaN;
        xpoint.z = NaN;
        xpoint.psi = NaN;
    end
elseif numel(I(TF)) == 1
    % 軸位置のrz座標を取得
    magaxis.r = rq(I(TF));
    magaxis.z = zq(I(TF));
    magaxis.psi = psi(I(TF));
    xpoint.r = NaN;
    xpoint.z = NaN;
    xpoint.psi = NaN;
else
    magaxis.r = NaN;
    magaxis.z = NaN;
    magaxis.psi = NaN;
    xpoint.r = NaN;
    xpoint.z = NaN;
    xpoint.psi = NaN;
end


end