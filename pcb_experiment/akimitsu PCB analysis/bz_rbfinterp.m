function vq = bz_rbfinterp(rpos, zpos, grid2D, bz, ok, t)
%%入力
% bz[samplingnumber×ch],rpos(ch),zpos(ch),grid2D{rq(n×n),zq(n×n)},ok(128(true or false)),t[us]
%%出力
% vq(n×n) ：bzがn×nのグリッドにrbf補間されたもの

%%スムージングと関数の選択
smoothval=0.06;
func='multiquadric';%Gaussian, Linear, Cubic,multiquadric, Thinplate から選べる
const=0.05;%should be ~distance of the points
%const = 0.3;%125chの場合zqは0.0027m刻み,rqは0.0028m刻み
%%無視するチャンネルを除いたbzの散布データ（okのチャンネルのみ残す）
x=zpos(ok);
y=rpos(ok);
z=double(bz(t,ok))';

%%補間
vq= rbfinterp([grid2D.zq(:)'; grid2D.rq(:)'], rbfcreate([x' ; y'], z','RBFFunction',func,'RBFConstant',const,'RBFSmooth', smoothval ));%,'RBFSmooth', smoothval,'RBFConstant',const));
vq = reshape(vq, size(grid2D.zq));
end