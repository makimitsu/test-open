function [bz, ok_bz, ok_bz_plot] = ng_replace(bz, ok_bz, sheet_date)
%bz生信号を前後の信号で線形補完して置き換え
%   日付ごとに置き換えるverを変える
%   ok_bz:【input】生きているch→【output】生きているch＋補間したch
%   ok_bz_plot：【output】もともと生きているchのみ
if sheet_date<230103
    bz=bz;
    ok_bz=ok_bz;
elseif sheet_date<230110
    r1=0.021;
    r2=0.0535;
    bz(:,8)=(bz(:,7)+bz(:,9))./2;
    bz(:,16)=(bz(:,15)+bz(:,17))./2;
    bz(:,22)=(bz(:,21)+bz(:,23))./2;
    bz(:,26)=(bz(:,25)+bz(:,27))./2;
    bz(:,28)=(bz(:,27)+bz(:,29))./2;
    bz(:,31)=(r1.*bz(:,21)+r2.*bz(:,41))./(r1+r2);
    bz(:,32)=(r1.*bz(:,22)+r2.*bz(:,42))./(r1+r2);
    bz(:,33)=(r1.*bz(:,23)+r2.*bz(:,43))./(r1+r2);
    bz(:,39)=(bz(:,38)+bz(:,40))./2;
    bz(:,63)=(bz(:,62)+bz(:,64))./2;
    bz(:,95)=(bz(:,94)+bz(:,96))./2;
    ok_bz_plot=ok_bz;
    ok_bz([8 16 22 26 28 31 32 33 39 63 95])=true;
elseif sheet_date<230119
    r1=0.021;
    r2=0.0535;
    bz(:,16)=(bz(:,15)+bz(:,17))./2;
    bz(:,22)=(bz(:,21)+bz(:,23))./2;
    bz(:,26)=(bz(:,25)+bz(:,27))./2;
    bz(:,31)=(r1.*bz(:,21)+r2.*bz(:,41))./(r1+r2);
    bz(:,32)=(r1.*bz(:,22)+r2.*bz(:,42))./(r1+r2);
    bz(:,33)=(r1.*bz(:,23)+r2.*bz(:,43))./(r1+r2);
    bz(:,35)=(bz(:,34)+bz(:,36))./2;
    bz(:,39)=(bz(:,38)+bz(:,40))./2;
    bz(:,63)=(bz(:,62)+bz(:,64))./2;
    bz(:,95)=(bz(:,94)+bz(:,96))./2;
    ok_bz_plot=ok_bz;
    ok_bz([16 22 26 31 32 33 35 39 63 95])=true;
elseif sheet_date<230127
    r1=0.021;%0.1275
    r2=0.0535;%0.
    bz(:,22)=(bz(:,21)+bz(:,23))./2;
    bz(:,26)=(bz(:,25)+bz(:,27))./2;
    bz(:,28)=(bz(:,27)+bz(:,29))./2;
    bz(:,31)=(r1.*bz(:,21)+r2.*bz(:,41))./(r1+r2);
    bz(:,32)=(r1.*bz(:,22)+r2.*bz(:,42))./(r1+r2);
    bz(:,33)=(r1.*bz(:,23)+r2.*bz(:,43))./(r1+r2);
    bz(:,35)=(bz(:,34)+bz(:,36))./2;
    bz(:,39)=(bz(:,38)+bz(:,40))./2;
    %bz(:,58)=(bz(:,57)+bz(:,59))./2;
    bz(:,63)=(bz(:,62)+bz(:,64))./2;
    bz(:,95)=(bz(:,94)+bz(:,96))./2;
    ok_bz_plot=ok_bz;
    ok_bz([22 26 28 31 32 33 35 39 63 95])=true;
else
    r1=0.021;
    r2=0.0535;
    bz(:,22)=(bz(:,21)+bz(:,23))./2;
    bz(:,31)=(r1.*bz(:,21)+r2.*bz(:,41))./(r1+r2);

    %new
    bz(:,16)=(bz(:,15)+bz(:,17))./2;
    bz(:,27)=(r1.*bz(:,17)+r2.*bz(:,37))./(r1+r2);
    bz(:,26)=(bz(:,25)+bz(:,27))./2;
    bz(:,42)=(bz(:,41)+bz(:,43))./2;
    bz(:,32)=(r1.*bz(:,22)+r2.*bz(:,42))./(r1+r2);
    
   
    bz(:,61)=(2.*bz(:,51)+bz(:,81))./3;
    bz(:,71)=(bz(:,51)+2.*bz(:,81))./3;
    bz(:,68)=(bz(:,67)+bz(:,69))./2;
    bz(:,100)=(bz(:,41)+bz(:,43))./2;
    
    %

    bz(:,33)=(r1.*bz(:,23)+r2.*bz(:,43))./(r1+r2);
    bz(:,35)=(bz(:,34)+bz(:,36))./2;
    bz(:,39)=(bz(:,38)+bz(:,40))./2;
    bz(:,58)=(bz(:,57)+bz(:,59))./2;
    bz(:,63)=(bz(:,62)+bz(:,64))./2;
    bz(:,95)=(bz(:,94)+bz(:,96))./2;
    ok_bz_plot=ok_bz;
    %ok_bz([22 31 33 35 39 58 63 95])=true;
    ok_bz([16 27 26 32 42 61 71 100 68 22 31 33 35 39 58 63 95])=true;
end
end