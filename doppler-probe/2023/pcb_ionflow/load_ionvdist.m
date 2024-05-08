function [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,ICCD,pathname)
filename = [pathname.vdistdata,'/',num2str(date),'/shot',num2str(ICCD.shot),'_',num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.mat'];
if exist(filename,"file")
    load(filename,'V_i','absV','T_i','F','W','P','Lambda','Vx','Vy','ppoints','Angle')
else
    warning(strcat(filename,' does not exist.'));
    V_i = char.empty;
    absV = char.empty;
    T_i = char.empty;
    F = char.empty;
    W = char.empty;
    P = char.empty;
    Lambda = char.empty;
    Vx = char.empty;
    Vy = char.empty;
    ppoints = char.empty;
    Angle = char.empty;
    return
end
end
