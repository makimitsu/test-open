function [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname)
filename = [pathname.flowdata,'/',num2str(date),'/shot',num2str(ICCD.shot),'_',num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.mat'];
if exist(filename,"file")
    load(filename,'V_i','absV','T_i')
else
    warning(strcat(filename,' does not exist.'));
    V_i = char.empty;
    absV = char.empty;
    T_i = char.empty;
    return
end
end

