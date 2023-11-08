function [z_IDS,r_IDS,Ti_IDS,Ti_max_IDS,Ti_min_IDS,Em_IDS] = load_Doppler288ch(date,IDS288ch,pathname)
filename = [pathname.tempdata,'/',num2str(date),'/shot',num2str(IDS288ch.shot),'_',num2str(IDS288ch.trg),'us_w=',num2str(IDS288ch.exp_w),'.mat'];
if exist(filename,"file")
    load(filename,'z_IDS','r_IDS','Ti_IDS','Ti_max_IDS','Ti_min_IDS','Em_IDS')
else
    warning(strcat(filename,' does not exist.'));
    z_IDS = char.empty;
    r_IDS = char.empty;
    Ti_IDS = char.empty;
    Ti_max_IDS = char.empty;
    Ti_min_IDS = char.empty;
    Em_IDS = char.empty;
    return
end
end
