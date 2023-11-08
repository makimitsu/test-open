function [grid2D,data2D,ok_z,ok_r] = load_pcb200ch_with_tfoffset(date,dtacq,pathname)
filename = [pathname.processeddata,'/',num2str(date),'/processeddata_dtacq',num2str(dtacq.num), ...
        '_shot',num2str(dtacq.shot),'_tfshot0.mat'];
if exist(filename,"file")
    load(filename,'grid2D','data2D','ok_z','ok_r')
else
    warning(strcat(filename,' does not exist.'));
    grid2D = char.empty;
    data2D = char.empty;
    ok_z = char.empty;
    ok_r = char.empty;
    return
end
end
