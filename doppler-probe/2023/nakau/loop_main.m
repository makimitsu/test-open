for t = 1:4
    for ndata = 1:4
%         main(462+4*t,ndata,1,1,true,true,true,false,0)
        main(462+4*t,ndata,1,1,false,true,true,false,0)
%         main_96(462+4*t,ndata,1,1,true,true,true,false,0)
%         main_96(462+4*t,ndata,1,1,false,true,true,false,0)
        a = 4*(t-1)+ndata
    end
end
