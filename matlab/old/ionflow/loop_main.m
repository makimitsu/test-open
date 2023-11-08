for t = 1:9
    for ndata = 1:5
        main(463+2*t,ndata,1,1,false,true,true,false,0)
        a = 5*(t-1)+ndata
    end
end
