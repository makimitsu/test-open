function cell = get_ch_dist(dir,date)

if date >= 191209 && date < 191218
    cell = readcell(dir,'Sheet','Bz_191209');
elseif date >= 191218 && date < 200121
    cell = readcell(dir,'Sheet','Bz_191218');
elseif date >= 200121
    cell = readcell(dir,'Sheet','Bz_200121');
else
    disp('Wrong date!!!!');
    return;
end

end