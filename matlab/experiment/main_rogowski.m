date = 230927;%230523;
shot = 12;%8;
ch = [1 3 4 7 8];
calib(ch) = [116.6647,225.71,174.19,63.9568,223.2319];

rgw2txt(date,shot)
rogowski_plot(date,shot,ch,calib)