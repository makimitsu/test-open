setenv('SXR_IMAGE_DIR','G:\My Drive\X-ray\Data\TIF');
setenv('SXR_RECONSTRUCTED_DIR','G:\My Drive\X-ray\Data\SXROUT');
setenv('SXR_MATRIX_DIR','G:\My Drive\X-ray\Data\SXROUT');
% --
% load('parameters.mat');
% load('G:\My Drive\X-ray\Data\SXROUT\NLF_NLR\231216\shot21\2.mat');
filter = 3;% どの視点を試験するか
% --

prompt = {'doSave:','doFilter:','doNLR:'};
definput = {'','',''};
if exist('doSave','var')
    definput{1} = num2str(doSave);
end
if exist('doFilter','var')
    definput{2} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{3} = num2str(doNLR);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
doSave = logical(str2double(cell2mat(answer(1))));
doFilter = logical(str2double(cell2mat(answer(2))));
doNLR = logical(str2double(cell2mat(answer(3))));
switch filter
    case 1
        EE = EE1;
        plot_sxr_multi(doSave,doFilter,doNLR,M,K,gm2d1,U1,s1,v1,EE,range);
    case 2
        EE = EE2;
        plot_sxr_multi(doSave,doFilter,doNLR,M,K,gm2d1,U2,s2,v2,EE,range);
    case 3
        EE = EE3;
        plot_sxr_multi(doSave,doFilter,doNLR,M,K,gm2d1,U3,s3,v3,EE,range);        
    case 4
        EE = EE4;
        plot_sxr_multi(doSave,doFilter,doNLR,M,K,gm2d1,U4,s4,v4,EE,range);
end
