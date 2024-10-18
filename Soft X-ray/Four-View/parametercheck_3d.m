function [N_projection, N_grid3d, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck_3d(newProjectionNumber, newGridNumber)
% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も範囲に関しては読み込む
parameterFile = sprintf('/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view/parameters_3d%d%d.mat', newProjectionNumber, newGridNumber);

if isfile(parameterFile)
        disp(strcat('Loading parameter from ',which(parameterFile)))
        load(parameterFile,'N_projection','N_grid3d','gm3d','U3d','s3d','v3d','M3d','K3d');
            
        if newProjectionNumber ~= N_projection || newGridNumber ~= N_grid3d
            disp('Different parameters - Start calculation!');
            get_parameters_3d(newProjectionNumber,newGridNumber,parameterFile);
            load(parameterFile,'N_projection','N_grid3d','gm3d','U3d','s3d','v3d','M3d','K3d');
        end
else
    disp('No parameters - Start calculation!');
    get_parameters_3d(newProjectionNumber,newGridNumber,parameterFile);
    load(parameterFile,'N_projection','N_grid3d','gm3d','U3d','s3d','v3d','M3d','K3d');
end

end