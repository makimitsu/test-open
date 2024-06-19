%NIFSのデータをドライブにコピーする

function copyFolderIfNotExist(destFolder, sourceFolder)
    if ~isfolder(destFolder)
        try
            disp('NIFS data copying to Drive');
            copyfile(sourceFolder, destFolder);
        catch
            error('failed to copy NIFS data to Drive');
        end
    end
end