function deleteFolderContents(dirName)


diri = dir(dirName);

for idx = 1:length(diri)
    if ~diri(idx).isdir
        delete([dirName diri(idx).name]);
    end
    
end
   

end