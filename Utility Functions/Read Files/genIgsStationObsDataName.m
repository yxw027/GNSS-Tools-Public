function filenameOut = genIgsStationObsDataName(statCode,Year,dayNum,settings)

% function to look for and produce name of obs data from IGS mgex stations
% (files contained locally)
obsMatDir = settings.mgxObsDir;
dirName = [ obsMatDir int2str(Year) '/' num2str(dayNum,'%03d') '/'];

if ~exist(dirName,'dir')
    % directory wasn't found, so the file can't be in there!
    filenameOut = [];
else
    % look in the directory for the file we want
    diri = dir(dirName);
    
    filenames = {diri.name};
    
    indsStatCode = strFindCell(filenames,statCode);
    if ~isempty(indsStatCode)
        % found it!
        filenameOut = [dirName diri(indsStatCode(1)).name];
    else
        % didn't find anything in there!
       filenameOut = []; 
    end
end

end