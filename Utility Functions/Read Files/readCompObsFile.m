function [obsStruc, time_ref, time, week, date, rxPos0, interval, antoff, ...
    antmod, rxmod, obsColumns] = readCompObsFile(filename,constellations,settings)

% function to open and read a compressed observation file

if nargin < 3
    settings = [];
end

% check if this is fully uncompressed observation file.
[dirPath,filenamei,exti] = fileparts(filename);

zipExts = {'.GZ','.Z'};
hatanakaExt = {'.CRX'};


% check if anything needs to be done to this file (unzip, decompress).  if
% so, setup a temporary folder to work in.
if contains(upper(filename),zipExts) || contains(upper(filename),hatanakaExt)
    % do we need to unzip first?
    if ~isempty(settings)
        tempDir = [settings.tempDir 'tempDir' num2str(randi(1000)) '\'] ;
    else
        tempDir = [cd '\tempDir' num2str(randi(1000)) '\'];
    end
    
    % create the directory
    mkdir(tempDir);
    
    if contains(upper(filename),zipExts)
        unzipFile(filename,tempDir);
        
        % look in the new folder for the unzipped file
        diri = dir(tempDir);
        
        indFile = find(~[diri.isdir]);
        
        filename = [tempDir diri(indFile).name];
    else
        % move it to the temp dir either way
        movefile(filename,tempDir);
    end
    
    if contains(upper(filename),hatanakaExt)
        % decompress the rinex
        crx2rnx(filename);
        [~,justNameCrx] = fileparts(filename);
        filename = [tempDir justNameCrx '.rnx'];
    end
else
    tempDir = [];
end


% ok actually open the file
[obsStruc, time_ref, time, week, date, rxPos0, interval, ...
    antoff, antmod, rxmod, obsColumns] = load_RINEX_obs_fullObs(filename,constellations,0,1);

% clean up the mess you made you dirty function
if ~isempty(tempDir)
    deleteFolderContents(tempDir);
    rmdir(tempDir);
end

end















