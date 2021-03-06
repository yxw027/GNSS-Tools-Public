function [svPos,svVel,iod] = svPosFromProd(prns, epochs,orbClockData,settings,...
    pPosInds,pPosPoly,constInds,FLAG_APC_OFFSET,atxData,sunPos,dttx)

if nargin < 8
    FLAG_APC_OFFSET = 0;
end

if nargin < 9
    atxData = [];
end

if nargin < 10
    sunPos = [];
end

if nargin < 11
    dttx = [];
end

iod = [];

if strcmp(orbClockData.orbMode,'PRECISE')
    % this should be adjusted to allow for all of the inputs
    [svPos, svVel] = PPosInterp(prns, epochs,orbClockData.PEph.position,...
        orbClockData.PEph.PRN,orbClockData.PEph.epochs,settings,NaN,...
        [],[],constInds,orbClockData.PEph.constellation,FLAG_APC_OFFSET,...
        atxData,sunPos,dttx);
        
else
    % put broadcast orbit stuff here
    svPos = nan(length(prns),3);
    svVel = nan(length(prns),3);
    iod   = nan(length(prns),1);
    
    % do broadcast stuff
    constUn = unique(constInds);
    
    % convert all epochs to week number and TOW
    [weeks,tows] = epochs2gps(epochs);
    
    if ismember(constUn,1)
        % gps
        indsi = find(constInds == 1);
        pos = SVPos(orbClockData.brdcGPS,prns(indsi),weeks(indsi),tows(indsi),'GPS');
        svPos(indsi,:) = [pos.x pos.y pos.z];
        svVel(indsi,:) = [pos.x_dot pos.y_dot pos.z_dot];
        
        iod(indsi) = pos.IODC;
    end
    
    if ismember(constUn,2)
        % glonass
        indsi = find(constInds == 2);
        pos = SVPos(orbClockData.brdcGLO,prns(indsi),weeks(indsi),tows(indsi),'GLO');
        svPos(indsi,:) = [pos.x pos.y pos.z];
        svVel(indsi,:) = [pos.x_dot pos.y_dot pos.z_dot];
    end
    
end






end