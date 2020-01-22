function cbias = clockBiasFromProd(prns,constInds,epochs,orbClockData)

% check whether to use precise clock


if strcmp(orbClockData.clkMode,'PRECISE')
    cbias = CClockInterp(prns,constInds,epochs,orbClockData.PClock);
    
    
    
else
    cbias = nan(size(prns));
   % put broadcast clock stuff here 
      % do broadcast stuff
    constUn = unique(constInds);
    
    % convert all epochs to week number and TOW
    [weeks,tows] = epochs2gps(epochs);
    
    if ismember(constUn,1)
        % gps    
        indsi = find(constInds == 1);
        pos = SVPos(orbClockData.brdcGPS,prns(indsi),weeks(indsi),tows(indsi),'GPS');
        cbias(indsi) = pos.clock_bias-pos.TGD*0;
    end
    
    if ismember(constUn,2)
        % gps    
        indsi = find(constInds == 2);
        pos = SVPos(orbClockData.brdcGLO,prns(indsi),weeks(indsi),tows(indsi),'GLO');
        cbias(indsi) = pos.clock_bias-pos.TGD*0;
    end
    
    
end











end