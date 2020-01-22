function [obsGnss dcbCorr] = preprocessGnssObs(obsGnssRaw,obsInds,signalInds,obsDes,...
    ifPairs,settings,PARAMS,varargin)



%%
epochs = obsGnssRaw.epochs;
prns   = obsGnssRaw.PRN;
constInds = obsGnssRaw.constInds;

nPrn = length(prns);

jdu = unique(floor(epochs2jd(epochs)-0.5)+0.5);
[dayList,YearList] = jd2doy(jdu);

nEpochs = length(epochs);


nObs = size(obsDes,2);
obsOut = zeros(nObs,nEpochs,nPrn);
obsTypes = cell(nObs,nPrn);
for pdx = 1:nPrn
    prni = prns(pdx);
    consti = constInds(pdx);
    for odx = 1:nObs
        obsSeti = obsDes{consti,odx};
        for osdx = 1:length(obsSeti)
            obsTypei = obsSeti{osdx};
            if ~isfield(obsGnssRaw.meas,obsTypei)
                continue
            end
            if ~any(obsGnssRaw.meas.(obsTypei)(pdx,:))
                continue
            end
            % Save it off.
            obsOut(odx,:,pdx) = obsGnssRaw.meas.(obsTypei)(pdx,:);
            obsTypes{odx,pdx} = obsTypei;
            
            % Don't need to look for additional data
            break
        end
    end
end

freqDes2 = nan(size(obsTypes));
for idx = 1:size(freqDes2,1)
    for jdx = 1:size(freqDes2,2)
        if ~isempty(obsTypes{idx,jdx})
            freqDes2(idx,jdx) = str2double(obsTypes{idx,jdx}(2));
        end
    end
end
[freqs, freqInds] = mapSignalFreq2(freqDes2(obsInds == 1 | obsInds == 2,:)',prns,constInds,jdu(1));

[freqsDop, freqIndsDop] = mapSignalFreq2(freqDes2(obsInds == 4,:)',prns,constInds,jdu(1));


prph12   = squeeze(obsOut(obsInds == 1 | obsInds == 2,:,:));
prphType = squeeze(obsTypes(obsInds == 1 | obsInds == 2,:));
prphSig  = signalInds(obsInds == 1 | obsInds == 2);
prphInd  = obsInds(obsInds == 1 | obsInds == 2);

snr12    = squeeze(obsOut(obsInds == 3 ,:,:));
snrType  = squeeze(obsTypes(obsInds == 3,:));
snrSig   = signalInds(obsInds == 3);

dop12    = squeeze(obsOut(obsInds == 4,:,:));
dopType  = squeeze(obsTypes(obsInds == 4,:));
dopSig   = signalInds(obsInds == 4);

% convert carrier phase from cycles to meters
for idx = 1:size(prph12,1)
    if prphInd(idx) == 2
        prph12(idx,:,:) = squeeze(prph12(idx,:,:)).*settings.c./freqs(:,idx)';
    end
end
% convert doppler from cycles/second to meters/second
for idx = 1:size(dop12,1)
    dop12(idx,:,:) = squeeze(dop12(idx,:,:)).*settings.c./freqs(:,min(find(prphSig == dopSig(idx))))';
end


%% Load dcb data
YearListDcb = YearList;
dayListDcb = dayList;
epochDcb = jd2epochs(doy2jd(YearListDcb,dayListDcb))+100;

% Use IGS if possible
dcbType = 3; % 1 = CODE, 0/2 = SU, 3 = DLR
% if ~isempty(statCode)
%     [~,filenameDcb] = loadDcb(YearListDcb,dayListDcb,settings,1,4);
%     
%     % If it wasn't found locally, check online
%     if ~exist(filenameDcb{1},'file')
%         disp('Downloading new DCB data...')
%         ftpHelper(6,YearListDcb,dayListDcb,settings);
%     end
%     
%     % Try to parse something
%     [dcbData,filenameDcb] = loadDcb(YearListDcb,dayListDcb,settings,0,4);
%     
%     if ~isempty(dcbData)
%         dcbType = 1;
%     end
% end

% Still need to pull SU estimates
if dcbType == 0
    [~,filenameDcb] =  loadDcb(YearListDcb,dayListDcb,settings,1,5);
    
    if ~exist(filenameDcb{1},'file')
        % generate a dcb file
        
        consts = settings.multiConst;
        
        dcbData = genDcbEstTECMap(dayListDcb,YearListDcb,'STFU',consts,settings);
    else
        [dcbData,filenameDcb] =  loadDcb(YearListDcb,dayListDcb,settings,0,5);
    end
    dcbType = 2;
elseif dcbType == 1
    % CODE
   dcbData = loadDcb(YearListDcb,dayListDcb,settings,0,4);
elseif dcbType == 3
    % DLR
    dcbData = loadDcb(YearListDcb,dayListDcb,settings,0,2);
end
    


%% Correct L1C-L1P for ISC (using GPS and Galileo MGEX precise products, this should be the only necessary change)
%     dcbData2 = [];
if dcbType == 1
    % CODE
    dcbCorr = zeros(size(prphType));
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if  consti == 1
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,{'AJAC'});
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if strcmp(obsi{1},'C1P') && consti == 2 && ~isempty(statCode)
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
            
            if strcmp(obsi{1},'C2P') && consti == 2 && ~isempty(statCode)
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
        end
    end
    
    dcbCorr = settings.dcbUse*dcbCorr*settings.c;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    % Apply the corrections to the observations
    prph12i = prph12-permute(repmat(dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 3
    % DLR 
    % these are relative corrections that need to be further referenced to
    % the L1P-L2P combination (need the TGD term)
    dcbCorr = zeros(size(prphType));
    
%     for idx = 1:size(prphType,1)
%         for jdx = 1:size(prphType,2)
%             prni = prns(jdx);
%             consti = constInds(jdx);
%             obsi = prphType(idx,jdx);
%             
%             if isempty(obsi{1}) | ~strcmp(obsi{1}(1),'C')
%                 continue;
%             end
%             freqi = str2num(obsi{1}(2));
%             
%             if  consti == 1
%                 biasi = 0;
%                  
%                 % If this isn't C1C, need a correction to get to C1C
%                 if ~strcmp(obsi{1},'C1C')
%                     if strcmp(obsi{1}(2),'2')
%                         % this is L2 signal- adjust to C2W then to C1C
%                         if strcmp(obsi{1},'C2W')
%                             dcb_c2x_c2w = 0;
%                         else
%                             dcb_c2x_c2w = findDcbElement(prni,consti,{'C2W'},obsi,epochDcb,dcbData);
%                         end
%                         
%                         dcb_c2w_c1c = findDcbElement(prni,consti,{'C1C'},{'C2W'},epochDcb,dcbData);
%                            
%                         dcb_xyz_c1c = -dcb_c2x_c2w+dcb_c2w_c1c;
%                     end
%                         
%                 else
%                    dcb_xyz_c1c = 0; 
%                 end
%                 
%                 % C1C to C1W correction
%                 dcb_c1c_c1w = findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData);
%                 
%                 % Add the TGD term no matter what- gets from L1P to L1P/L2P
%                 dcb_c1c_c2w = findDcbElement(prni,consti,{'C1C'},{'C2W'},epochDcb,dcbData);
%                 dcb_c1w_c2w = dcb_c1c_c2w-dcb_c1c_c1w;
%                 gam12 = 1575.42^2./1227.6^2;
%                 tgd = dcb_c1w_c2w./(1-gam12);
%                 
%                 % sum up the corrections
%                 biasi = dcb_xyz_c1c+dcb_c1c_c1w+tgd;
%                 
%                 dcbCorr(idx,jdx) = biasi;
%             end
%             
%             if consti == 2
%                 biasi = 0;
%                  
%                 % If this isn't C1C, need a correction to get to C1C
%                 if ~strcmp(obsi{1},'C1C')
%                     
%                 else
%                    dcb_xyz_c1c = 0; 
%                 end
%                 
%                 % C1C to C1W correction
%                 dcb_c1c_c1w = findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
%                 
%                 % Add the TGD term no matter what- gets from L1P to L1P/L2P
%                 dcb_c1c_c2w = findDcbElement(prni,consti,{'C1C'},{'C2P'},epochDcb,dcbData);
%                 dcb_c1w_c2w = dcb_c1c_c2w-dcb_c1c_c1w;
%                 gam12 = 1602^2./1246^2;
%                 tgd = dcb_c1w_c2w./(1-gam12);
%                 
%                 % sum up the corrections
%                 biasi = dcb_xyz_c1c+dcb_c1c_c1w+tgd;
%                 
%                 dcbCorr(idx,jdx) = biasi;
%                 dcbCorr(idx,jdx) = biasi;
%             end
%             
%         end
%     end


    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if isempty(obsi{1}) | ~strcmp(obsi{1}(1),'C')
                continue;
            end
            freqi = str2num(obsi{1}(2));
            
            if  consti == 1
                switch obsi{1}
                    case 'C1C'
                        biasi = findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData);
                    case 'C2W'
                        biasi = 0;
                    case 'C2S'
                        biasi = -findDcbElement(prni,consti,{'C2W'},{'C2S'},epochDcb,dcbData);
                    otherwise 
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                 switch obsi{1}
                    case 'C1C'
                        biasi = findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
                    case 'C2P'
                        biasi = 0;
                    case 'C2C'
                        biasi = -findDcbElement(prni,consti,{'C1C'},{'C2C'},epochDcb,dcbData);
                    otherwise 
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
        end
    end
    
    
    dcbCorr = dcbCorr*settings.c;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    
    % Apply the corrections to the observations
    prph12i = prph12-PARAMS.dcbUse*permute(repmat(dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 5
    dcbCorr = zeros(size(prphType));
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if strcmp(obsi{1},'C1C') && consti == 1
                %                     biasi = findDcbElement(prni,consti,obsi,{'C1W'},epochDcb,dcbData);
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if strcmp(obsi{1},'C2S') && consti == 1
                %                     biasi = findDcbElement(prni,consti,obsi,{'C1W'},epochDcb,dcbData);
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if strcmp(obsi{1},'C1P') && consti == 2 && ~isempty(statCode)
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
            
            if strcmp(obsi{1},'C2P') && consti == 2 && ~isempty(statCode)
                biasi = findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
        end
    end
    
    
    dcbCorr(isnan(dcbCorr)) = 0;
    % Apply the corrections to the observations
    prph12i = prph12-settings.c*permute(repmat(dcbUse*dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 2
    % STANFORD 
    dcbCorr = zeros(size(prphType));
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if ~isempty(obsi{1})
                if consti == 1 && strcmp(obsi(1),'C2S')
                    obsi(1) = {'C2X'};
                end
                
                if consti == 2
                    'fdafad';
                end
                
                biasi = settings.c*findDcbElement(prni,consti,obsi(1),{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
        end
    end
    
%     dcbCorr(:,sum(~(isnan(dcbCorr) | dcbCorr == 0)) < 2) = nan;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    fullCorr = dcbCorr;
    
    prph12i = prph12-PARAMS.dcbUse*permute(repmat(fullCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
else
    disp('No DCB''s available :(')
    prph12i = prph12;
end

%% Add dual frequency measurements
for idx = 1:size(ifPairs,1)
    freq1 = freqs(:,prphSig == ifPairs(idx,1) & prphInd == 1)';
    freq2 = freqs(:,prphSig == ifPairs(idx,2) & prphInd == 1)';
    freqInd1 = freqInds(:,prphSig == ifPairs(idx,1) & prphInd == 1);
    freqInd2 = freqInds(:,prphSig == ifPairs(idx,2) & prphInd == 1);
    
    pr1 = squeeze(prph12i(prphSig == ifPairs(idx,1) & prphInd == 1,:,:));
    pr2 = squeeze(prph12i(prphSig == ifPairs(idx,2) & prphInd == 1,:,:));
    pr1(pr1 == 0) = nan;
    pr2(pr2 == 0) = nan;
    ph1 = squeeze(prph12i(prphSig == ifPairs(idx,1) & prphInd == 2,:,:));
    ph2 = squeeze(prph12i(prphSig == ifPairs(idx,2) & prphInd == 2,:,:));
    ph1(ph1 == 0) = nan;
    ph2(ph2 == 0) = nan;
    
    prifi = (pr1.*freq1.^2-pr2.*freq2.^2)./(freq1.^2-freq2.^2);
    phifi = (ph1.*freq1.^2-ph2.*freq2.^2)./(freq1.^2-freq2.^2);
    freqif = (freq1.^2-freq2.^2)./(freq1-freq2);
    
    prphSig = [prphSig 100*ifPairs(idx,1)+ifPairs(idx,2) 100*ifPairs(idx,1)+ifPairs(idx,2)];
    prphInd = [prphInd 1 2];
    
    prph12i = cat(1,prph12i,permute(prifi, [3 1 2]));
    prph12i = cat(1,prph12i,permute(phifi, [3 1 2]));
    
    freqs = [freqs freqif' freqif'];
    freqInds = [freqInds freqInd1*10+freqInd2 freqInd1*10+freqInd2];
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 1,:),prphType(prphSig == ifPairs(idx,2) & prphInd == 1,:));
    typeAdd(cellfun(@length,typeAdd) <= 3) = {[]};
    prphType = [prphType; typeAdd];
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 2,:),prphType(prphSig == ifPairs(idx,2) & prphInd == 2,:));
    typeAdd(cellfun(@length,typeAdd) <= 3) = {[]};
    prphType = [prphType; typeAdd];
end

prph12i(isnan(prph12i)) = 0;
dop12(isnan(dop12)) = 0;
snr12(isnan(snr12)) = 0;

%% Arrange the lock time information (cycle slip detection) for easier use
if ~isempty(obsGnssRaw.tLock)
    lockTime = nan(size(prph12i));
    % Loop through each satellite
    for pdx = 1:length(prns)
        indSat = find(obsGnssRaw.PRN == prns(pdx) & obsGnssRaw.constInds == constInds(pdx));
        % Loop through each signal
        for jdx = 1:size(prph12i,1)
            %             indSignal
            sigNamei = prphType{jdx,pdx};
            
            if isempty(sigNamei) %| ~strcmp(sigNamei(1),'L')
                % no obs available- skip this
                continue;
            end
            
            if length(sigNamei) == 3
                % single frequency
                lockTimei = obsGnssRaw.tLock.(sigNamei)(indSat,:);
            else
               % dual frequency- take the minimum of both 
               lockTimei = min([obsGnssRaw.tLock.(sigNamei(1:3))(indSat,:);
                   obsGnssRaw.tLock.(sigNamei(4:6))(indSat,:)]);
            end
            
            lockTime(jdx,:,pdx) = lockTimei;
        end
    end

    
else
    lockTime = [];
end



%% Put everything into the output structure
obsGnss.PRN       = prns;
obsGnss.constInds = constInds;
obsGnss.epochs    = epochs;
obsGnss.freqs     = freqs';

obsGnss.range.obs         = permute(prph12i,[1 3 2]);
obsGnss.range.rnxCode     = prphType;
obsGnss.range.ind         = repmat(prphInd',1,size(obsGnss.range.obs,2));
obsGnss.range.sig         = repmat(prphSig',1,size(obsGnss.range.obs,2));
obsGnss.range.freqs       = freqs';
obsGnss.range.PRN         = repmat(prns,size(obsGnss.range.obs,1),1);
obsGnss.range.constInds   = repmat(constInds,size(obsGnss.range.obs,1),1);
obsGnss.range.lockTime    = permute(lockTime,[1 3 2]);

obsGnss.doppler.obs       = permute(dop12,[1 3 2]);
obsGnss.doppler.rnxCode   = repmat(dopType',1,size(obsGnss.doppler.obs,2));
obsGnss.doppler.sig       = repmat(dopSig',1,size(obsGnss.doppler.obs,2));
obsGnss.doppler.freqs     = freqsDop';
obsGnss.doppler.PRN       = repmat(prns,size(obsGnss.doppler.obs,1),1);
obsGnss.doppler.constInds = repmat(constInds,size(obsGnss.doppler.obs,1),1);

obsGnss.snr.obs           = permute(snr12,[1 3 2]);
obsGnss.snr.rnxCode       = repmat(snrType',1,size(obsGnss.snr.obs,2));
obsGnss.snr.sig           = repmat(snrSig',1,size(obsGnss.snr.obs,2));
obsGnss.snr.PRN           = repmat(prns,size(obsGnss.snr.obs,1),1);
obsGnss.snr.constInds     = repmat(constInds,size(obsGnss.snr.obs,1),1);


end
