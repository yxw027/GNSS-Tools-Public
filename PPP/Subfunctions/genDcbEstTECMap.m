%%
function dcbDataOut = genDcbEstTECMap(dayNum,Year,statCode,consts,settings)
% generate DCB estimates given input measurements and TEC map.


close all force
nSkipSatEpochs = 1;


constellations = initConstellation(consts(1),consts(2),consts(3),consts(4),0);

settings.multiConst = consts;
settings.constellation = 'MULTI';

prns = constellations.PRN;
constInds = constellations.constInds';
nPrn = length(prns);

% reference pairs for each constellation
refpairs = {'C1W' 'C2W';
            'C1P' 'C2P'};


%% Read in the obs file
disp('Reading observation data...')
if ~exist('obsStruc','var')
    AttachToMice;
        
    filenameObs = genIgsStationObsDataName(statCode,Year,dayNum,settings);
    
    if isempty(filenameObs)
        % try to download the file
        ftpHelper(8,Year,dayNum,settings,{statCode})
        
        % look for it again
        filenameObs = genIgsStationObsDataName(statCode,Year,dayNum,settings);
        
        if isempty(filenameObs)
            error(':(')
        end
    end
        
    % g r e c
%     consts = [useGPS useGLO useGAL useBDS];
    constellations = initConstellation(consts(1),consts(2),consts(3),consts(4),0);
    
    % open the obs file (may be compressed)
    [obsStruc, time_ref, epochs, week, date, rxPos0, interval, antoff, ...
        antMod, rxmod, obsColumns] = readCompObsFile(filenameObs,constellations,settings);
        
    usrPos = repmat(rxPos0,1,length(epochs));
    
    jdu = unique(floor(epochs2jd(epochs)-0.5)+0.5);
    [dayList,YearList] = jd2doy(jdu);
    
    nEpochs = length(epochs);
    
    % collect all of the code phase measurements
    obsNamesFull = fields(obsStruc);
    obsOut   = zeros(0,nEpochs,size(obsStruc.(obsNamesFull{1}),1));
    obsTypes = {};
    freqDes = zeros(0,size(obsStruc.(obsNamesFull{1}),1));
    for idx = 1:length(obsNamesFull)
        obsNamei = obsNamesFull{idx};
        
        if ~strcmp(obsNamei(1),'C')
            % only interested in code phase- sorry buddy!
            continue;
        end
        
        % collect this one
        obsOut = cat(1,obsOut,permute(obsStruc.(obsNamei),[3 2 1]));
        obsTypes = [obsTypes; repmat({obsNamei},1,size(obsStruc.(obsNamesFull{1}),1))];
        
        freqDes = [freqDes; repmat(str2double(obsNamei(2)),1,size(obsStruc.(obsNamesFull{1}),1))];
    end
    
    nSignals = size(obsOut,1);
    obsOut(obsOut == 0) = nan;
    
    % only keep the constellations and PRNs of interest
%     indsKeepObs = ismember(obsDm.constInds,constInds);
%     obsOut = obsOut(:,:,indsKeepObs);
%     obsTypes = obsTypes(:,indsKeepObs);
%     freqDes = freqDes(:,indsKeepObs);
    
    [freqs, freqInds] = mapSignalFreq2(freqDes',prns,constInds,jdu(1));
    
    prph12   = obsOut;
    prphType = obsTypes;
    prphSig  = ones(size(obsTypes));
    prphInd  = repmat([1:size(obsOut,1)]',1,size(obsTypes,2));
    
    
    %% Check if we need to change our reference signals
%     refpairs = {'C1W' 'C2W';
%             'C1P' 'C2P'};
    if ~ismember('C1P',obsNamesFull)
        refpairs{1,1} = 'C1W';
    end
    if ~ismember('C2P',obsNamesFull)
        refpairs{1,2} = 'C2W';
    end
    
    
    %% Add dual frequency measurements
    prDf = zeros(nPrn,nEpochs);
    refSig = zeros(nPrn,nEpochs);
    freqDf = zeros(nPrn,1);
    for idx = 1:nPrn
        % find the observation we care about
        indObs1 = strFindCell(obsTypes(:,idx),refpairs(constInds(idx),1));
        indObs2 = strFindCell(obsTypes(:,idx),refpairs(constInds(idx),2));
        
        freq1 = freqs(idx,indObs1);
        freq2 = freqs(idx,indObs2);
        
        pr1 = squeeze(obsOut(indObs1,:,idx));
        pr2 = squeeze(obsOut(indObs2,:,idx));
        pr1(pr1 == 0) = nan;
        pr2(pr2 == 0) = nan;
        
        prDf(idx,:) = (pr1.*freq1.^2-pr2.*freq2.^2)./(freq1.^2-freq2.^2);
        freqDf(idx) = (freq1.^2-freq2.^2)./(freq1-freq2);
        
        refSig(idx,:) = pr1;
    end
    
    %% Load dcb data
    %     % Check if we need it in the first place
    YearListDcb = YearList;
    dayListDcb = dayList;
    epochDcb = jd2epochs(doy2jd(YearListDcb,dayListDcb))+100;
    [~,filenameDcb] = loadDcb(YearListDcb,dayListDcb,settings,1,2);
    
    %     if ~isempty(statCode) % non-igs station receiver
    %         [~,filenameDcb2] = loadDcb(YearListDcb,dayListDcb,settings,1,4);
    %         if isempty(filenameDcb2); filenameDcb = [];end
    %     end
    %     if isempty(filenameDcb)
    %         disp('Downloading new DCB data...')
    %         ftpHelper(6,YearListDcb,dayListDcb,settings);
    %     end
    %     disp('Parsing DCB data...')
    [dcbData,filenameDcb] = loadDcb(YearListDcb,dayListDcb,settings,0,4);
    %     if isempty(statCode)
    %         [dcbData2,filenameDcb] = loadDcb(YearListDcb,dayListDcb,settings,0,2);
    %     end
    
    %% Antenna phase center offset data (satellite and receiver)
    filenameAtx = [cd '\Utility Functions\igs14.atx'];
    
    atxData = ReadATX(filenameAtx);
    % Find antenna data for the receiver of interest
    atxRx = atxData( find(~cellfun(@isempty,strfind({atxData.block},antMod))));
    
    % remove data that isn't for satellites
    atxData = atxData([atxData(:).type] ~= 0);
    
    %% Load precise ephemeris data (check locally then download if necessary)
    [dayListPad,YearPad] = jd2doy(unique([jdu jdu+1 jdu-1]));
    
    [~,~,PFileNameFull] = loadPEph(YearPad,dayListPad,settings,1,atxData,0);
    
    dlFlag = 0;
    for idx = 1:length(PFileNameFull)
        if ~exist(PFileNameFull{idx},'file')
            dlFlag = 1;
        end
    end
    if dlFlag
        disp('Downloading new precise ephemeris...')
        % Download necessary orbit products
        for idx = find(consts)
            ftpHelper(1,YearPad,dayListPad,settings,settings.orbCenter{idx});
        end
    end
    disp('Loading precise ephemeris...')
    
    % Read in new or existing products
    [Peph,~,PFileNameFull] = loadPEph(YearPad,dayListPad,settings,0,atxData,0);
    
    %%
    if 1
        [~,~,~,IFileNameFull] = loadIonex(YearList,dayList,settings,1);
        dlFlag = 0;
        for idx = 1:length(IFileNameFull)
            if ~exist(IFileNameFull{idx},'file')
                dlFlag = 1;
            end
        end
        if dlFlag
            % Download necessary orbit products
            ftpHelper(15,YearPad,dayListPad,settings,settings.clkCenter{idx});
        end
        ionoData = loadIonex(YearList,dayList,settings,0);
    else
        ftpHelper(15,YearList,dayList,settings);
        ionoData = loadCodeKlobIono(YearList,dayList,settings);
    end
end

%% compute the measured L1-L.IF difference
% indsDfCode = find(prphSig > 100 & prphInd == 1);
% nDf = length(indsDfCode);
% ionoMeas = prph12(1,:,:)-prph12(indsDfCode,:,:);
measAvail = squeeze(any(prph12));

refFreq = 1575.42e6;


%% satellite positions and az/el (rough. might do just every 10 seconds... only care for az/el)
if ~exist('az','var')
    az = nan(nPrn,nEpochs);
    el = nan(nPrn,nEpochs);
    latPp = nan(nPrn,nEpochs);
    lonPp = nan(nPrn,nEpochs);
    obliq = nan(nPrn,nEpochs);
    ionoDelay = nan(nSignals,nPrn,nEpochs);
    
    tecMap = squeeze(ionoData.tecMap);
    
    satPos = nan(nPrn,nEpochs,3);
    tdxPrev = 1;
    for tdx= 1:nSkipSatEpochs:nEpochs
        sIndsi = find(measAvail(tdx,:));
        
        satPos(sIndsi,tdx,:) = PPosInterp(prns(sIndsi)', epochs(tdx)*ones(size(sIndsi)),Peph.position,Peph.PRN,...
            Peph.epochs,settings,NaN,[],[],constInds(sIndsi),Peph.constellation,1,atxData);
        
        [el(sIndsi,tdx),az(sIndsi,tdx)] = pos2elaz(usrPos(:,tdx)', squeeze(satPos(sIndsi,tdx,:)));
        
        usrLlh = xyz2llh(usrPos(:,tdx)');
        
        for sdx = 1:length(sIndsi)
            [latPp(sIndsi(sdx),tdx), lonPp(sIndsi(sdx),tdx),obliq(sIndsi(sdx),tdx)] = iono_pierce_point(usrLlh(1)*pi/180, ...
                usrLlh(2)*pi/180,az(sIndsi(sdx),tdx)',el(sIndsi(sdx),tdx)');
            
            teci =  -interpn(ionoData.latVec,ionoData.lonVec,ionoData.epochs,tecMap,latPp(sIndsi(sdx),tdx)*180/pi,lonPp(sIndsi(sdx),tdx)*180/pi,epochs(tdx));
            
            ionoDelay(:,sIndsi(sdx),tdx) = obliq(sIndsi(sdx),tdx)*teci*40.3*10^15./freqs(sIndsi(sdx),:).^2;
        end
        
        satPos(sIndsi,tdxPrev:tdx,:) = repmat(satPos(sIndsi,tdx,:),1,tdx-tdxPrev+1,1);
        el(sIndsi,tdxPrev:tdx) = repmat(el(sIndsi,tdx),1,tdx-tdxPrev+1,1);
        az(sIndsi,tdxPrev:tdx) = repmat(az(sIndsi,tdx),1,tdx-tdxPrev+1,1);
        ionoDelay(:,sIndsi,tdxPrev:tdx) = repmat(ionoDelay(:,sIndsi,tdx),1,1,tdx-tdxPrev+1);
        
        tdxPrev = tdx;
    end
end


%% iono corrections based on TEC map
ionoMeas = prph12-permute(repmat(prDf,1,1,nSignals),[3 2 1]);

% ionoMeas = prph12-permute(repmat(refSig+squeeze(ionoDelay(1,:,:)),1,1,nSignals),[3 2 1]);

measAvail = squeeze(any(ionoMeas));

ionoMeas = permute(ionoMeas,[1 3 2]);



% take the delta between the predicted and measured iono dleay

% dIono = ionoMeas+repmat(permute(ionoDelay,[3 2 1]),2,1,1);

dIono = ionoMeas+ionoDelay;

% plot(squeeze(dIono(4,33:end,1:30:end))')
% plot(squeeze(dIono(10,1:32,1:30:end))')


dcbEst = nanmean(dIono,3);
ionoStd = nanstd(dIono,1,3);

figure(1); clf; hold on;
plot(dcbEst','o','linewidth',2);
legend(prphType(:,1))

%% actually make sure doing this aligns stuff
% use the DCB estimates to align everything to the P1-P2 combination.
% prAligned = permute(prph12,[1 3 2])+ionoDelay;
prAligned = permute(prph12,[1 3 2])+ionoDelay-repmat(dcbEst,1,1,nEpochs);

prDiff2 = permute(repmat(prDf,1,1,nSignals),[3 1 2])-prAligned;

if 1
   figure(10); clf; hold on;
   plot(squeeze(prDiff2(1,:,:))')
end

prDiffMean2 = nanmean(prDiff2,3);
figure(11); clf; hold on;
hist(prDiff2(:),100)


% make dual frequency combinations of these now
indsDf = [1 8];
prDf2 = (squeeze(prAligned(indsDf(1),:,:)).*freqs(:,indsDf(1)).^2-...
    squeeze(prAligned(indsDf(2),:,:)).*freqs(:,indsDf(2)).^2)./(freqs(:,indsDf(1)).^2-freqs(:,indsDf(2)).^2);

dfDiff = prDf-prDf2;
dfDiffMean = nanmean(dfDiff,2);


%% write to a dcb structure
startEpoch = epochs(1);
endEpoch   = epochs(end);

[startDay, startYr] = jd2doy(epochs2jd(startEpoch));
[endDay, endYr]     = jd2doy(epochs2jd(endEpoch));
startDay = round(startDay);
endDay   = round(endDay);

jdi = epochs2jd(startEpoch);
indMeas = 1;
dcbDataOut.biasMode = 'ABS';
dcbDataOut.sites = [];
constLetters = 'GRECJ';
for idx = 1:size(dcbEst,2)
    
    prni = prns(idx);
    consti = constInds(idx);
    svni  = prn2svn(prni,jdi,consti);
    
    for jdx = 1:size(dcbEst,1)
        if isnan(dcbEst(jdx,idx))
            continue;
        end
        
        dcbDataOut.SVNs(indMeas,1) = svni;
        dcbDataOut.PRNs(indMeas,1) = prni;
        dcbDataOut.domes(indMeas,1) = NaN;
        dcbDataOut.obs1(indMeas,1) = obsTypes(jdx,idx);
        dcbDataOut.obs2(indMeas,1) = {'ABS'};
        dcbDataOut.startYr(indMeas,1) = startYr;
        dcbDataOut.startDy(indMeas,1) = startDay;
        dcbDataOut.startEpoch(indMeas,1) = startEpoch;
        dcbDataOut.endYr(indMeas,1) = endYr;
        dcbDataOut.endDy(indMeas,1) = endDay;
        dcbDataOut.endEpoch(indMeas,1) = endEpoch;
        dcbDataOut.bias(indMeas,1) = dcbEst(jdx,idx)/settings.c*1e9;
        dcbDataOut.stdDev(indMeas,1) = NaN;
        dcbDataOut.satFlag(indMeas,1) = 1;
        dcbDataOut.const(indMeas,1) = {constLetters(consti)};
        dcbDataOut.constInd(indMeas,1) = consti;
        
        indMeas = indMeas+1;
    end
end

% save this one
DfileNameFormat = ['%4d/%03d/' statCode 'MGXRAP_%4d%03d0000_01D_01D_DCB.mat'];

% File location
DpathName = settings.dcbMgexDir;

% Filename
DFileName = sprintf(DfileNameFormat, startYr, startDay,startYr, startDay);

dcbData = dcbDataOut;
save([DpathName DFileName],'dcbData')

%%
%{
prnPlot = 8;
figure(10); clf; hold on;
plot(squeeze(ionoMeas(prnPlot,:,:)),'.')

plot(-ionoDelay(prnPlot,:)+ionoMean(prnPlot),'r','linewidth',2)


figure(11); clf; hold on;
% plot(ionoMean,'o')
% plot(-dcbCorr(1,:)*settings.c*10-6,'ro')
svns =  prn2svn(1:32,2017.10,1);
errorbar(svns,ionoMean(1:32),ionoStd(1:32),'o')
errorbar(svns,ionoMean2(1:32),ionoStd2(1:32),'ro')
xlabel('SVN')
ylabel('Iono correction - measured iono diff (TGD lol :()')


figure(12); clf; hold on;
% plot(ionoMean,'o')
% plot(-dcbCorr(1,:)*settings.c*10-6,'ro')
svns =  prn2svn(1:24,2017.10,2);
errorbar(svns,ionoMean(33:56),ionoStd(33:56),'o')
errorbar(svns,ionoMean2(33:56),ionoStd2(33:56),'ro')
xlabel('SVN')
ylabel('Iono correction - measured iono diff (TGD lol :()')


% xlim([0.5 32.5])


dcbTgd = [ionoMean; ionoMean2];

% fit something lol
% % A = [dcbCorr(1,1:32)'*settings.c ones(32,1)];
% % b = [ionoMean(1:32)-dcbCorr(1,1:32)*settings.c]';
% %
% % x = A\b;
%
%
% diff = ionoMean+(dcbCorr(1,:)*settings.c*x(1)-x(2));
% figure(12); clf; hold on;
% plot(diff(1:32),'o')




%}







