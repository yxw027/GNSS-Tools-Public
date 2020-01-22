% Inertial test script
close all force;
runNo = 3;

RANGEA = true;
switch runNo
    case 1
        filenameGnss = [];
        filenameObs = [boxDir '\Research\ION 2019\Data\rawdata_inertial.ASC'];
        epochStart     = 1222702044+100;
        epochEnd       = epochStart+60*15;
        truthFile = [boxDir '\Research\ION 2019\Data\nrcan solution\rawdata_inertial.pos'];
        
        accScale = 1./512;
        gyroScale = 1./131.2*pi/180;
        
        imuOrder = [1 2 3];
        IMU_ARM = [-0.262 0.944 1.575]';
        RANGEA = true;
        
    case 2
        filenameGnss = [];
        filenameObs = [boxDir '\Novatel Data\IGS Data\STFUHR\STFU00USA_S_20172802345_15M_01S_MO.mat'];
        truthFile =  [-2700404.33125106 ;        -4292605.32023686   ;       3855137.49176658]; % igs
        
        truthFile =  [-2700404.33125106 ;  -4292605.32023686   ;   3855137.49176658];  %nrcan

        epochStart = 1191380400+30;
        epochEnd   = epochStart+30;
        PARAMS.IMU_ARM = [0 0 0]';        
    case 3
%         filenameGnss = [boxDir '\Novatel Data\7500_OpenSky\Nov_1500_01Mar18_OpenSky.18O'];
        filenameObs = [boxDir '\Novatel Data\7500_OpenSky\Nov_1500_01Mar18_OpenSky_IMU.GPS'];
        
        filenameGnss = [];
        
        truthFile = [boxDir '\Novatel Data\7500_OpenSky\Truth_IGS08_OpenSky.txt'];
        
        truthNrcan = [boxDir 'Novatel Data\7500_OpenSky\Nov_1500_01Mar18_OpenSky2_2019.pos'];

        epochStart = 1203968280;
        epochEnd = epochStart+60*60;
        
        accScale = 1./512;
        gyroScale = 1./131.2*pi/180;
        
        imuOrder = [-2 1 3];
        %         accScale = 2/100;
        %         gyroScale = 1./(100*100);
        IMU_ARM = [0.328 0.327 1.310]';
        RANGEA = true;
        
    case 4
        filenameGnss = [boxDir '\Novatel Data\7500_Midtown-Suburban\Nov_1500_01Mar18_MidTown-Suburban.18O'];
        filenameObs = [boxDir '\Novatel Data\7500_Midtown-Suburban\Nov_1500_01Mar18_MidTown-Suburban_IMU.GPS'];
        
        filenameGnss = [];
        
        truthFile = [boxDir '\Novatel Data\7500_Midtown-Suburban\Truth_IGS08_MidTown-Suburban.txt'];
        
        epochStart = 1203955500;
        epochEnd= epochStart+60*60;
        
        accScale = 1./512;
        gyroScale = 1./131.2*pi/180;
        
        imuOrder = [-2 1 3]; % original
        
        IMU_ARM = [0.328 0.327 1.310]';
        RANGEA = true;
    case 5 % highway
        filenameObs = [boxDir '\Novatel Data\7500_Highway\Nov_1500_01Mar18_Highway_IMU.GPS'];
        filenameGnss = [boxDir '\Novatel Data\7500_Highway\Nov_1500_01Mar18_Highway.18O'];
        filenameGnss = [];
        
        gloFiti = [];
        statCode = [];
        truthFile = [boxDir '\Novatel Data\7500_Highway\Truth_IGS08_Highway.txt'];
        
        epochStart = 1203964680+30;
        epochEnd= epochStart+60*60;
        
        accScale = 1./512;
        gyroScale = 1./131.2*pi/180;
        
        imuOrder = [-2 1 3]; % original
        
        IMU_ARM = [0.328 0.327 1.310]';
    case 6 % downtown
        filenameGnss = [boxDir '\Novatel Data\7500_Downtown\Nov_1500_01Mar18_Downtown.18O'];
        filenameObs = [boxDir '\Novatel Data\7500_Downtown\Nov_1500_01Mar18_Downtown_IMU.GPS'];
        
        filenameGnss = [];
        
        truthFile = [boxDir '\Novatel Data\7500_Downtown\Truth_IGS08_Downtown.txt'];
        
        epochStart = 1203959160+1;
        epochEnd= epochStart+15*60;
        accScale = 1./512;
        gyroScale = 1./131.2*pi/180;
        
        imuOrder = [-2 1 3]; % original
        
        IMU_ARM = -1*[-0.262 0.944 1.575]';
end

useImu = true;

settings = initSettings;
settings.constellation = 'MULTI';
settings.multiConst = [1 1 0 0 0];
constellations = initConstellation(settings.multiConst(1),settings.multiConst(2),...
    settings.multiConst(3),settings.multiConst(4),settings.multiConst(5));

PARAMS = pppParams;
PARAMS.IMU_ARM = IMU_ARM;
PARAMS.measUse.excludeThresh = PARAMS.measUse.excludeThresh;


PARAMS.states.RX_DCB_GLO = true;
PARAMS.states.RX_DCB_GPS = true;
PARAMS.states.MP_CODE = false;
PARAMS.states.MP_CARR = false;

PARAMS.Q.AMB = 0.005;

PARAMS.solSep.nMaxSubset = 100;


% parse the inertial measurements
if ~exist('gnssImuFull','var')
    
    [~,~,extObs] = fileparts(filenameObs);
    
    switch extObs
        case {'.ASC','.GPS'}
            
            gnssImuFull = readSPANdata(filenameObs,'RANGEA',RANGEA,...
                'RAWIMUA',true,'constellations',constellations,...
                'correctHalfCycle',true);
            
            if ~isempty(filenameGnss)
                [obsStruc, epochs, time, week, date, rxPos0, interval, ...
                    antOffset, antMod, rxmod] = load_RINEX_obs_fullObs(filenameGnss,constellations);
                
                obsGnssRaw.meas      = obsStruc;
                obsGnssRaw.PRN       = constellations.PRN;
                obsGnssRaw.constInds = constellations.constInds;
                obsGnssRaw.epochs    = time;
            else
                obsGnssRaw.meas      = gnssImuFull.RANGEA.obsData;
                obsGnssRaw.PRN       = constellations.PRN;
                obsGnssRaw.constInds = constellations.constInds;
                obsGnssRaw.epochs   = gnssImuFull.RANGEA.epochs;
            end
            
            tLock = gnssImuFull.RANGEA.tLock;
            tLock.epochs = obsGnssRaw.epochs;
            
            obsGnssRaw.tLock = tLock;
            
            if ~isempty(gnssImuFull.RAWIMUA.headerWeek)
                obsImu = gnssImuFull.RAWIMUA;
            else
                obsImu = gnssImuFull.RAWIMUSXA;
            end
            
            obsImu.acc  = obsImu.acc*accScale;
            obsImu.gyro = obsImu.gyro*gyroScale;
            
            obsImu.acc = [sign(imuOrder(1))*obsImu.acc(:,abs(imuOrder(1))) ...
                sign(imuOrder(2))*obsImu.acc(:,abs(imuOrder(2))) ...
                sign(imuOrder(3))*obsImu.acc(:,abs(imuOrder(3)))];
            
            obsImu.gyro = [sign(imuOrder(1))*obsImu.gyro(:,abs(imuOrder(1))) ...
                sign(imuOrder(2))*obsImu.gyro(:,abs(imuOrder(2))) ...
                sign(imuOrder(3))*obsImu.gyro(:,abs(imuOrder(3)))];
            
            antMod = [];
            antOffset = [0 0 0]';
            
        case '.mat'
            datai = load(filenameObs);
            
            % Make sure to only keep the constellations of interest
            indsKeep = ismember([datai.constellations.PRN' datai.constellations.constInds'],...
                [constellations.PRN' constellations.constInds'],'rows');
            measTypes = fields(datai.obsStruc);
            for idx = 1:length(measTypes)
                obsStruc.(measTypes{idx}) = datai.obsStruc.(measTypes{idx})(indsKeep,:);
            end
            
            obsGnssRaw.meas      = obsStruc;
            obsGnssRaw.PRN       = constellations.PRN(indsKeep);
            obsGnssRaw.constInds = constellations.constInds(indsKeep);
            obsGnssRaw.epochs    = datai.time';
            
            obsImu = [];
    end
    
end

%% pre-process the observation data- only keep what we want to use and combine to make dual frequency measurements

obsInds    = repmat([1 2 3 4],1,3); % 1 = code, 2 = carrier, 3 = snr, 4 = doppler
signalInds = kron(1:3,ones(1,4));
%           1                               2                               3                               4                                5
obsDes  = {{'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2S'} {'L2S'} {'S2S'} {'D2S'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}  };

% obsDes  = {{'C1W'} {'L1W'} {'S1W'} {'D1W'} {'C2B'} {'L2B'} {'S2B'} {'D2B'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%     {'C1P'} {'L1P'} {'S1P'} {'D1P'} {'C2B'} {'L2B'} {'S2B'} {'D2B'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%     {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%     {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%     {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}  };



% signal pairs to include as iono-free combinations
ifPairs = [1 3;
    1 2];

% load precise orbit and clock
[obsGnss, dcbCorr0] = preprocessGnssObs(obsGnssRaw,obsInds,signalInds,obsDes,ifPairs,settings,PARAMS);

% add dual frequency to the dcb corrections
dcbCorr = nan(10,size(dcbCorr0,2));
dcbCorr0(dcbCorr0 == 0) = nan;
dcbCorr(1:6,:) = dcbCorr0;
dcbCorr(7,:)   = (dcbCorr0(1,:).*obsGnss.freqs(1,:).^2-dcbCorr0(3,:).*obsGnss.freqs(3,:).^2)./(obsGnss.freqs(1,:).^2-obsGnss.freqs(3,:).^2);
dcbCorr(9,:)   = (dcbCorr0(1,:).*obsGnss.freqs(1,:).^2-dcbCorr0(5,:).*obsGnss.freqs(5,:).^2)./(obsGnss.freqs(1,:).^2-obsGnss.freqs(5,:).^2);

%% load precise orbit and clock data
if ~exist('corrData','var')
    epoch0 = obsGnss.epochs(1);
    [doy,year] = jd2doy(epochs2jd(epoch0));
    doy = floor(doy);
    
    corrData = svOrbitClock;
    % load precise ephemeris
    corrData.initPEph(year,doy,settings);
    % load precise clock data
    corrData.initPClock(year,doy,settings);
    
    % load tec map
    corrData.initIonoData(year,doy,settings);
    
    % load antenna phase center data
    filenameAtx = [cd '\Utility Functions\igs14.atx'];
    
    corrData.initAtxData(filenameAtx);
end


%% setup the exact measurements to be passed in
indsGnssMeas = find(obsGnss.epochs >= epochStart & obsGnss.epochs < epochEnd);
% indsGnssMeas = indsGnssMeas(1:10:end);
obsGnssi = obsGnss;
obsGnssi.epochs      = obsGnssi.epochs(indsGnssMeas);
obsGnssi.range.obs   = obsGnssi.range.obs(:,:,indsGnssMeas);
obsGnssi.doppler.obs = obsGnssi.doppler.obs(:,:,indsGnssMeas);
obsGnssi.snr.obs     = obsGnssi.snr.obs(:,:,indsGnssMeas);

if ~isempty(obsImu) && useImu
    indsImuMeas = find(obsImu.epochs >= min(obsGnssi.epochs) & obsImu.epochs < max(obsGnssi.epochs));
    imuFields = fields(obsImu);
    obsImui = obsImu;
    for fieldi = imuFields'
        obsImui.(fieldi{1}) = obsImui.(fieldi{1})(indsImuMeas,:);
    end
else
    obsImui = [];
    obsImui.epochs = [];
end

PARAMS.c = settings.c;

if ~useImu
    PARAMS.IMU_ARM = PARAMS.IMU_ARM*0;
end

obsGnssi.range.obs(:,39,:) = 0;
obsGnssi.doppler.obs(:,39,:) = 0;
obsGnssi.snr.obs(:,39,:) = 0;

%% run the estimator!
[outStruc,outSs] = pppImuSs(obsGnssi,obsImui,corrData,PARAMS,settings);

%% Output plots
outStruc.plotResidSummary;

outStruc.plotResids;

outStruc.plotAgainstTruth(truthFile,PARAMS);

outStruc.midRunPlot;

outStruc.plotRemoved;

outSs.plotSs(truthFile);


































