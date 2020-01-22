function [outStruc,outSolSep] = pppImuSs(obsGnss,obsImu,corrData,PARAMS,settings,varargin)

%% Parse inputs
p = inputParser;

p.addParameter('truthPos',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
truthPos   = res.truthPos;  % Truth position input

%% Sort all of the GNSS and IMU measurements
% 1 = GNSS, 2 = IMU
%                      EPOCH      |        GNSS/IMU            |  INDEX WITHIN TYPE
obsInfo = sortrows([obsGnss.epochs 1*ones(size(obsGnss.epochs)) [1:size(obsGnss.epochs)]';
    obsImu.epochs  2*ones(size(obsImu.epochs))  [1:size(obsImu.epochs)]'],1);

indStart = min(find(obsInfo(:,2) == 1));
obsInfo = obsInfo(indStart:end,:);

nEpochs = size(obsInfo,1);

%% strip the first GNSS measurement and do a least squares solution to initialize
obsi = stripMeas({obsGnss obsImu},obsInfo(1,2),obsInfo(1,3),PARAMS);

[state,dstate] = lsSolGnss(obsi,corrData,PARAMS,settings);

pos0 = state(1:3);
vel0 = dstate(1:3);
epoch0 = obsi.epochs;
clockBias0 = state(4:end);
clockDrift0 = dstate(4:end);

%%
% Initialize attitude - euler angles with body frame
xi  = -vel0./norm(vel0);
z0i = pos0./norm(pos0);
yi  = -cross(xi,z0i);
zi = cross(xi,yi);
R_b_e = [xi yi zi];

pos0 = pos0-R_b_e*PARAMS.IMU_ARM;

% Initialize bias estimates
imuBiasStates = zeros(6,1); % accelerometer(1:3) and gyro(1:3)
imuBiasStates(1:3) = [-0.3 0.1 -0.4];

user = inertialFilter;
user.initialize(epoch0,PARAMS,'pos',pos0,'vel',vel0,'imuBiasState',imuBiasStates,...
    'R_b_e',R_b_e,'constsUnique',unique(obsGnss.constInds),'clockBias',clockBias0,...
    'clockDrift',clockDrift0,'PRN',obsi.PRN,'constInds',obsi.constInds,...
    'range',obsGnss.range);

% Initialize output structure- this is really just for the all in view.
outStruc = inertialSave(obsInfo(:,1),length(clockBias0),obsGnss,...
    obsGnss.epochs,user.INDS_STATE,obsInfo(:,2),'gnssEpochsOnly',false);
epochLastPlot = obsInfo(1,1);

outStruc.saveState(user,PARAMS,'tdx', 1)

% Initialize the solution separation
solSep0 = solSep;
solSep0(1).ekf        = user;
solSep0(1).aivFlag    = true;
solSep0(1).measOut    = [];
solSep0(1).subsetType = 0;
solSep0(1).pSat       = 1;
solSep0(1).ID         = 0;

ss = solSepWrap(PARAMS);
ss.ssList(1) = solSep0;

% Initialize the output data for the solution separation
outSolSep = solSepSave(obsInfo(obsInfo(:,2) == 1));

% Initialize the waitbar
runTimeStart = tic;
pctDone = 0;
h = waitbar(0,'0 Percent Complete');

PARAMS_subsets = copy(PARAMS);
PARAMS_subsets.measUse.excludeThresh = PARAMS_subsets.measUse.excludeThresh*Inf;

%% Run the loop
for tdx = 2:nEpochs
    % Pull the measurement from the full list
    obsi = stripMeas({obsGnss obsImu},obsInfo(tdx,2),obsInfo(tdx,3),PARAMS);
    
    epochi = obsInfo(tdx,1);
    
    measType = obsInfo(tdx,2);
    
    updated = true;
    
    switch measType
        case 1
            % GNSS measurements
            manageSubsets(ss,obsi,PARAMS);
            
            obs0 = obsi;
            
            % Manage all states- they are kept the same here
            manageStatesMulti([ss.ssList.ekf],epochi,obsi,PARAMS,outStruc);
            
            if exist('accMeasi','var')
                for idx = 1:length(ss.ssList)
                    % Propagate position, velocity, attitude given acc and gyro
                    % measurements
                    ss.ssList(idx).ekf.inertialNavEquationsEcef(epochi,accMeasi,gyroMeasi);
                end
            end
            
            aivInfoShare = [];
            for idx = 1:length(ss.ssList)
                % strip measurements from each subset
                %                 obsi = solSepMeasMask(obs0,ss.ssList(idx).measOut);
                
                obsi = obs0;
                
                if ~any(obsi.range.obs(:))
                    updated = false;
                    continue;
                end
                
                % Do the tightly coupled update
                if idx == 1
                    [aivInfoShare, measMatRemoved,measMatLow] = ss.ssList(idx).ekf.tcUpdate(epochi,obsi,corrData,...
                        PARAMS,settings,outStruc);
                    %                     obs0 = aivRemovedMask(obs0,measMatRemoved);
                else
                    ss.ssList(idx).ekf.tcUpdate(epochi,obsi,corrData,...
                        PARAMS_subsets,settings,[],aivInfoShare,ss.ssList(idx).measOut);
                end
            end
            
            if updated
                % Compute protection levels
                pl = ss.computePL(PARAMS);
                
                % Save info
                outSolSep.saveNewValue(epochi,ss,pl,PARAMS);
            end
        case 2
            % IMU
            % Pull measurements from this epoch
            accMeas0 = obsi.acc'*(9.83);
            gyroMeas0 = obsi.gyro';
            
            for idx = 1:length(ss.ssList)
                % Correct measurements given bias estimates
                accMeasi  = accMeas0-ss.ssList(idx).ekf.imuBiasStates(1:3);
                gyroMeasi = gyroMeas0-ss.ssList(idx).ekf.imuBiasStates(4:6);
                
                % Propagate position, velocity, attitude given acc and gyro
                % measurements
                ss.ssList(idx).ekf.inertialNavEquationsEcef(epochi,accMeasi,gyroMeasi);
            end
            
            pl = [];
    end
    
    if updated
        outStruc.saveState(ss.ssList(1).ekf,PARAMS,'tdx',tdx,'pl',pl)
    end
    % Update the waitbar
    if mod(floor(tdx/nEpochs*100),1) == 0 && floor(tdx/nEpochs*100) > pctDone
        tElapsed = toc(runTimeStart);
        tRemaining = tElapsed*(nEpochs-tdx)./tdx;
        pctDone =  floor(tdx/nEpochs*100);
        waitbar(pctDone/100,h,[num2str(pctDone) '% Complete, ' num2str(tRemaining/60,'%5.2f') ' Minutes Remaining'])
    end
end

close(h);
outStruc.closeSaveState


end














