function initialize(obj,epoch,PARAMS,varargin)

p = inputParser;

p.addParameter('pos',[]);
p.addParameter('vel',[]);
p.addParameter('R_b_e',[]);
p.addParameter('imuBiasStates',[]);
p.addParameter('constsUnique',1);
p.addParameter('clockBias',0);
p.addParameter('clockDrift',0);
p.addParameter('PRN',[]);
p.addParameter('constInds',[]);
p.addParameter('range',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;

pos           = res.pos;
vel           = res.vel;
R_b_e         = res.R_b_e;
imuBiasStates = res.imuBiasStates;
constsUnique  = res.constsUnique;
clockDrift    = res.clockDrift;
clockBias     = res.clockBias;
PRN           = res.PRN(:);
constInds     = res.constInds(:);
rangeStruc    = res.range;

%% Populate various fields based on what we have
obj.pos = pos;

obj.vel = vel;

obj.R_b_e = R_b_e;

obj.imuBiasStates = imuBiasStates;

obj.clockBias  = clockBias;
obj.clockDrift = clockDrift;

obj.epochLastInertialUpdate = epoch;
obj.epochLastGnssUpdate     = epoch;

%% Setup state, covariance, and state index mapping
obj.INDS_STATE = [];
obj.INDS_STATE.ATTITUDE = 1:3;
obj.INDS_STATE.VEL = 4:6;
obj.INDS_STATE.POS = 7:9;
obj.INDS_STATE.ACC_BIAS = 10:12;
obj.INDS_STATE.W_BIAS = 13:15;

% Add the correct number of clock biases and drifts
lastInd = 15;
obj.INDS_STATE.CLOCK_BIAS = (lastInd+1):(lastInd+length(constsUnique));
obj.INDS_STATE.CLOCK_BIAS_CONSTS = constsUnique;
lastInd = lastInd+length(constsUnique);
obj.INDS_STATE.CLOCK_DRIFT = (lastInd+1):(lastInd+length(constsUnique));
obj.INDS_STATE.CLOCK_DRIFT_CONSTS = constsUnique;
lastInd = lastInd+length(constsUnique);

if PARAMS.states.trop
    obj.INDS_STATE.TROP = lastInd+1;
    
    lastInd = lastInd+1;
else
   obj.INDS_STATE.TROP = []; 
end

if PARAMS.states.RX_DCB
    % For each constellation and signal (excluding one reference), we need
    % a receiver differential code bias state
    if isempty(rangeStruc)
        error('Ranging measuremnets required for RX DCB state initialization')
    end
    
    sigConstInds = unique([rangeStruc.sig(:) rangeStruc.constInds(:)],'rows');
    constUn = unique(rangeStruc.constInds(:));
    for idx = 1:length(constUn)
        if constUn(idx) == 1 && PARAMS.states.RX_DCB_GPS
            % If using separate DCBs for each GLONASS satellite, no need
            % for a signal specific one
            continue;
        end
        
        if constUn(idx) == 2 && PARAMS.states.RX_DCB_GLO
            % If using separate DCBs for each GLONASS satellite, no need
            % for a signal specific one
            continue;
        end
        % Remove a reference signal from each constellation- just using the largest one for now.
        % this should change based on whereever measurement masking happens
        refSigi = max(sigConstInds(sigConstInds(:,2) == constUn(idx),1));
        
        sigConstInds(ismember(sigConstInds,[refSigi constUn(idx)],'rows'),:) = [];
    end
    
    % Add the state indices
    nStateDcb = size(sigConstInds,1);
    
    obj.INDS_STATE.RX_DCB.INDS = lastInd+[1:nStateDcb];
    obj.INDS_STATE.RX_DCB.sig  = sigConstInds(:,1)';
    obj.INDS_STATE.RX_DCB.constInds = sigConstInds(:,2)';
    
    lastInd = lastInd+nStateDcb;
else
    obj.INDS_STATE.RX_DCB.INDS = [];
    obj.INDS_STATE.RX_DCB.sig  = [];
    obj.INDS_STATE.RX_DCB.constInds = [];
end

obj.INDS_STATE.FLEX_STATE_MIN = lastInd+1;
obj.INDS_STATE.FLEX_STATES = zeros(0,1);
obj.INDS_STATE.FLEX_STATES_INFO = zeros(0,4); % PRN | CONST | TYPE | SIG INDICATOR 

% Initialize covariance
cov = zeros(lastInd);
cov(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE)       = eye(3) * PARAMS.SIGMA0.ATTITUDE^2;
cov(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL)                 = eye(3) * PARAMS.SIGMA0.VEL^2;
cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS)                 = eye(3) * PARAMS.SIGMA0.POS^2;
cov(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS)       = eye(3) * PARAMS.SIGMA0.ACC_BIAS^2;
cov(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS)           = eye(3) * PARAMS.SIGMA0.W_BIAS^2;
cov(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_BIAS)   = eye(length(obj.INDS_STATE.CLOCK_BIAS))*100^2;
cov(obj.INDS_STATE.CLOCK_DRIFT,obj.INDS_STATE.CLOCK_DRIFT) = eye(length(obj.INDS_STATE.CLOCK_DRIFT))*100^2;
cov(obj.INDS_STATE.TROP,obj.INDS_STATE.TROP)               = PARAMS.SIGMA0.TROP^2;
cov(obj.INDS_STATE.RX_DCB.INDS,obj.INDS_STATE.RX_DCB.INDS) = eye(length(obj.INDS_STATE.RX_DCB.INDS))*PARAMS.SIGMA0.RX_DCB^2;

obj.cov = cov;

% Initialize state
state = zeros(lastInd,1);

obj.state = state;

%% Initialize phase windup if we have PRN and constInds
if ~isempty(PRN) && ~isempty(constInds)
    nSv = length(PRN);
    obj.phWind.phaseOffset = zeros(size(PRN));
    obj.phWind.PrnConstInd = [PRN constInds];
end

end















