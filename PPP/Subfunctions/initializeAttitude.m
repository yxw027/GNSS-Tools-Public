function [attitude0, dcm] = initializeAttitude(usrVelEnu,PARAMS)

% Initialize attitude given the initial estimate of position and  velocity 
% (might need something else to do this though)
% This also assumes that z axis is pointing upward

% Form of the output attitude is dependent on PARAM.state.attitudeMode

%% first just building attitude unit vectors

% matching what novatel does
% x back
% y right
% z up

% z = usrPos0./norm(usrPos0);
% y1 = usrVel0./norm(usrVel0);
% x = cross(y1,z);
% y = cross(z,x);

xi = -usrVelEnu./norm(usrVelEnu);

z0i = [0 0 1]';
yi = -cross(xi,z0i);
zi = cross(xi,yi);

dcm = [xi';
       yi'; 
       zi']';

%% translate to desired attitude output (could just be DCM)
switch PARAMS.states.attitudeMode
    case 'QUATERNION'
        attitude0 = dcm2quaternion(dcm);
        
    case 'EULER'
        attitude0 = dcm2euler123(dcm);
                
    case 'DCM'
        attitude0 = dcm;
end

end