function dop = dopplerModel(vel,svVel,clockDrift,A,gyroMeas,att,IMU_ARM,Omega_ie)

%  vel = velRx-R_b_e*crossProdMatrix(gyroMeas)*PARAMS.IMU_ARM + ...
%         Omega_ie*R_b_e*PARAMS.IMU_ARM;
% R_b_e = euler2dcm123(att);

% velRx = vel+R_b_e*crossProdMatrix(gyroMeas)*IMU_ARM - ...
%         Omega_ie*R_b_e*IMU_ARM;

velRx = velModel(vel,svVel,A,gyroMeas,att,IMU_ARM,Omega_ie);

dVel = velRx-svVel;

dop = -dot(dVel,-A)-clockDrift;

end