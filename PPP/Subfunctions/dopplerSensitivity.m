function [hVel,hAtt,hBias] = dopplerSensitivity(vel,svVel,clockDrift,A,gyroMeas,R_b_e,IMU_ARM,Omega_ie)


% 
dx = 1e-8;

% velocity sensitivity
hVel = nan(3,1);
for idx = 1:3
    veli = vel;
    veli(idx) = vel(idx)+0.5*dx;
    
    doppP = dopplerModel(veli,svVel,clockDrift,A,gyroMeas,R_b_e,IMU_ARM,Omega_ie);
    
    veli = vel;
    veli(idx) = vel(idx)-0.5*dx;
    
    doppM = dopplerModel(veli,svVel,clockDrift,A,gyroMeas,R_b_e,IMU_ARM,Omega_ie);
    
    
    hVel(idx,:) = (doppP-doppM)./dx;
    
end

% attitude sensitivity
hAtt = nan(3,1);
for idx = 1:3
%     atti = att;
%     atti(idx) = atti(idx)+0.5*dx;

    wi = zeros(3,1);
    wi(idx) = 0.5*dx;
    R_b_ei = (eye(3) + crossProdMatrix(wi)) * R_b_e;


    doppP = dopplerModel(vel,svVel,clockDrift,A,gyroMeas,R_b_ei,IMU_ARM,Omega_ie);
    
    wi = zeros(3,1);
    wi(idx) = -0.5*dx;
    R_b_ei = (eye(3) + crossProdMatrix(wi)) * R_b_e;
    
    doppM = dopplerModel(vel,svVel,clockDrift,A,gyroMeas,R_b_ei,IMU_ARM,Omega_ie);
    
    
    hAtt(idx,:) = (doppP-doppM)./dx;
    
end
% bias sensitivity
hBias = nan(3,1);
for idx = 1:3
%     atti = att;
%     atti(idx) = atti(idx)+0.5*dx;

    gyroMeasi = gyroMeas;
    gyroMeasi(idx) = gyroMeas(idx)+0.5*dx;
%     R_b_ei = (eye(3) + crossProdMatrix(wi)) * R_b_e;


    doppP = dopplerModel(vel,svVel,clockDrift,A,gyroMeasi,R_b_e,IMU_ARM,Omega_ie);
    
    gyroMeasi = gyroMeas;
    gyroMeasi(idx) = gyroMeas(idx)-0.5*dx;
    
    doppM = dopplerModel(vel,svVel,clockDrift,A,gyroMeasi,R_b_e,IMU_ARM,Omega_ie);
    
    
    hBias(idx,:) = (doppP-doppM)./dx;
    
end



end