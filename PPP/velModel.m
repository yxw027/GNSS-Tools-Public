function velRx = velModel(vel,svVel,A,gyroMeas,R_b_e,IMU_ARM,Omega_ie)

% R_b_e = euler2dcm(att);

velRx = vel+R_b_e*crossProdMatrix(gyroMeas)*IMU_ARM - ...
        Omega_ie*R_b_e*IMU_ARM;



end