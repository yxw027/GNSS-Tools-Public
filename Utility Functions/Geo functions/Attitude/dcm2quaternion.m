function q = dcm2quaternion(dcm)

% simple function to translate from direction cosines matrix to quaternion

c11 = dcm(1,1);
c22 = dcm(2,2);
c33 = dcm(3,3);
c12 = dcm(1,2);
c13 = dcm(1,3);
c21 = dcm(2,1);
c23 = dcm(2,3);
c31 = dcm(3,1);
c32 = dcm(3,2);

% exact translation is from Inertial Navigation Systems with Geodetic
% Applications by Jekeli, page 18

q1 = 1/2*sqrt(1+c11+c22+c33);
q2 = 1/(4*q1)*(c23-c32);
q3 = 1/(4*q1)*(c31-c13);
q4 = 1/(4*q1)*(c12-c21);

q = [q1 q2 q3 q4]';


end