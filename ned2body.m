function [var_xb, var_yb, var_zb] = ned2body(var_xn, var_yn, var_zn, phi, theta, psi )
% Phi = Roll, Theta = pitch, Psi = yaw
% Will rotate var_n (NED) vector to var_b frame (Body)


% Rotation matrix 3-2-1
rotationMat(1,1) = cos(theta)*cos(psi);
rotationMat(1,2) = cos(theta)*sin(psi);
rotationMat(1,3) = -sin(theta);
rotationMat(2,1) =  sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
rotationMat(2,2) = sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi);
rotationMat(2,3) = sin(phi)*cos(theta);
rotationMat(3,1) = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
rotationMat(3,2) = cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
rotationMat(3,3) = cos(phi)*cos(theta);

var_n = [var_xn; var_yn; var_zn];

var_b = rotationMat * var_n;

var_xb = var_b(1);
var_yb = var_b(2);
var_zb = var_b(3);