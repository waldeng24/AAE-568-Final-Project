% [xPred, P] = kalmanFilter(F,x,B,u,P,H,y,Q,R) returns state estimate, x 
% and state covariance, P for dynamic system:
%           x[k+1] = f(x[k], u[k]) + w[k]
%           y[k]   = h(x[k]) + v[k]
% where w ~ N(0,Q) meaning w is Gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is Gaussian noise with covariance R
% Inputs:   
%    F: state dynamics matrix
%    x: current state estimate
%    B: commands matrix
%    u: actual commands
%    P: current state estimate covariance
%    H: measurement matrix
%    y: observation / measurement
%    Q: process noise covariance 
%    R: measurement noise covariance
% Output:
%    xPred: updated state estimate
%    P:    updated state covariance
%
function [xPred, P] = kalmanFilter(F,x,B,u,P,H,y,Q,R)

% Prediction portion
xHat = F*x + B*u;
Ppriori = F*P*transpose(F) + Q;

% Update step
z = y-H*x;
S = H*Ppriori*transpose(H) + R;
K = Ppriori*transpose(H)*inv(S);

% Output values
xPred = xHat + K*z;
P = (eye(9) - K*H)*Ppriori;


end