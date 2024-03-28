%% DEBUG
clc; clear all; close all;
run("createTrajectory.m");

%% Main code
% Conversion constants
d2r = pi/180;
r2d = 180/pi;
gToMps = 9.80665;
toMilli = 1000;
fromMilli = 0.001;
ft2m = 0.305;
hr2sec = 3600;

% Configuration 1
timeStep   = 0.1; % Seconds
updateTime = 0;
endSim     = false; % Termination flag
satUpdateFreq = 0.1; % Sec
imuUpdateFreq = 0.2; % Sec
satUpdateTime = 0.1; % Next time to update
imuUpdateTime = 0.2; % Next time to update

% % Configuration 2
% timeStep   = 0.1; % Seconds
% updateTime = 0;
% endSim     = false; % Termination flag
% satUpdateFreq = 2; % Sec
% imuUpdateFreq = 4; % Sec
% satUpdateTime = 0.5; % Next time to update
% imuUpdateTime = 1; % Next time to update

% Initialization
iter = 1;
t = 0;
eulerBody(iter,:) = [0 0 0]; % roll, pitch, yaw - radians
eulerRateBody(iter,:) = [0 0 0];
posNED(iter,:) = waypointsNED(1,:);
posLLA(iter,:)  = waypointsLLA(1,:);
posBody(iter,:) = transpose(bodyToNED(eulerBody(iter,:)))*transpose(posNED(iter,:));
wpTo   = 2;
wpFrom = 1;
unitDirHeadingNED = (waypointsNED(wpTo,:) - waypointsNED(wpFrom,:)) ./ norm(waypointsNED(wpTo,:) - waypointsNED(wpFrom,:));
velNED(iter,:)    = unitDirHeadingNED .* spec.speedMetersPerSec;
velBody(iter,:)   = transpose(bodyToNED(eulerBody(iter,:)))*transpose(velNED(iter,:));
accelBody(iter,:) = [0 0 9.80665];
accelNED(iter,:) = bodyToNED(eulerBody(iter,:))*transpose(accelBody(iter,:)); % Steady level flight

% Kalman Characteristics
q = 3; % Process uncertainty
Q = q^2*eye(6);
ripH = 1   ; % IMU Measurement uncertainty (horizontal)
ripV = 9   ; % IMU Measurement uncertainty (vertcal)
rgpH = 1.6 ; % GPS Measurement uncertainty (horizontal)
rgpV = 3   ; % GPS measurement uncertainty (vertical)
rgv  = 0.1 ; % GPS Velocity uncertainty (all directions)
Rimu = [ripH^2 0 0 0 0 0;
        0 ripH^2 0 0 0 0;
        0 0 ripV^2 0 0 0;
        zeros(3,3) 0.1^2*eye(3)];
Rgps = [rgpH^2 0 0 0 0 0;
        0 rgpH^2 0 0 0 0;
        0 0 rgpV^2 0 0 0;
        zeros(3,3) rgv^2*eye(3)];
P = Rimu.*Rgps; % Initial covariance

% Kalman filter outputs and track for performance analysis
filterEstX(iter,:) = [posNED(iter,:)'; velNED(iter,:)'];
filterCovPosX(iter,:) = P(1,1);
filterCovPosY(iter,:) = P(2,2);
filterCovPosZ(iter,:) = P(3,3);
filterCovVelX(iter,:) = P(4,4);
filterCovVelY(iter,:) = P(5,5);
filterCovVelZ(iter,:) = P(6,6);
x = [posNED(iter,:)'; velNED(iter,:)']; % Filter state

% Create IMU and GPS sensors
gps = gpsSensor;
paramsAccel = accelparams;
paramsAccel.MeasurementRange = 30.0 * gToMps;
paramsAccel.ConstantBias = repelem(3.0 * fromMilli * gToMps, 3);
paramsAccel.RandomWalk = repelem(0.132 * ft2m, 3);
paramsGyro = gyroparams;
paramsGyro.MeasurementRange = 1000*d2r;
paramsGyro.ConstantBias = repelem(2 * d2r * hr2sec, 3);
paramsGyro.RandomWalk = repelem(0.125 * d2r, 3);
imu = imuSensor(...
    'accel-gyro', ...
    'Accelerometer', paramsAccel, ...
    'Gyroscope', paramsGyro);

% Create true flight path taken
while ((endSim == false) && (iter < 100000))
    % Update time step and iteration
    iter = iter + 1;
    t = t + timeStep;

    % Simulate measurement from the sensor - calculate velocity and
    % position readings in the body and rotate to NED
    if ( t >= imuUpdateTime )
        %[accelReading, eulerReading] = simulateIMU(accelBody(iter-1,:), eulerBody(iter-1,:));
        [accelReading, eulerReading] = imu( accelBody(iter-1,:), eulerBody(iter-1,:) );
        eulerReading = eulerReading .* r2d;

        velReadBody = velBody(iter-1,:) + accelReading .* timeStep;
        posReadBody = posBody(iter-1,:) + velReadBody .*timeStep;
        velReadNED = bodyToNED(eulerBody(iter-1,:))*transpose(velReadBody);
        posReadNED = bodyToNED(eulerBody(iter-1,:))*transpose(posReadBody);
        eulerRateRead = eulerRateBody(iter-1,:) + eulerReading .* timeStep;
        eulerAngRead = eulerBody(iter-1,:) + eulerRateRead .* timeStep;

        % Kalman update
        y = [posReadNED; velReadNED];
        u = [0; 0; 0; velNED(iter-1,:)'];
        F = [ 1 0 0 t-updateTime 0 0;
            0 1 0 0 t-updateTime 0;
            0 0 1 0 0 t-updateTime;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 1 ];
        B = [0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 t-updateTime 0 0;
            0 0 0 0 t-updateTime 0;
            0 0 0 0 0 t-updateTime];
        H = [1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 1];

        [x, P] = kalmanFilter(F, x, B, u, P, H, y, Q, Rimu);
        updateTime = t;
        imuUpdateTime = imuUpdateTime + imuUpdateFreq;
    end

    % Simulate measurement from the GPS signal
    % add measurments of X and Y with some error and noise
    if ( t >= satUpdateTime )
        %[gpsPosRead, gpsVelRead] = simulateGPS2(posNED(iter-1,:), velNED(iter-1,:), posLLA(1,:));
        [gpsPosReadLLA, gpsVelReadMpS] = gps(posNED(iter-1,:), velNED(iter-1,:));
        [gpsPosRead(1) gpsPosRead(2) gpsPosRead(3)] = geodetic2ned(gpsPosReadLLA(1), gpsPosReadLLA(2), gpsPosReadLLA(3), 0, 0, 0, wgs84Ellipsoid);
        y = [gpsPosRead'; gpsVelReadMpS'];
        u = [0; 0; 0; velNED(iter-1,:)'];
        F = [ 1 0 0 t-updateTime 0 0;
            0 1 0 0 t-updateTime 0;
            0 0 1 0 0 t-updateTime;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 1 ];
        B = [0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 t-updateTime 0 0;
            0 0 0 0 t-updateTime 0;
            0 0 0 0 0 t-updateTime];
        H = [1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 1];

        [x, P] = kalmanFilter(F, x, B, u, P, H, y, Q, Rgps);
        updateTime = t;
        satUpdateTime = satUpdateTime + satUpdateFreq;
    end

    % If a time step goes by and no measurement was made, propagate the
    % state and the covariance
    if ( t > updateTime )
        x(1:3) = x(1:3)' + x(4:6)'.*timeStep;
        P = P + Q;
    end

    % Log the filter state estimate and covariance matrix
    filterEstX(iter,:) = x';
    filterCovPosX(iter,:) = P(1,1);
    filterCovPosY(iter,:) = P(2,2);
    filterCovPosZ(iter,:) = P(3,3);
    filterCovVelX(iter,:) = P(4,4);
    filterCovVelY(iter,:) = P(5,5);
    filterCovVelZ(iter,:) = P(6,6);

    % Propagate the position based on the velocity from before
    posNED(iter,:) = posNED(iter-1,:) + velNED(iter-1,:) .* timeStep;
    [lat, lon, h] = ned2geodetic(posNED(iter,1), posNED(iter,2), posNED(iter,3), ...
        waypointsLLA(1,1), waypointsLLA(1,2), waypointsLLA(1,3), wgs84Ellipsoid);
    posLLA(iter,:) = [lat, lon, h];
    posBody(iter,:) =  transpose(bodyToNED(eulerBody(iter-1,:)))*transpose(posNED(iter-1,:));
    velNED(iter,:)   = unitDirHeadingNED .* spec.speedMetersPerSec;
    velBody(iter,:)   = transpose(bodyToNED(eulerBody(iter-1,:)))*transpose(velNED(iter-1,:));
    % switch wpFrom
    %     case 1
    eulerBody(iter,:) = [0 0 0];
    %     case 2
    %          eulerBody(iter,:) = [0 0 90];
    %     case 3
    %          eulerBody(iter,:) = [0 0 180];
    %     case 4
    %          eulerBody(iter,:) = [0 0 270];
    % end
    eulerRateBody(iter,:) = [0 0 0];
    accelBody(iter,:) = [0 0 9.80665];
    accelNED(iter,:) = bodyToNED(eulerBody(iter,:))*transpose(accelBody(iter,:)); % Steady level flight

    % If you've passed the waypoint of interest, head towards the next one
    display( num2str(distance(posLLA(iter-1,1), posLLA(iter-1,2), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid)) );
    if ( distance(posLLA(iter-1,1), posLLA(iter-1,2), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid) < (spec.speedMetersPerSec * timeStep) )
        wpFrom = wpTo;
        wpTo = min(wpTo + 1, height(waypointsNED));
    end

    % Find direction of travel for the next iteration
    unitDirHeadingNED = (waypointsNED(wpTo,:) - posNED(iter,:)) ./ norm(waypointsNED(wpTo,:) - posNED(iter,:));
    velNED(iter,:) = unitDirHeadingNED .* spec.speedMetersPerSec .* timeStep;

    % If you passed the final waypoint, terminate simulation
    if ( wpTo == height(waypointsNED) && ...
            distance(posLLA(iter-1,1), posLLA(iter-1,2), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid) < (spec.speedMetersPerSec * timeStep) )
        endSim = true;
    end

end

%% DEBUG
outData = 1:100:iter;
timeVec = timeStep.*[0:1:iter];
timeVec = timeVec(outData);

% True position waypoint based
figure; hold on; grid on;
plot3(waypointsNED(:,1), waypointsNED(:,2), -waypointsNED(:,3), 'b')
scatter3(waypointsNED(1,1), waypointsNED(1,2), -waypointsNED(1,3), 'ko', 'filled')
scatter3(waypointsNED(2:end-1,1), waypointsNED(2:end-1,2), -waypointsNED(2:end-1,3), 'bo')
scatter3(posNED(outData,1), posNED(outData,2), -posNED(outData,3), 'r.');
view(3)
title('NED Position Visualization');
xlabel('xNorth (meters)');
ylabel('Negative yEast (meters)');
zlabel('Negative zDown (meters)');

% Filtered position for tracking - 3D
figure; hold on; grid on;
scatter3(posBody(outData,1), posBody(outData,2), posBody(outData,3), 'r.');
plot3(filterEstX(outData,1), filterEstX(outData,2), filterEstX(outData,3), 'b')
view(3)
legend('True Body Pos', 'KF Estimation');
title('Body Position Filter Estimate Comparison');
xlabel('xNorth (meters)');
ylabel('yEast (meters)');
zlabel('Negative zDown (meters)');

% Filtered position for tracking - 2D
figure; hold on; grid on;
scatter(posBody(outData,1), posBody(outData,2), 'r.');
plot(filterEstX(outData,1), filterEstX(outData,2), 'b')
legend('True Body Pos', 'KF Estimation');
title('Body Position Filter Estimate Comparison');
xlabel('xNorth (meters)');
ylabel('yEast (meters)');

% Each individual position compared to it's error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('North Position Comparison')
plot(timeVec, posNED(outData,1), 'b');
plot(timeVec, filterEstX(outData,1), 'r');
plot(timeVec, filterEstX(outData,1) + sqrt(filterCovPosX(outData)), 'k');
plot(timeVec, filterEstX(outData,1) - sqrt(filterCovPosX(outData)), 'k');
xlabel('Time (s)');
ylabel('X Pos North (m)');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(1) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosX(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(filterCovPosX), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma(1,1)');
ax(2) = gca;

figure;
subplot(3,1,1:2)
hold on; grid on;
title('East Position Comparison')
plot(timeVec, posNED(outData,2), 'b');
plot(timeVec, filterEstX(outData,2), 'r');
plot(timeVec, filterEstX(outData,2) + sqrt(filterCovPosY(outData)), 'k');
plot(timeVec, filterEstX(outData,2) - sqrt(filterCovPosY(outData)), 'k');
xlabel('Time (s)');
ylabel('Y Pos East (m)');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(3) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosY(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovPosY)), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma (2,2)');
ax(4) = gca;

figure;
subplot(3,1,1:2)
hold on; grid on;
title('Down Position Comparison')
plot(timeVec, posNED(outData,3), 'b');
plot(timeVec, filterEstX(outData,3), 'r');
plot(timeVec, filterEstX(outData,3) + sqrt(filterCovPosZ(outData)), 'k');
plot(timeVec, filterEstX(outData,3) - sqrt(filterCovPosZ(outData)), 'k');
xlabel('Time (s)');
ylabel('Z Pos Down (m');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(5) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosZ(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovPosZ)), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma (3,3)');
ax(6) = gca;

figure;
subplot(3,1,1:2)
hold on; grid on;
title('North Velocity Comparison')
plot(timeVec, velNED(outData,1), 'b');
plot(timeVec, filterEstX(outData,4), 'r');
plot(timeVec, filterEstX(outData,4) + sqrt(filterCovVelX(outData)), 'k');
plot(timeVec, filterEstX(outData,4) - sqrt(filterCovVelX(outData)), 'k');
xlabel('Time (s)');
ylabel('X Vel North (m');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(7) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelX(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelX)), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma (4,4)');
ax(8) = gca;

figure;
subplot(3,1,1:2)
hold on; grid on;
title('East Velocity Comparison')
plot(timeVec, velNED(outData,2), 'b');
plot(timeVec, filterEstX(outData,5), 'r');
plot(timeVec, filterEstX(outData,5) + sqrt(filterCovVelY(outData)), 'k');
plot(timeVec, filterEstX(outData,5) - sqrt(filterCovVelY(outData)), 'k');
xlabel('Time (s)');
ylabel('Y Vel East (m');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(9) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelY(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelY)), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma (5,5)');
ax(10) = gca;

figure;
subplot(3,1,1:2)
hold on; grid on;
title('Down Velocity Comparison')
plot(timeVec, velNED(outData,3), 'b');
plot(timeVec, filterEstX(outData,6), 'r');
plot(timeVec, filterEstX(outData,6) + sqrt(filterCovVelZ(outData)), 'k');
plot(timeVec, filterEstX(outData,6) - sqrt(filterCovVelZ(outData)), 'k');
xlabel('Time (s)');
ylabel('Z Vel Down (m');
legend('Truth', 'KF Est', '1\sigma Uncert')
ax(11) = gca;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelZ(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelZ)), 2), 'k')
legend('1\sigma Value', 'Mean')
ylabel('\sigma (6,6)');
ax(12) = gca;

figure;
subplot(3,1,1);
hold on; grid on; grid minor;
plot(timeVec, posNED(outData,1)-filterEstX(outData,1), 'b');
title('Position Errors')
ylabel('North Pos Err (m)');
ax(13) = gca;
subplot(3,1,2);
hold on; grid on; grid minor;
plot(timeVec, posNED(outData,2)-filterEstX(outData,2), 'b');
ylabel('East Pos Err (m)');
ax(14) = gca;
subplot(3,1,3);
hold on; grid on; grid minor;
plot(timeVec, posNED(outData,3)-filterEstX(outData,3), 'b');
ylabel('Down Pos Err (m)');
xlabel('Time (s)');
ax(15) = gca;

figure;
subplot(3,1,1);
hold on; grid on; grid minor;
plot(timeVec, velNED(outData,1)-filterEstX(outData,4), 'b');
title('Velocity Errors')
ylabel('North Vel Err (m)');
ax(16) = gca;
subplot(3,1,2);
hold on; grid on; grid minor;
plot(timeVec, velNED(outData,2)-filterEstX(outData,5), 'b');
ylabel('East Vel Err (m)');
ax(17) = gca;
subplot(3,1,3);
hold on; grid on; grid minor;
plot(timeVec, velNED(outData,3)-filterEstX(outData,6), 'b');
ylabel('Down Vel Err (m)');
xlabel('Time (s)');
ax(18) = gca;


linkaxes(ax, 'x')

%% Save figures that have been generated
h =  findobj('type','figure');
n = length(h);
for i = 1:1:n
    saveas(figure(i), sprintf('save_for_later_figure%i', i), 'png' )
    saveas(figure(i), sprintf('save_for_later_figure%i', i), 'fig' )
end