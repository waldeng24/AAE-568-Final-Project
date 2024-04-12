% Conversion constants
mph2fps = 5280/60/60;
mph2mps = 0.44704;

% Specifications 
spec.speedMilesPerHour = 15;
spec.speedFeetPerSec = spec.speedMilesPerHour * mph2fps;
spec.speedMetersPerSec = spec.speedMilesPerHour * mph2mps;
spec.thetaDotDegPerSec = 3;

% Create waypoints will be flown
waypointsLLA = [ ...
    32.2283506765183,  -110.94445142758335, 100; ...
    32.22806024925261, -110.95655355453485, 100; ...
    32.23579255859184, -110.95648918151916, 100; ...
    32.23568365739473, -110.94445142758335, 100; ...
    32.2283506765183,  -110.94445142758335, 100; ...
    ];

[X, Y, Z] = geodetic2ecef(wgs84Ellipsoid, waypointsLLA(:,1), waypointsLLA(:,2), waypointsLLA(:,3));
waypointsECEF = [X, Y, Z];

[xNorth,yEast,zDown] = geodetic2ned(waypointsLLA(:,1), waypointsLLA(:,2), waypointsLLA(:,3), ...
    waypointsLLA(1,1), waypointsLLA(1,2), waypointsLLA(1,3), wgs84Ellipsoid, 'degrees');
waypointsNED = [xNorth yEast zDown];

numWaypoints = height(waypointsLLA);

traj = waypointTrajectory( ...
    waypointsLLA, ...
    GroundSpeed=repelem(spec.speedMetersPerSec, numWaypoints));

initialPositionLLA = waypointsLLA(1,:);