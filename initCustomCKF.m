function filter = initCustomCKF(detection)
% INITCUSTOMCKF Initializes a trackingCKF for Cartesian Input
%
% ADAPTATION:
%   - Handles Radar data that arrives as [x; y; z] (Scenario Frame)
%   - Uses a Linear Measurement Model inside the Cubature framework.
%   - Maintains robust Variable-DT motion logic.

    % 1. Parse Measurement (Cartesian Input)
    meas = detection.Measurement; % [x; y; z]
    
    % Initial State: [x; vx; y; vy; z; vz]
    % We map measurement directly to state positions.
    % We assume 0 velocity initially.
    initState = [meas(1); 0; meas(2); 0; meas(3); 0];
    
    % 2. Tuned Initial Covariance (P)
    % Position Variance: 10000 (100m std dev) - Trust measurement reasonably well.
    % Velocity Variance: 40000 (200m/s std dev) - High uncertainty to catch fast planes.
    posVar = 10000;
    velVar = 40000;
    initCov = diag([posVar, velVar, posVar, velVar, posVar, velVar]); 
    
    % 3. Measurement Noise (R)
    % If the scenario sends Cartesian data, detection.MeasurementNoise 
    % is already the correct 3x3 Cartesian covariance matrix.
    if isempty(detection.MeasurementNoise)
        R = diag([50 50 50].^2); % Default fallback (50m error)
    else
        R = detection.MeasurementNoise;
    end

    % 4. Create Filter
    % We use @cartMeasFcn instead of @radarMeasFcn
    filter = trackingCKF(@cvStateFcn, @cartMeasFcn, initState, ...
        'StateCovariance', initCov, ...
        'ProcessNoise', eye(6), ...            
        'MeasurementNoise', R, ...
        'HasAdditiveProcessNoise', false, ...  
        'HasAdditiveMeasurementNoise', true); 
end

% ==========================================
% CUSTOM MOTION MODEL (Robust Variable DT)
% ==========================================
function xNext = cvStateFcn(x, w, varargin)
    % Signature: f(x, w, dt)
    % Accepts varargin to absorb 'dt' safely.
    
    if ~isempty(varargin)
        dt = varargin{1};
    else
        dt = 0.1;
    end
    
    % Process Noise Tuning
    q_intensity = 40; % m/s^2
    
    % Transition Matrix (Constant Velocity)
    F = [1 dt 0  0 0  0;
         0  1 0  0 0  0;
         0  0 1 dt 0  0;
         0  0 0  1 0  0;
         0  0 0  0 1 dt;
         0  0 0  0 0  1];
     
    % Dynamic Noise Matrix Q(dt)
    Q_block = [(dt^3)/3, (dt^2)/2; 
               (dt^2)/2,  dt     ] * q_intensity;
           
    try
        L_block = chol(Q_block, 'lower');
    catch
        L_block = zeros(2);
    end
    
    L = blkdiag(L_block, L_block, L_block);
    xNext = F * x + L * w;
end

% ==========================================
% CUSTOM MEASUREMENT MODEL (Cartesian)
% ==========================================
function z = cartMeasFcn(x, varargin)
    % Signature: h(x, measurementParams)
    % Maps State [x; vx; y; vy; z; vz] -> Measurement [x; y; z]
    
    % Simple Linear Extraction
    z = [x(1); x(3); x(5)];
end