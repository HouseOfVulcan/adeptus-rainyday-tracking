function [x_out, P_out, z_pred, S_inv, nis] = kf_cubature_cv(x_in, P_in, z_meas, dt, is_radar, q_accel, r_pos, q_scale)
% KF_CUBATURE_CV  Robust Cubature Kalman Filter for 3D Tracking
%
% This function implements a Cubature Kalman Filter (CKF) for a Constant Velocity (CV) model.
% It is superior to the Extended Kalman Filter (EKF) for radar tracking because it handles
% non-linearities (Polar <-> Cartesian) using statistical sigma points rather than
% potentially unstable Jacobians.
%
% =========================================================================
% INPUTS
% =========================================================================
% x_in      : [6x1] PRIOR STATE VECTOR [x; vx; y; vy; z; vz]
%             The "best guess" of the target's state from the previous time step.
%
% P_in      : [6x6] PRIOR COVARIANCE MATRIX
%             Represents the uncertainty of x_in. 
%             - Diagonal elements are variances (uncertainty^2).
%             - Off-diagonal elements are correlations between states.
%
% z_meas    : [3x1] MEASUREMENT VECTOR
%             The raw data received from the sensor for this specific detection.
%             - If is_radar=true : [Azimuth (deg); Elevation (deg); Range (m)]
%             - If is_radar=false: [x (m); y (m); z (m)]
%
% dt        : [Scalar] TIME STEP (seconds)
%             Time elapsed since the last update for this specific track.
%             (Note: In GNN, this can vary if a track misses a detection).
%
% is_radar  : [Boolean] SENSOR MODE FLAG
%             - true : Performs non-linear conversion (Spherical <-> Cartesian).
%             - false: Assumes linear measurement (Cartesian <-> Cartesian).
%
% q_accel   : [Scalar] PROCESS NOISE INTENSITY (m/s^2)
%             The "Tunable Parameter" of the filter.
%             - Describes the magnitude of random accelerations (maneuvers) expected.
%             - High Value: Filter is responsive but noisy (follows turns well).
%             - Low Value : Filter is smooth but lags during turns.
%
% r_pos     : [3x1] MEASUREMENT NOISE STANDARD DEVIATION
%             The sensor's accuracy. 
%             - Example: [1.0; 1.0; 10.0] means 1 deg accuracy in angle, 10m in range.
%
% q_scale   : [Scalar] (Optional, Default=1.0) ADAPTIVE NOISE SCALER
%             A multiplier to temporarily boost process noise.
%             - Use q_scale > 1.0 when the GNN detects a potential maneuver or 
%               when the NIS is high, allowing the filter to "catch up" quickly.
%
% =========================================================================
% OUTPUTS
% =========================================================================
% x_out     : [6x1] POSTERIOR STATE (Updated Estimate after seeing measurement)
% P_out     : [6x6] POSTERIOR COVARIANCE (Updated Uncertainty)
% z_pred    : [3x1] PREDICTED MEASUREMENT (What we expected the sensor to see)
% S_inv     : [3x3] INVERSE INNOVATION COVARIANCE (Needed for GNN Cost Matrix)
% nis       : [Scalar] NORMALIZED INNOVATION SQUARED
%             - This is the "Mahalanobis Distance" squared.
%             - Used for GATING: If nis > 9.0 (approx), reject the match.
%             - Used for SCORING: This is the cost value for the assignment matrix.
% =========================================================================

    % --- 0. INITIALIZATION & DEFAULTS ---
    if nargin < 8, q_scale = 1.0; end % Default to no scaling
    
    n = size(x_in, 1); % State Dimensions (6)
    m = size(z_meas, 1); % Measurement Dimensions (3)

    %% 1. DEFINE DYNAMIC MODEL (Constant Velocity)
    % The State Transition Matrix (F)
    % Projects position based on velocity: Pos_new = Pos_old + Vel * dt
    F = [1 dt 0  0  0  0; 
         0  1  0  0  0  0;
         0  0  1 dt  0  0;
         0  0  0  1  0  0;
         0  0  0  0  1 dt;
         0  0  0  0  0  1];

    % The Process Noise Matrix (Q)
    % Models the uncertainty added by target maneuvers over time dt.
    % We apply q_scale^2 here to allow adaptive maneuvering logic.
    current_q = (q_accel * q_scale)^2; 
    
    % Discrete time noise block (Piecewise White Noise Acceleration Model)
    q_block = [(dt^3)/3, (dt^2)/2; 
               (dt^2)/2,  dt     ] * current_q;
           
    Q = blkdiag(q_block, q_block, q_block); % Apply to X, Y, and Z axes

    % The Measurement Noise Matrix (R)
    % Sensor variance matrix (Standard Deviation squared)
    R = diag(r_pos.^2);

    %% 2. PREDICTION STEP (Time Update)
    % Project the state and uncertainty forward in time.
    x_pred = F * x_in;
    P_pred = F * P_in * F' + Q;
    
    % Stability Fix: Enforce symmetry and positive definiteness.
    % Floating point errors can make P asymmetric, which crashes Cholesky decomp.
    P_pred = (P_pred + P_pred') / 2 + eye(n) * 1e-9; 

    %% 3. GENERATE CUBATURE POINTS (Sigma Points)
    % We generate 2*n (12) statistical points around the mean to represent the 
    % probability distribution through the non-linear transformation.
    nPts = 2 * n;
    
    % Standard Cubature weights
    xi = [eye(n), -eye(n)] * sqrt(n);
    
    % Calculate Matrix Square Root of Covariance
    try
        % Cholesky is fast and standard
        S_chol = chol(P_pred, 'lower');
    catch
        % Robust Fallback: If P is slightly degenerate, use Eigendecomposition
        % This prevents the simulation from crashing if a track becomes perfect.
        [V, D] = eig(P_pred);
        S_chol = V * sqrt(max(D, 0));
    end
    
    % Create the cloud of points: Mean + (Spread * Geometry)
    pts = repmat(x_pred, 1, nPts) + S_chol * xi;

    %% 4. MEASUREMENT PROPAGATION
    % Convert every Cubature Point into "Measurement Space"
    % (e.g., Transform [x,y,z] points into [Az,El,Range] points)
    
    Z_pts = zeros(m, nPts);
    
    for i = 1:nPts
        px = pts(1,i); 
        py = pts(3,i); 
        pz = pts(5,i);
        
        if is_radar
            % --- NON-LINEAR RADAR MODEL ---
            
            % Protection: Avoid divide-by-zero if target is at range 0 (rare but fatal)
            rng_val = sqrt(px^2 + py^2 + pz^2);
            rng_val = max(rng_val, 1e-6); 
            
            % Calculate Azimuth and Elevation
            az_val = atan2d(py, px);
            el_val = asind(pz / rng_val);
            
            Z_pts(:, i) = [az_val; el_val; rng_val];
        else
            % --- LINEAR POSITION MODEL ---
            Z_pts(:, i) = [px; py; pz];
        end
    end
    
    % Predicted Measurement (The mean of the transformed points)
    z_pred = sum(Z_pts, 2) / nPts;

    %% 5. INNOVATION COVARIANCE & NIS CALCULATION
    % Calculate how the points spread around the predicted mean (Covariance).
    
    Pzz = R;            % Initialize with Sensor Noise
    Pxz = zeros(n, m);  % Initialize Cross-Covariance
    
    for i = 1:nPts
        % Residual: Point - Mean
        dz = Z_pts(:,i) - z_pred;
        
        % CRITICAL: Handle Angle Wrapping for Radar
        % (359 degrees is close to 1 degree, not far away)
        if is_radar
            dz(1) = mod(dz(1) + 180, 360) - 180;
            dz(2) = mod(dz(2) + 90, 180) - 90;
        end
        
        dx = pts(:,i) - x_pred;
        
        Pzz = Pzz + (dz * dz') / nPts;
        Pxz = Pxz + (dx * dz') / nPts;
    end
    
    % S_inv: Inverse Innovation Covariance
    % We compute this here so the GNN doesn't have to re-compute it repeatedly.
    % Using (A+A')/2 ensures symmetry for stable inversion.
    S_inv = inv((Pzz + Pzz') / 2);

    % --- CALCULATE INNOVATION (RESIDUAL) ---
    innovation = z_meas - z_pred;
    
    % Wrap Innovation angles if Radar
    if is_radar
        innovation(1) = mod(innovation(1) + 180, 360) - 180;
        innovation(2) = mod(innovation(2) + 90, 180) - 90;
    end

    % --- CALCULATE NIS (NORMALIZED INNOVATION SQUARED) ---
    % NIS = Innovation' * InverseCovariance * Innovation
    % This is the primary metric for GATING and GNN SCORING.
    nis = innovation' * S_inv * innovation;

    %% 6. CORRECTION STEP (Measurement Update)
    % Compute the Kalman Gain
    K = Pxz * S_inv;
    
    % Update State Estimate
    x_out = x_pred + K * innovation;
    
    % Update Covariance Estimate (Joseph Form subtraction)
    P_out = P_pred - K * Pzz * K';
    
    % Final Stability Check: Ensure P remains symmetric
    P_out = (P_out + P_out') / 2;
    
end