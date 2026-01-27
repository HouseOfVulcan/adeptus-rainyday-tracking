% CUBATUREKALMANFILTER.m
% =========================================================================
% CUBATURE KALMAN FILTER (CUSTOM MATH LIBRARY)
% =========================================================================
% PURPOSE:
%   A stateless, static math library providing the "Predict" and "Correct"
%   steps for the tracker. It implements the Spherical-Radial Cubature 
%   Kalman Filter (CKF), which is a non-linear estimator superior to the 
%   Extended Kalman Filter (EKF) for complex tracking problems.
%
% HOW IT CONNECTS:
%   - Called by: GNNTracker.m (exclusively).
%   - Input:  Current State vector (6x1) and Covariance matrix (6x6).
%   - Output: Updated State and Covariance estimates.
%
% KEY LOGIC:
%   - Stateless Design: The functions are pure math inputs/outputs. It does 
%     not store history, making it robust and easy to restart.
%   - Numerical Safety: Includes fallback logic (using eigendecomposition) 
%     if the covariance matrix becomes asymmetric or non-positive definite, 
%     preventing mathematical crashes during long simulations.
%   - NIS Calculation: Computes the Normalized Innovation Squared (Mahalanobis 
%     Distance), which is the critical metric used by the tracker to deciding 
%     which radar dots belong to which track.
% =========================================================================
% CUBATURE KALMAN FILTER PARAMETER & TUNING GUIDE
% =========================================================================
% This library is mathematically deterministic, but its performance depends 
% entirely on the quality of the two noise matrices passed into it:
%
% 1. PROCESS NOISE (q_intensity)
%    - Definition: Represents the "trust" in the Constant Velocity model.
%    - Tuning:
%       * Low Value (e.g., 5.0): Tells the filter "This object flies in a 
%         straight line." Result: Smooth tracks, but laggy on turns.
%       * High Value (e.g., 4000.0): Tells the filter "This object maneuvers 
%         violently." Result: Fast reaction to turns, but jittery path.
%
% 2. MEASUREMENT NOISE (R)
%    - Definition: Represents the "trust" in the sensor data.
%    - Tuning:
%       * Small Diagonals (e.g., 50^2): "The sensor is precise." The filter 
%         will jump to every detection point.
%       * Large Diagonals (e.g., 150^2): "The sensor is noisy." The filter 
%         will ignore individual jitters and average the path over time.
%
% 3. STABILITY CONTROLS
%    - Symmetric Enforcement: P = (P + P') / 2. This line prevents the 
%      covariance matrix from becoming asymmetric due to floating-point 
%      rounding errors, which would crash the Cholesky decomposition.
% =========================================================================
classdef CubatureKalmanFilter
    
    methods(Static)
        
        function [x_new, P_new] = predict(x, P, dt, q_intensity)
            % PREDICT: Projects state forward in time.
            %
            % MODEL: Constant Velocity (CV)
            %   X_k+1 = F * X_k + Noise
            %
            % INPUTS:
            %   x           : [6x1] State Vector [x; vx; y; vy; z; vz]
            %   P           : [6x6] Covariance Matrix
            %   dt          : [Scalar] Time step (seconds)
            %   q_intensity : [Scalar] Process Noise Spectral Density
            %
            % MATH DETAILS:
            %   We use the "Discrete White Noise Acceleration" model.
            %   Q is calculated dynamically based on dt^3/3 (pos) and dt (vel).
            
            % 1. Transition Matrix F (Newton's Laws)
            F = [1 dt 0  0 0  0;
                 0  1 0  0 0  0;
                 0  0 1 dt 0  0;
                 0  0 0  1 0  0;
                 0  0 0  0 1 dt;
                 0  0 0  0 0  1];
             
            % 2. Process Noise Q(dt)
            Q_block = [(dt^3)/3, (dt^2)/2; 
                       (dt^2)/2,  dt     ] * q_intensity;
            Q = blkdiag(Q_block, Q_block, Q_block);
            
            % 3. Standard Linear Prediction
            x_new = F * x;
            P_new = F * P * F' + Q;
            
            % 4. Force Symmetry (Stability)
            P_new = (P_new + P_new') / 2;
        end

        function [x_post, P_post] = correct(x_pred, P_pred, z_meas, R)
            % CORRECT: Updates state based on measurement.
            %
            % ALGORITHM: Unscented/Cubature Transform
            %   Uses Sigma Points to propagate the Gaussian through the 
            %   measurement model.
            
            % 1. Generate Sigma Points & Statistics
            [z_pred, S_inv, Pxz, Pzz] = CubatureKalmanFilter.getSigmaStats(x_pred, P_pred, R);
            
            % 2. Kalman Gain
            K = Pxz * S_inv;
            innovation = z_meas - z_pred;
            
            % 3. Update
            x_post = x_pred + K * innovation;
            P_post = P_pred - K * Pzz * K';
            
            % 4. Force Symmetry
            P_post = (P_post + P_post') / 2;
        end
        
        function nis = calcNIS(trackState, trackCov, z_meas, R)
            % CALCNIS: Calculates Normalized Innovation Squared (Mahalanobis)
            %
            % FORMULA: NIS = Innovation' * Inv(InnovationCov) * Innovation
            %
            % USAGE:
            %   Used for Gating (Chi-Square test) and GNN Cost Matrices.
            %   A dimensionless number representing "Sigma^2".
            
            [z_pred, S_inv, ~, ~] = CubatureKalmanFilter.getSigmaStats(trackState, trackCov, R);
            innovation = z_meas - z_pred;
            nis = innovation' * S_inv * innovation;
        end

        function [z_mean, S_inv, Pxz, Pzz] = getSigmaStats(x, P, R)
            % GETSIGMASTATS: Core Cubature Transform Logic
            % Generates 2n Sigma Points and propagates them.
            
            n = 6; m = 3; nPts = 2 * n;
            
            % 1. Matrix Square Root (Robust)
            try
                S_chol = chol(P, 'lower');
            catch
                % Fallback: Eigendecomposition if P is not perfectly positive definite
                [V, D] = eig(P);
                S_chol = V * sqrt(max(D, 0));
            end
            
            % 2. Generate Points: Mean +/- Sqrt(P)
            xi = [eye(n), -eye(n)] * sqrt(n);
            pts = repmat(x, 1, nPts) + S_chol * xi;
            
            % 3. Measurement Model (Linear Extraction)
            Z_pts = zeros(m, nPts);
            for i = 1:nPts
                Z_pts(:, i) = [pts(1, i); pts(3, i); pts(5, i)];
            end
            
            % 4. Reconstruct Mean and Covariance
            z_mean = sum(Z_pts, 2) / nPts;
            Pzz = R; Pxz = zeros(n, m);
            for i = 1:nPts
                dz = Z_pts(:, i) - z_mean;
                dx = pts(:, i) - x;
                Pzz = Pzz + (dz * dz') / nPts;
                Pxz = Pxz + (dx * dz') / nPts;
            end
            
            % 5. Inverse Innovation Covariance
            S_inv = inv((Pzz + Pzz') / 2);
        end
    end
end