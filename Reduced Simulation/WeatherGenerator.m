% WEATHERGENERATOR.m
% =========================================================================
% WEATHER GENERATION (Configuration Factory)
% =========================================================================
% PURPOSE:
%   The "Central Database" of simulation settings. It acts as a static 
%   factory that returns a configuration structure based on the requested 
%   difficulty level ('Clear', 'Rainy', 'Storm', etc.).
%
% HOW IT CONNECTS:
%   - Called by: MainReducedSim.m (at startup).
%   - Outputs: 'config' struct containing:
%       1. Sensor Params (Pd, R, Clutter) -> Sent to DetectionGenerator.
%       2. Physics Params (Turbulence) -> Sent to TruthGenerator.
%       3. Tracker Tuning (Gate, Q) -> Sent to GNNTracker.
%
% UNIQUE FEATURE:
%   - Adaptive Tuning: This file doesn't just make the weather harder; 
%     it also automatically re-tunes the GNNTracker (tightening gates, 
%     increasing process noise) to give the system a fighting chance 
%     against the harsher conditions.
% =========================================================================
% WEATHER & SENSOR TUNING GUIDE
% =========================================================================
% 1. ENVIRONMENTAL PARAMETERS (The "Problem")
%   - Pd (0.0 - 1.0): Probability of Detection. Lower values simulate 
%      sensor "blinking", where targets disappear for frames at a time.
%   - R (Covariance): Measurement noise. High values (e.g., 150^2) make 
%      dots jitter wildly around the true path.
%   - FalseAlarmRate: Clutter density. 1e-10 is roughly the "Survivable Limit"
%      before the tracker gets overwhelmed by ghosts.
%   - TurbulenceSigma: The strength of wind gusts pushing the planes.
%
% 2. TRACKER ADAPTATION (The "Solution")
%   - Gate: We tighten the gate (45 -> 25) in Clutter to prevent attaching 
%      false alarms, but this increases the risk of losing fast targets.
%   - Q_Maneuver: We increase process noise in Storms to tell the Kalman 
%      filter "expect heavy turbulence, don't trust the straight line."
%   - ConfirmTime: We wait longer (1.5s -> 2.0s) in rain to ensure a 
%      track is real before displaying it to the user.
%   - ConfirmDensity: (M-of-N) The percentage of hits required within the
%      window to confirm. Lower this in low-Pd environments.
% =========================================================================
classdef WeatherGenerator
    
    methods(Static)
        function config = getParams(scenarioType)
            % --- BASELINES ---
            config.Pd = 0.99; 
            config.R = diag([50, 50, 50].^2); 
            config.FalseAlarmRate = 1e-12;
            config.TurbulenceSigma = 0; 
            config.AttenuationFactor = 0;
            
            % --- DEFAULT TUNING ---
            tuning.DeleteTime     = 2.5;
            tuning.Gate           = 45;     
            tuning.Q_Cruise       = 20;      
            tuning.Q_Maneuver     = 4500;   
            tuning.ManeuverThresh = 4;      
            tuning.MaxSpeed       = 1200.0;
            tuning.ConfirmTime        = 1.99; % N = confirmtime / dt
            tuning.PotentialCoastTime = 1.0; 
            tuning.MaxPotentialTracks = 100;
            tuning.ConfirmDensity     = 0.79; % M = ceil(density * N)
            
            switch lower(scenarioType)
                case 'clear'
                    % Leave Defaults

                case 'rainy'
                    config.Pd = 0.98; 
                    config.FalseAlarmRate = 2e-10; 
                    config.R = diag([70, 70, 70].^2); 
                    config.AttenuationFactor = 2e-5; 
                    config.TurbulenceSigma = 40.0;    
                    
                    tuning.DeleteTime     = 2.5;
                    tuning.Gate           = 40;     
                    tuning.Q_Cruise       = 20;      
                    tuning.Q_Maneuver     = 9000;   
                    tuning.ManeuverThresh = 2.6;      
                    tuning.MaxSpeed       = 1200.0;
                    tuning.ConfirmTime        = 1.99; % N = confirmtime / dt
                    tuning.PotentialCoastTime = 0.9; 
                    tuning.MaxPotentialTracks = 100;
                    tuning.ConfirmDensity     = 0.84; % M = ceil(density * N)   
                    
                case 'fog'
                    config.Pd = 0.95; 
                    config.R = diag([60, 60, 60].^2);
                    config.FalseAlarmRate = 5e-11;
                    config.AttenuationFactor = 1e-5; 
                    config.TurbulenceSigma = 0;      

                    tuning.DeleteTime     = 2.5;
                    tuning.Gate           = 40;     
                    tuning.Q_Cruise       = 20;      
                    tuning.Q_Maneuver     = 4500;   
                    tuning.ManeuverThresh = 4;      
                    tuning.MaxSpeed       = 1200.0;
                    tuning.ConfirmTime        = 1.99; % N = confirmtime / dt
                    tuning.PotentialCoastTime = 1.0; 
                    tuning.MaxPotentialTracks = 100;
                    tuning.ConfirmDensity     = 0.79; % M = ceil(density * N)
                    
                case 'storm'
                    config.Pd = 0.90; 
                    config.R = diag([100, 100, 100].^2); 
                    config.FalseAlarmRate = 1e-10;        
                    config.AttenuationFactor = 5e-5; 
                    config.TurbulenceSigma = 80.0;   
                    
                    tuning.DeleteTime     = 2.5;
                    tuning.Gate           = 25;     
                    tuning.Q_Cruise       = 5000;      
                    tuning.Q_Maneuver     = 8000;   
                    tuning.ManeuverThresh = 6;      
                    tuning.MaxSpeed       = 1200.0;
                    tuning.ConfirmTime        = 3; % N = confirmtime / dt
                    tuning.PotentialCoastTime = 3; 
                    tuning.MaxPotentialTracks = 100;
                    tuning.ConfirmDensity     = 0.4; % M = ceil(density * N)
                    
                otherwise
                    warning('Unknown weather type. Using Clear defaults.');
            end
            
            config.TrackerTuning = tuning;
           
        end
    end
end