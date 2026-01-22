classdef GNNTracker < handle
    % =========================================================================
    % ADAPTIVE-Q CUBATURE KALMAN FILTER (CKF) GNN TRACKER
    % =========================================================================
    %
    % SYSTEM OVERVIEW:
    %   This is a "Global Nearest Neighbor" (GNN) tracker that uses a 
    %   Constant Velocity (CV) Cubature Kalman Filter (CKF).
    %
    %   It employs a unique "Adaptive Process Noise" (Adaptive-Q) strategy
    %   to solve the "Tracker's Dilemma" (Stability vs. Agility) without
    %   requiring the complexity of an IMM (Interacting Multiple Model) filter.
    %
    % KEY ALGORITHMS:
    %   1. PREDICTION: 
    %      Predicts all tracks forward using High Process Noise (Q). This keeps
    %      the validation gate wide, ensuring we don't lose maneuvering targets
    %      (like Spirals) during the Association step.
    %
    %   2. GNN ASSOCIATION:
    %      Uses the Hungarian Algorithm (Munkres) to globally minimize the 
    %      Normalized Innovation Squared (NIS) cost between all tracks and detections.
    %      Includes a "Coarse Gate" (Euclidean) for speed and a "Fine Gate" (NIS)
    %      for statistical validity.
    %
    %   3. RETROACTIVE Q-CORRECTION:
    %      After association, it checks the NIS of the match.
    %      - If NIS is Low (Target is straight): It artificially tightens the
    %        covariance, simulating as if we had predicted with Low Q.
    %        (Result: Smooth, stable tracks).
    %      - If NIS is High (Target is turning): It keeps the High Q prediction.
    %        (Result: Responsive, non-breaking tracks during turns).
    %
    %   4. LIFECYCLE MANAGEMENT:
    %      - Promotion: M-of-N logic (3 hits in N scans) to confirm a track.
    %      - Deletion:  Coasting time logic (delete after X seconds of misses).
    %
    % =========================================================================
    
    properties
        Tracks          % Struct Array of Active Tracks
        NextTrackID     % ID Counter
        PotentialTracks % Struct Array of Tentative (Unconfirmed) Tracks
        
        % --- TUNING PARAMETERS ---
        Gate            % [NIS] Statistical Gate Threshold. 
                        % Determines "How many sigmas" a measurement can be from 
                        % the track before being rejected.
                        
        DeleteTime      % [Seconds] Coasting Duration.
                        % How long to keep a track alive without updates.
                        
        MaxSpeed        % [m/s] Maximum physically possible speed.
                        % Used for "Coarse Gating" optimization.
        
        % --- ADAPTIVE PROCESS NOISE SETTINGS ---
        Q_Cruise        % [m^2/s^3] Spectral Density for Straight Flight.
                        % Low value = High trust in model, smooth track.
                        
        Q_Maneuver      % [m^2/s^3] Spectral Density for Maneuvers.
                        % High value = Low trust in model, accepts turns.
                        
        ManeuverThresh  % [NIS] Threshold to trigger Adaptive Logic.
                        % If error < this, force smoothing (Cruise mode).
    end
    
    methods
        function obj = GNNTracker(dt)
            if nargin < 1, dt = 1.0; end
            
            obj.Tracks = struct('ID', {}, 'State', {}, 'Covariance', {}, ...
                                'Misses', {}, 'Status', {}, 'History', {});
            obj.PotentialTracks = struct('State', {}, 'Covariance', {}, ...
                                         'History', {}, 'Misses', {}); 
            obj.NextTrackID = 1;
            
            % --- CONFIGURATION ---
            % Gate = 45.0 (NIS). Sqrt(45) approx 6.7 Sigma.
            % We use a wide gate because we rely on the Adaptive Q to tighten
            % the track AFTER association, rather than restricting entry.
            obj.Gate = 60.0;           
            
            % Deletion: 2.5s. Allows a track to coast through a brief
            % period of dropped detections (e.g., occlusion).
            obj.DeleteTime      = 2.5; 
            
            % MaxSpeed: 1200 m/s (Mach 3.5). Safe upper limit for aircraft.
            obj.MaxSpeed        = 1500.0; 
            
            % --- ADAPTIVE Q TUNING ---
            obj.Q_Cruise       = 5.0;    % Very smooth (Rail-straight)
            obj.Q_Maneuver     = 6000.0; % Very agile (Catch 9G turns)
            obj.ManeuverThresh = 4.0;    % Approx 2-Sigma. 
                                         % If error is within 2-Sigma, assume straight.
        end
        
        function update(obj, detections, dt)
            % -----------------------------------------------------------
            % STEP 1: PREDICT (High Agility Mode)
            % -----------------------------------------------------------
            % Predict all tracks using Q_Maneuver. This inflates the covariance
            % P, effectively widening the gate. This ensures that if a target
            % turns sharply, it is still likely to fall within the gate.
            
            for i = 1:length(obj.Tracks)
                [obj.Tracks(i).State, obj.Tracks(i).Covariance] = ...
                    CubatureKalmanFilter.predict(obj.Tracks(i).State, ...
                                                 obj.Tracks(i).Covariance, ...
                                                 dt, obj.Q_Maneuver);
            end
            for i = 1:length(obj.PotentialTracks)
                [obj.PotentialTracks(i).State, obj.PotentialTracks(i).Covariance] = ...
                    CubatureKalmanFilter.predict(obj.PotentialTracks(i).State, ...
                                                 obj.PotentialTracks(i).Covariance, ...
                                                 dt, obj.Q_Maneuver);
            end
            
            % -----------------------------------------------------------
            % STEP 2: GLOBAL ASSOCIATION
            % -----------------------------------------------------------
            % Match predictions to measurements using Hungarian Algorithm.
            [matches, unassignedT, unassignedD] = obj.global_association(obj.Tracks, detections, obj.Gate, dt);
            
            % -----------------------------------------------------------
            % STEP 3: ADAPTIVE UPDATE (The "Magic" Step)
            % -----------------------------------------------------------
            for k = 1:size(matches, 1)
                tIdx = matches(k,1); dIdx = matches(k,2);
                z = detections{dIdx}.Measurement;
                R = detections{dIdx}.MeasurementNoise;
                t = obj.Tracks(tIdx);
                
                % Check: How good was this match really?
                nis = CubatureKalmanFilter.calcNIS(t.State, t.Covariance, z, R);
                
                % Retroactive Correction:
                % If the match was tight (Low NIS), our High-Q prediction was
                % pessimistic. We scale down the covariance to simulate 
                % what would have happened if we used Q_Cruise.
                if nis < obj.ManeuverThresh
                    t.Covariance = t.Covariance * 0.5; 
                end
                
                % Standard Kalman Correction
                [t.State, t.Covariance] = CubatureKalmanFilter.correct(t.State, t.Covariance, z, R);
                t.Misses = 0;
                obj.Tracks(tIdx) = t;
            end
            
            % Increment miss counters for unmatched tracks
            for i = 1:length(unassignedT)
                obj.Tracks(unassignedT(i)).Misses = obj.Tracks(unassignedT(i)).Misses + 1;
            end
            
            % -----------------------------------------------------------
            % STEP 4: TRACK INITIATION (M/N Logic)
            % -----------------------------------------------------------
            % Unassigned detections attempt to match with Potential Tracks.
            leftoverDets = detections(unassignedD);
            [potMatches, unassignedPot, unAssDetsFinal] = ...
                obj.global_association(obj.PotentialTracks, leftoverDets, obj.Gate, dt);
            
            for k = 1:size(potMatches, 1)
                pIdx = potMatches(k,1); dIdx = potMatches(k,2);
                [obj.PotentialTracks(pIdx).State, obj.PotentialTracks(pIdx).Covariance] = ...
                    CubatureKalmanFilter.correct(obj.PotentialTracks(pIdx).State, ...
                                                 obj.PotentialTracks(pIdx).Covariance, ...
                                                 leftoverDets{dIdx}.Measurement, ...
                                                 leftoverDets{dIdx}.MeasurementNoise);
                
                % Record "Hit" in history
                obj.PotentialTracks(pIdx).History(end+1) = true;
                obj.PotentialTracks(pIdx).Misses = 0;
                
                % CONFIRMATION: If we have >= 3 hits, promote to Track
                if sum(obj.PotentialTracks(pIdx).History) >= 3
                    obj.promote_track(pIdx);
                end
            end
            
            % Prune missed potential tracks immediately (hard logic)
            obj.PotentialTracks(unassignedPot) = [];
            
            % Create new potentials from remaining raw detections
            for i = 1:length(unAssDetsFinal)
                obj.create_potential(leftoverDets{unAssDetsFinal(i)});
            end
            
            % -----------------------------------------------------------
            % STEP 5: CLEANUP
            % -----------------------------------------------------------
            % Remove confirmed tracks that have coasted too long.
            maxMisses = ceil(obj.DeleteTime / dt);
            keep = [obj.Tracks.Misses] <= maxMisses;
            obj.Tracks = obj.Tracks(keep);
            
            % Deduplicate tracks that have merged onto the same target
            obj.merge_duplicate_tracks();
        end
        
        function [matches, unT, unD] = global_association(obj, tracks, dets, gate, dt)
            % GLOBAL_ASSOCIATION: Hungarian Algorithm Wrapper
            nT = length(tracks); nD = length(dets);
            costMat = inf(nT, nD);
            
            % Coarse Gating Radius (Euclidean optimization)
            searchRadius = obj.MaxSpeed * dt * 1.5; 
            
            for t=1:nT
                for d=1:nD
                    % 1. Coarse Gate (Fast Distance Check)
                    if norm(tracks(t).State([1,3,5]) - dets{d}.Measurement) > searchRadius
                        continue; 
                    end
                    % 2. Fine Gate (Statistical NIS Check)
                    nis = CubatureKalmanFilter.calcNIS(tracks(t).State, tracks(t).Covariance, ...
                        dets{d}.Measurement, dets{d}.MeasurementNoise);
                    
                    if nis < gate, costMat(t,d) = nis; end
                end
            end
            
            % Solve Assignment Problem
            [matches, unT, unD] = matchpairs(costMat, gate);
            
            % Filter out any matches that were technically 'infinite' cost
            valid = costMat(sub2ind(size(costMat), matches(:,1), matches(:,2))) < inf;
            matches = matches(valid,:);
            unT = setdiff(1:nT, matches(:,1))';
            unD = setdiff(1:nD, matches(:,2))';
        end
        
        function promote_track(obj, pIdx)
            % Moves a track from 'Potential' list to 'Confirmed' list
            pot = obj.PotentialTracks(pIdx);
            newTrk.ID = obj.NextTrackID;
            newTrk.State = pot.State; newTrk.Covariance = pot.Covariance;
            newTrk.Misses = 0; newTrk.Status = 'Confirmed'; newTrk.History = [];
            obj.Tracks(end+1) = newTrk;
            obj.NextTrackID = obj.NextTrackID + 1;
            obj.PotentialTracks(pIdx).Misses = 999; % Mark for deletion
        end
        
        function create_potential(obj, det)
            % Initializes a new potential track from a raw detection
            meas = det.Measurement;
            % Initialize State: [x; 0; y; 0; z; 0] (Zero velocity guess)
            newPot.State = [meas(1); 0; meas(2); 0; meas(3); 0];
            % Initialize Covariance: High uncertainty in Velocity (10000)
            newPot.Covariance = diag([100, 10000, 100, 10000, 100, 10000]);
            newPot.History = [true]; newPot.Misses = 0;
            obj.PotentialTracks(end+1) = newPot;
        end
        
        function merge_duplicate_tracks(obj)
            % Merges tracks that have converged to the same location.
            % (Simple Euclidean check < 50m)
             if length(obj.Tracks) < 2, return; end
             keepMask = true(length(obj.Tracks), 1);
             for i=1:length(obj.Tracks)
                 if ~keepMask(i), continue; end
                 for j=(i+1):length(obj.Tracks)
                     if ~keepMask(j), continue; end
                     if norm(obj.Tracks(i).State([1,3,5]) - obj.Tracks(j).State([1,3,5])) < 50
                         keepMask(j) = false; 
                     end
                 end
             end
             obj.Tracks = obj.Tracks(keepMask);
        end
    end
end