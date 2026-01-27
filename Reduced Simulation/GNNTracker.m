% GNNTRACKER.M
% =========================================================================
% GNN TRACKER (Global Nearest Neighbor)
% =========================================================================
% PURPOSE:
%   Acts as the central "intelligence" of the system. It manages the 
%   lifecycle of multiple target trajectories simultaneously, using Global 
%   Nearest Neighbor (GNN) association and Cubature Kalman Filtering (CKF) 
%   to estimate future target states.
%
% HOW IT CONNECTS:
%   - Input:  Filtered detections from DetectionGenerator.m.
%   - Math:   Calls CubatureKalmanFilter.m for prediction and correction steps.
%   - Output: Maintains a list of 'Confirmed' and 'Potential' track structures.
%
% KEY LOGIC:
%   - Track Lifecycle: Implements a "M-of-N" style logic where potential 
%     tracks must prove consistency (ConfirmThreshold) before being promoted.
%   - GNN Association: Uses a greedy matching algorithm (matchpairs) combined 
%     with Mahalanobis distance gating to assign detections to tracks.
%   - Merge Logic: Automatically consolidates redundant tracks that share 
%     the same physical space to prevent "track ghosting."
%   - Time-Invariance: All thresholds (deletion time, confirmation time) are 
%     automatically scaled by the simulation time step (dt).
% =========================================================================
% GNN TRACKER PARAMETER & TUNING GUIDE
% =========================================================================
% This section defines the behavioral constants that dictate how the 
% tracker interprets sensor data and manages target lifecycles.
%
% 1. ASSOCIATION & GATING
%   - Gate: The Chi-Square threshold for Mahalanobis distance. It defines 
%      the statistical "safety zone" around a prediction.
%   - MaxSpeed: The physical speed limit (m/s) used for coarse gating to 
%      quickly reject impossible detections.
%   - MergeThresholdSq: The squared distance limit used to identify and 
%      consolidate redundant tracks occupying the same space.
%
% 2. FILTRATION & NOISE (CKF)
%   - Q_Cruise: The process noise spectral density used during steady, 
%      linear flight.
%   - Q_Maneuver: Higher process noise used when a target is suspected 
%      of turning or accelerating.
%   - ManeuverThresh: The NIS (Normalized Innovation Squared) threshold 
%      that triggers the switch to high-gain maneuver tracking.
%
% 3. LIFECYCLE MANAGEMENT
%   - targetConfirmTime: The required duration of consistent detections 
%      (in seconds) to promote a Potential track to Confirmed.
%   - targetPotCoastTime: How long a Potential track can survive without 
%      a detection before being deleted.
%   - DeleteTime: The maximum "coasting" duration for a Confirmed track 
%      to persist during sensor dropouts or occlusions.
%
% 4. DYNAMIC SCALING
%   - ConfirmThreshold: The calculated number of consecutive hits needed 
%      based on the current dt (ConfirmTime / dt).
%   - PotentialMaxMisses: The frame-based limit for potential track 
%      deletion (PotentialCoastTime / dt).
% =========================================================================
classdef GNNTracker < handle
   
    properties
        Tracks, NextTrackID, PotentialTracks 
        Gate, DeleteTime, MaxSpeed        
        Q_Cruise, Q_Maneuver, ManeuverThresh    
        
        % Calculated Per-Frame Thresholds
        ConfirmThreshold 
        PotentialMaxMisses 
        MergeThresholdSq
    end
    
    methods
        function obj = GNNTracker(dt, tuningConfig)
            if nargin < 1, dt = 1.0; end
            
            obj.Tracks = struct('ID', {}, 'State', {}, 'Covariance', {}, ...
                                'Misses', {}, 'Status', {}, 'History', {});
            obj.PotentialTracks = struct('State', {}, 'Covariance', {}, ...
                                         'History', {}, 'Misses', {}, 'ConsecutiveHits', {}); 
            obj.NextTrackID = 1;

            % Defaults
            obj.Gate = 45.0; 
            obj.DeleteTime = 2.5; % Seconds
            obj.MaxSpeed = 1200.0; 
            obj.Q_Cruise = 5.0; 
            obj.Q_Maneuver = 4000.0; 
            obj.ManeuverThresh = 4.0;
            
            % Time-Based Defaults
            targetConfirmTime = 1.5; 
            targetPotCoastTime = 1.0;
            
            % Apply Config Overrides
            if nargin > 1 && ~isempty(tuningConfig)
                fnames = fieldnames(tuningConfig);
                for i = 1:length(fnames)
                    if isprop(obj, fnames{i})
                        obj.(fnames{i}) = tuningConfig.(fnames{i});
                    elseif strcmp(fnames{i}, 'ConfirmTime')
                        targetConfirmTime = tuningConfig.ConfirmTime;
                    elseif strcmp(fnames{i}, 'PotentialCoastTime')
                        targetPotCoastTime = tuningConfig.PotentialCoastTime;
                    end
                end
            end
            
            obj.MergeThresholdSq = (obj.Gate * 2)^2;

            % --- CONVERT SECONDS TO FRAMES ---
            % We use ceil() to ensure we wait at least the specified time.
            % We use max(3, ...) to ensure we physically have enough points for velocity.
            obj.ConfirmThreshold = max(3, ceil(targetConfirmTime / dt));
            
            % Potential tracks die if they miss for > X seconds
            obj.PotentialMaxMisses = ceil(targetPotCoastTime / dt);
        end
        
        function update(obj, detections, dt)
            
            % 1. Predict
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
            
            % 2. Association
            [matches, unassignedT, unassignedD] = obj.global_association(obj.Tracks, detections, obj.Gate, dt);
            
            % 3. Update Tracks
            for k = 1:size(matches, 1)
                tIdx = matches(k,1); dIdx = matches(k,2);
                z = detections{dIdx}.Measurement;
                R = detections{dIdx}.MeasurementNoise;
                t = obj.Tracks(tIdx);
                nis = CubatureKalmanFilter.calcNIS(t.State, t.Covariance, z, R);
                if nis < obj.ManeuverThresh, t.Covariance = t.Covariance * 0.5; end
                [t.State, t.Covariance] = CubatureKalmanFilter.correct(t.State, t.Covariance, z, R);
                t.Misses = 0;
                obj.Tracks(tIdx) = t;
            end
            for i = 1:length(unassignedT)
                obj.Tracks(unassignedT(i)).Misses = obj.Tracks(unassignedT(i)).Misses + 1;
            end

            % --- STEP 4: POTENTIAL TRACKS ---
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
                
                obj.PotentialTracks(pIdx).Misses = 0;
                obj.PotentialTracks(pIdx).ConsecutiveHits = obj.PotentialTracks(pIdx).ConsecutiveHits + 1;
                
                if obj.PotentialTracks(pIdx).ConsecutiveHits >= obj.ConfirmThreshold
                    obj.promote_track(pIdx);
                end
            end
            
            % Manage Misses for Potentials
            for i = 1:length(unassignedPot)
                idx = unassignedPot(i);
                obj.PotentialTracks(idx).ConsecutiveHits = 0; 
                obj.PotentialTracks(idx).Misses = obj.PotentialTracks(idx).Misses + 1;
            end
            
            % Track Deletion
            keepPot = [obj.PotentialTracks.Misses] <= obj.PotentialMaxMisses;
            obj.PotentialTracks = obj.PotentialTracks(keepPot);
            
            % Create New
            for i = 1:length(unAssDetsFinal)
                obj.create_potential(leftoverDets{unAssDetsFinal(i)});
            end
            
            % 5. Cleanup Confirmed Tracks
            maxMisses = ceil(obj.DeleteTime / dt);
            keep = [obj.Tracks.Misses] <= maxMisses;
            obj.Tracks = obj.Tracks(keep);
            obj.merge_duplicate_tracks();
        end
        
        function [matches, unT, unD] = global_association(obj, tracks, dets, gate, dt)
            nT = length(tracks); nD = length(dets);
            if nT == 0 || nD == 0
                matches=[]; unT=(1:nT)'; unD=(1:nD)'; return; 
            end
            costMat = inf(nT, nD);
            minGate = 4 * max(sqrt(diag(dets{1}.MeasurementNoise))); 
            searchRadius = max(obj.MaxSpeed * dt * 1.5, minGate);
            searchRadiusSq = searchRadius^2;
            hasCandidate = false;
            for t=1:nT
                tState = tracks(t).State; tPos = tState([1,3,5]); tCov = tracks(t).Covariance;
                for d=1:nD
                    dPos = dets{d}.Measurement;
                    dx=tPos(1)-dPos(1); dy=tPos(2)-dPos(2); dz=tPos(3)-dPos(3);
                    if (dx*dx + dy*dy + dz*dz) > searchRadiusSq, continue; end
                    nis = CubatureKalmanFilter.calcNIS(tState, tCov, dPos, dets{d}.MeasurementNoise);
                    if nis < gate
                        costMat(t,d) = nis; hasCandidate = true;
                    end
                end
            end
            if ~hasCandidate
                matches=[]; unT=(1:nT)'; unD=(1:nD)'; return;
            end
            [matches, unT, unD] = matchpairs(costMat, gate);
            valid = costMat(sub2ind(size(costMat), matches(:,1), matches(:,2))) < inf;
            matches = matches(valid,:);
            unT = setdiff(1:nT, matches(:,1))';
            unD = setdiff(1:nD, matches(:,2))';
        end
        
        function promote_track(obj, pIdx)
            pot = obj.PotentialTracks(pIdx);
            newTrk.ID = obj.NextTrackID;
            newTrk.State = pot.State; newTrk.Covariance = pot.Covariance;
            newTrk.Misses = 0; newTrk.Status = 'Confirmed'; newTrk.History = [];
            obj.Tracks(end+1) = newTrk;
            obj.NextTrackID = obj.NextTrackID + 1;
            obj.PotentialTracks(pIdx).Misses = 999; 
        end
        
        function create_potential(obj, det)
            if length(obj.PotentialTracks) > 50, return; end
            meas = det.Measurement;
            newPot.State = [meas(1); 0; meas(2); 0; meas(3); 0];
            newPot.Covariance = diag([100, 10000, 100, 10000, 100, 10000]);
            newPot.History = true; newPot.Misses = 0;
            newPot.ConsecutiveHits = 1; 
            obj.PotentialTracks(end+1) = newPot;
        end
        
        function merge_duplicate_tracks(obj)
             if length(obj.Tracks) < 2, return; end
             keepMask = true(length(obj.Tracks), 1);
             for i=1:length(obj.Tracks)
                 if ~keepMask(i), continue; end
                 p1 = obj.Tracks(i).State([1,3,5]);
                 for j=(i+1):length(obj.Tracks)
                     if ~keepMask(j), continue; end
                     p2 = obj.Tracks(j).State([1,3,5]);
                     dx=p1(1)-p2(1); dy=p1(2)-p2(2); dz=p1(3)-p2(3);
                     if (dx*dx + dy*dy + dz*dz) < obj.MergeThresholdSq 
                         keepMask(j) = false; 
                     end
                 end
             end
             obj.Tracks = obj.Tracks(keepMask);
        end
    end
end