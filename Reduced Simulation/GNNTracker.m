% GNNTRACKER.M
% =========================================================================
% GNN TRACKER (Global Nearest Neighbor) - ADAPTIVE & VECTORIZED
% =========================================================================
% PURPOSE:
%   State-of-the-art single-hypothesis tracker. It manages target lifecycles 
%   and uses "NIS-Driven Adaptive Process Noise" to balance high-precision 
%   smoothing during cruise with instant reaction to high-G maneuvers.
%
% CHANGES FROM PREVIOUS:
%   - FIXED: Added 'Age' field to PotentialTracks to prevent immediate 
%     deletion of new tracks (the "Infant Mortality" bug).
%   - LOGIC: Tracks are now immune to density pruning for the first few frames.
% =========================================================================
classdef GNNTracker < handle
   
    properties
        Tracks, NextTrackID, PotentialTracks 
        Gate, DeleteTime, MaxSpeed        
        Q_Cruise, Q_Maneuver, ManeuverThresh    
        
        ConfirmTime
        PotentialCoastTime
        ConfirmDensity

        % Logic Thresholds
        WindowSize, MinHits, PotentialMaxMisses 
        MergeThresholdSq, MaxPotentialTracks
        
        % Pre-calculated Physics Matrices (Cache)
        F                   % State Transition Matrix
        Q_CV_Mat            % Process Noise (Cruise)
        Q_Man_Mat           % Process Noise (Maneuver)
    end
    
    methods
        function obj = GNNTracker(dt, tuningConfig)
            if nargin < 1, dt = 0.1; end
            
            % --- INITIALIZATION ---
            obj.Tracks = struct('ID', {}, 'State', {}, 'Covariance', {}, ...
                                'Misses', {}, 'Status', {}, 'HitBits', {}, ...
                                'LastNIS', {});
            
            % Added 'Age' to track structure
            obj.PotentialTracks = struct('State', {}, 'Covariance', {}, ...
                                         'Misses', {}, 'HitBits', {}, ...
                                         'Age', {}); 
            obj.NextTrackID = 1;

            % --- DEFAULTS ---
            obj.Gate = 30.0;           
            obj.DeleteTime = 2.5;      
            obj.MaxSpeed = 1200.0;     
            obj.Q_Cruise = 20.0;       
            obj.Q_Maneuver = 4500.0;   
            obj.ManeuverThresh = 6.0;  
            
            obj.ConfirmTime = 2.0; 
            obj.PotentialCoastTime = 1.0;
            obj.ConfirmDensity = 0.8;
            
            % Apply Config Overrides
            if nargin > 1 && ~isempty(tuningConfig)
                fnames = fieldnames(tuningConfig);
                for i = 1:length(fnames)
                    if isprop(obj, fnames{i})
                        obj.(fnames{i}) = tuningConfig.(fnames{i});
                    end
                end
            end
            
            obj.MaxPotentialTracks = 100;
            
            % --- PRE-CALCULATE LOGIC ---
            obj.WindowSize = ceil(obj.ConfirmTime / dt); 
            obj.WindowSize = max(5, min(64, obj.WindowSize)); 
            obj.MinHits = ceil(obj.WindowSize * obj.ConfirmDensity); 
            obj.PotentialMaxMisses = ceil(obj.PotentialCoastTime / dt);
            obj.MergeThresholdSq = (100.0)^2; 

            % --- PRE-CALCULATE PHYSICS ---
            % State Transition (Newton)
            obj.F = [1 dt 0  0 0  0;
                     0  1 0  0 0  0;
                     0  0 1 dt 0  0;
                     0  0 0  1 0  0;
                     0  0 0  0 1 dt;
                     0  0 0  0 0  1];
                 
            % Pre-build both Q matrices
            q_block_cruise = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * obj.Q_Cruise;
            obj.Q_CV_Mat = blkdiag(q_block_cruise, q_block_cruise, q_block_cruise);
            
            q_block_man = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * obj.Q_Maneuver;
            obj.Q_Man_Mat = blkdiag(q_block_man, q_block_man, q_block_man);
        end
        
        function update(obj, detections, dt)
            
            % --- 1. ADAPTIVE PREDICTION STEP ---
            
            % A. Confirmed Tracks
            nC = length(obj.Tracks);
            if nC > 0
                allStates = [obj.Tracks.State];
                allStates = obj.F * allStates; % Vectorized
                
                for i = 1:nC
                    obj.Tracks(i).State = allStates(:, i);
                    
                    % ADAPTIVE LOGIC: Check previous frame's NIS
                    if obj.Tracks(i).LastNIS > obj.ManeuverThresh
                        Q_sel = obj.Q_Man_Mat;
                    else
                        Q_sel = obj.Q_CV_Mat;
                    end
                    
                    P = obj.Tracks(i).Covariance;
                    obj.Tracks(i).Covariance = obj.F * P * obj.F' + Q_sel;
                end
            end
            
            % B. Potential Tracks (Assume Cruise)
            nP = length(obj.PotentialTracks);
            if nP > 0
                allStatesP = [obj.PotentialTracks.State];
                allStatesP = obj.F * allStatesP;
                
                for i = 1:nP
                    obj.PotentialTracks(i).State = allStatesP(:, i);
                    P = obj.PotentialTracks(i).Covariance;
                    obj.PotentialTracks(i).Covariance = obj.F * P * obj.F' + obj.Q_CV_Mat;
                end
            end
            
            % --- 2. DYNAMIC MERGE THRESHOLD ---
            if ~isempty(detections)
                R_curr = detections{1}.MeasurementNoise;
                sigma_max = sqrt(max(diag(R_curr)));
                
                obj.MergeThresholdSq = (sigma_max)^2; 
            end
            
            % --- 3. ASSOCIATION (CONFIRMED) ---
            [matches, unassignedT, unassignedD] = ...
                obj.global_association(obj.Tracks, detections, obj.Gate, dt);
            
            % A. Update Matched Tracks
            for k = 1:size(matches, 1)
                tIdx = matches(k,1); dIdx = matches(k,2);
                
                % Store NIS for next frame's adaptation
                currentNIS = CubatureKalmanFilter.calcNIS(obj.Tracks(tIdx).State, ...
                                                          obj.Tracks(tIdx).Covariance, ...
                                                          detections{dIdx}.Measurement, ...
                                                          detections{dIdx}.MeasurementNoise);
                obj.Tracks(tIdx).LastNIS = currentNIS;

                [obj.Tracks(tIdx).State, obj.Tracks(tIdx).Covariance] = ...
                    CubatureKalmanFilter.correct(obj.Tracks(tIdx).State, ...
                                                 obj.Tracks(tIdx).Covariance, ...
                                                 detections{dIdx}.Measurement, ...
                                                 detections{dIdx}.MeasurementNoise);
                obj.Tracks(tIdx).Misses = 0;
            end
            
            % B. Increment Misses
            if ~isempty(unassignedT)
                for i = 1:length(unassignedT)
                    idx = unassignedT(i);
                    obj.Tracks(idx).Misses = obj.Tracks(idx).Misses + 1;
                    obj.Tracks(idx).LastNIS = obj.Tracks(idx).LastNIS * 0.9; % Decay NIS
                end
            end

            % --- 4. ASSOCIATION (POTENTIAL) ---
            leftoverDets = detections(unassignedD);
            [potMatches, unassignedPot, unAssDetsFinal] = ...
                obj.global_association(obj.PotentialTracks, leftoverDets, obj.Gate, dt);
            
            % A. Update Matched Potentials
            for k = 1:size(potMatches, 1)
                pIdx = potMatches(k,1); dIdx = potMatches(k,2);
                
                [obj.PotentialTracks(pIdx).State, obj.PotentialTracks(pIdx).Covariance] = ...
                    CubatureKalmanFilter.correct(obj.PotentialTracks(pIdx).State, ...
                                                 obj.PotentialTracks(pIdx).Covariance, ...
                                                 leftoverDets{dIdx}.Measurement, ...
                                                 leftoverDets{dIdx}.MeasurementNoise);
                
                obj.PotentialTracks(pIdx).Misses = 0;
                obj.PotentialTracks(pIdx).Age = obj.PotentialTracks(pIdx).Age + 1;
                
                % Bitwise History
                bits = bitshift(obj.PotentialTracks(pIdx).HitBits, 1);
                obj.PotentialTracks(pIdx).HitBits = bitset(bits, 1);
                
                % Promotion
                if obj.count_set_bits(obj.PotentialTracks(pIdx).HitBits) >= obj.MinHits
                    obj.promote_track(pIdx);
                end
            end
            
            % B. Update Unmatched Potentials
            for i = 1:length(unassignedPot)
                idx = unassignedPot(i);
                obj.PotentialTracks(idx).HitBits = bitshift(obj.PotentialTracks(idx).HitBits, 1);
                obj.PotentialTracks(idx).Misses = obj.PotentialTracks(idx).Misses + 1;
                obj.PotentialTracks(idx).Age = obj.PotentialTracks(idx).Age + 1;
            end
            
            % --- 5. TRACK MANAGEMENT ---
            
            % Prune Potential Tracks
            if ~isempty(obj.PotentialTracks)
                keepMask = false(length(obj.PotentialTracks), 1);
                
                % MinAge is roughly 25% of the window size (give it time to build history)
                minGraceAge = ceil(obj.WindowSize * 0.25); 
                
                for i = 1:length(obj.PotentialTracks)
                     m = obj.PotentialTracks(i).Misses;
                     age = obj.PotentialTracks(i).Age;
                     
                     if m <= obj.PotentialMaxMisses
                         % Keep if Young OR Density is Good
                         if age < minGraceAge
                             keepMask(i) = true;
                         else
                             % Only check density if mature
                             dens = obj.count_set_bits(obj.PotentialTracks(i).HitBits) / obj.WindowSize;
                             if dens >= 0.2
                                keepMask(i) = true;
                             end
                         end
                     end
                end
                obj.PotentialTracks = obj.PotentialTracks(keepMask);
            end
            
            % Create New Potentials
            for i = 1:length(unAssDetsFinal)
                obj.create_potential(leftoverDets{unAssDetsFinal(i)});
            end
            
            % Prune Confirmed Tracks
            if ~isempty(obj.Tracks)
                misses = [obj.Tracks.Misses];
                maxMisses = ceil(obj.DeleteTime / dt);
                obj.Tracks = obj.Tracks(misses <= maxMisses);
            end

            % Merge Duplicates
            obj.merge_duplicate_tracks();
        end
        
        function [matches, unT, unD] = global_association(obj, tracks, dets, gate, dt)
            nT = length(tracks); nD = length(dets);
            if nT == 0 || nD == 0
                matches=[]; unT=(1:nT)'; unD=(1:nD)'; return; 
            end
            
            % 1. Vectorized Euclidean Gating
            allTStates = [tracks.State];
            tPos = allTStates([1,3,5], :); 
            dPos = zeros(3, nD);
            for i=1:nD, dPos(:,i) = dets{i}.Measurement; end
            
            coarseGateSq = (obj.MaxSpeed * dt * 2.5)^2;
            
            % Fast Distance Matrix
            DistSq = (tPos(1,:)' - dPos(1,:)).^2 + ...
                     (tPos(2,:)' - dPos(2,:)).^2 + ...
                     (tPos(3,:)' - dPos(3,:)).^2;
                 
            % 2. Fine Gating
            costMat = inf(nT, nD);
            [candT, candD] = find(DistSq < coarseGateSq);
            
            for k = 1:length(candT)
                t = candT(k); d = candD(k);
                nis = CubatureKalmanFilter.calcNIS(tracks(t).State, ...
                                                   tracks(t).Covariance, ...
                                                   dets{d}.Measurement, ...
                                                   dets{d}.MeasurementNoise);
                if nis < gate
                    costMat(t,d) = nis;
                end
            end
            
            if all(isinf(costMat(:)))
                 matches=[]; unT=(1:nT)'; unD=(1:nD)'; return;
            end
            
            % 3. Hungarian Assignment
            [matches, unT, unD] = matchpairs(costMat, gate);
        end
        
        function promote_track(obj, pIdx)
            pot = obj.PotentialTracks(pIdx);
            newTrk.ID = obj.NextTrackID;
            newTrk.State = pot.State; 
            newTrk.Covariance = pot.Covariance;
            newTrk.Misses = 0; 
            newTrk.Status = 'Confirmed'; 
            newTrk.HitBits = 0; 
            newTrk.LastNIS = 0; 
            
            obj.Tracks(end+1) = newTrk;
            obj.NextTrackID = obj.NextTrackID + 1;
            obj.PotentialTracks(pIdx).Misses = 9999; % Mark for deletion
        end
        
        function create_potential(obj, det)
            if length(obj.PotentialTracks) >= obj.MaxPotentialTracks, return; end
            meas = det.Measurement;
            R_diag = diag(det.MeasurementNoise); 
            vVar = (obj.MaxSpeed / 3.0)^2; 
            
            newPot.State = [meas(1); 0; meas(2); 0; meas(3); 0];
            newPot.Covariance = diag([R_diag(1), vVar, R_diag(2), vVar, R_diag(3), vVar]);
            newPot.Misses = 0;
            newPot.HitBits = uint64(1); 
            newPot.Age = 1; % Initialize Age
            
            obj.PotentialTracks(end+1) = newPot;
        end
        
        function merge_duplicate_tracks(obj)
             n = length(obj.Tracks);
             if n < 2, return; end
             allStates = [obj.Tracks.State];
             positions = allStates([1,3,5], :);
             keepMask = true(n, 1);
             
             for i=1:n
                 if ~keepMask(i), continue; end
                 diffs = positions(:, (i+1):n) - positions(:, i);
                 distsSq = sum(diffs.^2, 1);
                 dups = find(distsSq < obj.MergeThresholdSq);
                 keepMask(i + dups) = false;
             end
             obj.Tracks = obj.Tracks(keepMask);
        end
        
        function count = count_set_bits(obj, bitmask)
            mask = bitshift(uint64(1), obj.WindowSize) - 1;
            maskedBits = bitand(bitmask, mask);
            count = sum(bitget(maskedBits, 1:obj.WindowSize));
        end
    end
end