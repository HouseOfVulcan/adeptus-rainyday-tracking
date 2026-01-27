classdef PerformanceEvaluator < handle
    % PERFORMANCE EVALUATOR (Updated: Stores Detections for Replay)
    
    properties
        History         
        Config          
    end
    
    methods
        function obj = PerformanceEvaluator()
            % Initialize History with Detections field
            obj.History = struct('Time', {}, 'Truths', {}, 'Tracks', {}, 'Detections', {});
            
            obj.Config.AssociationGate = 200.0; 
            obj.Config.QualityTolerance = 100.0;
            obj.Config.OSPA_c = 150.0; 
            obj.Config.OSPA_p = 2;     
        end
        
        function record_step(obj, time, truths, tracks, detections)
            stepData.Time       = time;
            stepData.Truths     = truths;
            stepData.Tracks     = tracks;
            stepData.Detections = detections; % <--- Storing the Green Dots
            obj.History(end+1)  = stepData;
        end

        function evaluate(obj, dt)
            % (Keep existing evaluate code unchanged)
            % ...
            
            % [Standard Printouts]
            % ...
            
            % [End of evaluate function]
             fprintf('\n[ EXECUTIVE SUMMARY (dt=%.3fs) ]\n', dt);
            % ... (Full previous print logic)
            % Use the previous code for the rest of this function
            
             % --- 1. CONSOLIDATE DATA ---
            [truthMap, trackMap] = obj.consolidate_histories();
            truthIDs = cell2mat(keys(truthMap));
            
            % --- 2. CALCULATE OSPA (The New Metric) ---
            [avgOSPA, ~, ~] = obj.calculate_ospa_metric();
            
            % --- 3. CALCULATE MOTA ---
            sysMetrics = obj.calculate_mota_metrics();
            
            % --- 4. ROBUST MATCHING ---
            usedTrackIDs = [];
            truthReport = containers.Map('KeyType', 'double', 'ValueType', 'any');
            
            allCandidates = [];
            for i = 1:length(truthIDs)
                tID = truthIDs(i);
                tData = truthMap(tID);
                cands = obj.find_robust_matches(tData, trackMap, dt);
                for c = 1:length(cands)
                    cands(c).TruthID = tID;
                    allCandidates = [allCandidates; cands(c)]; %#ok<AGROW>
                end
            end
            
            if ~isempty(allCandidates)
                [~, sortIdx] = sort([allCandidates.MatchedDuration], 'descend');
                allCandidates = allCandidates(sortIdx);
            end
            
            for k = 1:length(allCandidates)
                cand = allCandidates(k);
                if ismember(cand.TrackID, usedTrackIDs), continue; end
                usedTrackIDs(end+1) = cand.TrackID; %#ok<AGROW>
                if ~isKey(truthReport, cand.TruthID)
                    truthReport(cand.TruthID) = cand;
                else
                    existing = truthReport(cand.TruthID);
                    existing(end+1) = cand; %#ok<AGROW>
                    truthReport(cand.TruthID) = existing;
                end
            end
            
            % Report Card
            reportCard = [];
            for i = 1:length(truthIDs)
                tID = truthIDs(i);
                tData = truthMap(tID);
                truthDuration = tData.Time(end) - tData.Time(1);
                
                if isKey(truthReport, tID)
                    frags = truthReport(tID);
                    coverageMask = false(size(tData.Time));
                    totalPosRMSE=0; totalVelRMSE=0; totalSamples=0;
                    
                    for f = 1:length(frags)
                        [~, ~, iT] = intersect(frags(f).MatchedTimes, tData.Time);
                        coverageMask(iT) = true;
                        n = frags(f).NumSamples;
                        totalPosRMSE = totalPosRMSE + (frags(f).PosRMSE^2 * n);
                        totalVelRMSE = totalVelRMSE + (frags(f).VelRMSE^2 * n);
                        totalSamples = totalSamples + n;
                    end
                    
                    matchedTime = sum(coverageMask) * dt;
                    survPct = min(1.0, matchedTime / max(truthDuration, 1e-6));
                    finalPosRMSE = sqrt(totalPosRMSE / max(1, totalSamples));
                    finalVelRMSE = sqrt(totalVelRMSE / max(1, totalSamples));
                    
                    [~, mainIdx] = max([frags.MatchedDuration]);
                    mainTrackID = frags(mainIdx).TrackID;
                    status = sprintf('T%d', mainTrackID);
                    fragCount = length(frags);
                    
                    [~, sortT] = sort([frags.StartTime]);
                    ttft = frags(sortT(1)).StartTime - tData.Time(1);
                    
                    accFactor = exp(-(finalPosRMSE^2) / (2 * obj.Config.QualityTolerance^2));
                    qualScore = survPct * accFactor * 100;
                    
                    purity = frags(mainIdx).MatchedDuration / ...
                        max(1e-6, (trackMap(mainTrackID).Time(end) - trackMap(mainTrackID).Time(1)));
                    purity = min(100, purity * 100);
                else
                    survPct=0; fragCount=0; ttft=nan; finalPosRMSE=nan; finalVelRMSE=nan; qualScore=0;
                    purity=0; status = 'LOST';
                end
                rc.ID = tID; rc.Surv = survPct*100; rc.Frag = fragCount;
                rc.TTFT = ttft; rc.PosRMSE = finalPosRMSE; rc.VelRMSE = finalVelRMSE;
                rc.Qual = qualScore; rc.Purity = purity; rc.Status = status;
                reportCard = [reportCard; rc]; %#ok<AGROW>
            end
            
            % Print Report
            if ~isempty(reportCard)
                avgQual = mean([reportCard.Qual]);
                validPos = [reportCard.PosRMSE]; validPos = validPos(~isnan(validPos));
                avgPos = mean(validPos);
                avgSurv = mean([reportCard.Surv]);
                detectedCount = sum([reportCard.Surv] > 10.0);
            else
                avgQual=0; avgPos=0; avgSurv=0; detectedCount=0;
            end
            
            fprintf('Overall MOTA Score:     %5.1f%%   (System Health Index)\n', sysMetrics.MOTA_Pct);
            fprintf('OSPA Metric (Avg):      %5.1f m  (Lower is Better)\n', avgOSPA);
            fprintf('Targets Detected:       %d / %d   (Avg Survival: %.1f%%)\n', ...
                    detectedCount, length(truthIDs), avgSurv);
            fprintf('Average Track Quality:  %5.1f%%   (Combined Coverage & Accuracy)\n', avgQual);
            fprintf('Average Position Error: %5.1f m  (RMSE)\n', avgPos);
            fprintf('Total ID Switches:      %d        (Stability Count)\n', sysMetrics.ID_Switches);
            fprintf('Total False Alarms:     %d        (Clutter/Ghost Count)\n', sysMetrics.FalseAlarms);
            fprintf('-----------------------------------------------------------------------------------\n');
            
            fprintf('\n[ DETAILED TARGET BREAKDOWN ]\n');
            fprintf('%-4s | %-6s | %-6s | %-6s | %-9s | %-7s | %-7s | %-8s\n', ...
                    'ID', 'Surv %', 'Frags', 'TTFT', 'Pos RMSE', 'Qual %', 'Purity', 'Main Trk');
            fprintf('-----------------------------------------------------------------------------------\n');
            for k = 1:length(reportCard)
                rc = reportCard(k);
                fprintf('%-4d | %5.1f%% | %-6d | %5.1fs | %8.1fm | %6.1f%% | %6.1f%% | %s\n', ...
                        rc.ID, rc.Surv, rc.Frag, rc.TTFT, rc.PosRMSE, rc.Qual, rc.Purity, rc.Status);
            end
            
            % Ghosts
            confirmedGhosts = [];
            mapKeys = keys(trackMap);
            for i = 1:length(mapKeys)
                tID = mapKeys{i};
                if ~ismember(tID, usedTrackIDs)
                    gData = trackMap(tID);
                    if any(strcmp(gData.Status, 'Confirmed'))
                        confirmedGhosts(end+1) = tID; %#ok<AGROW>
                    end
                end
            end
            fprintf('\n[ GHOST TRACK ANALYSIS ]\n');
            if isempty(confirmedGhosts)
                fprintf('No False Tracks (Ghosts) Detected. Clean Run!\n');
            else
                fprintf('Detected %d False Confirmed Tracks (IDs: %s)\n', ...
                        length(confirmedGhosts), num2str(confirmedGhosts));
            end
            fprintf('===================================================================================\n');
        end
        
        function [avgOSPA, avgLoc, avgCard] = calculate_ospa_metric(obj)
            c = obj.Config.OSPA_c; p = obj.Config.OSPA_p;
            totalOSPA = 0; numSteps = 0;
            for k = 1:length(obj.History)
                step = obj.History(k);
                truths = step.Truths; tracks = step.Tracks;
                validTracks = [];
                for j=1:length(tracks)
                    if strcmp(tracks(j).Status,'Confirmed'), validTracks=[validTracks; tracks(j)]; end
                end
                m = length(truths); n = length(validTracks);
                if m==0 && n==0, continue; end
                numSteps = numSteps + 1;
                
                distMat = zeros(m, n);
                for i=1:m
                    tPos = truths(i).State([1,3,5]);
                    for j=1:n
                        rPos = validTracks(j).State([1,3,5]);
                        d = min(c, norm(tPos - rPos));
                        distMat(i,j) = d^p;
                    end
                end
                [M, ~] = matchpairs(distMat, 1e20);
                if ~isempty(M)
                    matchedSum = sum(distMat(sub2ind(size(distMat), M(:,1), M(:,2))));
                else
                    matchedSum = 0;
                end
                cardPenalty = (c^p) * abs(m-n);
                ospa_val = ((matchedSum + cardPenalty) / max(m,n)) ^ (1/p);
                totalOSPA = totalOSPA + ospa_val;
            end
            if numSteps>0, avgOSPA = totalOSPA/numSteps; else, avgOSPA = c; end
            avgLoc=0; avgCard=0;
        end
        
        function candidates = find_robust_matches(obj, tData, trackMap, dt)
            allTracks = keys(trackMap); candidates = [];
            for k = 1:length(allTracks)
                trkID = allTracks{k}; trkData = trackMap(trkID);
                [~, iT, iTrk] = intersect(tData.Time, trkData.Time);
                if isempty(iT), continue; end
                
                posDiff = tData.Pos(iT, :) - trkData.Pos(iTrk, :);
                distErr = sqrt(sum(posDiff.^2, 2));
                validMask = distErr < obj.Config.AssociationGate;
                numValid = sum(validMask);
                
                if numValid >= 3 
                    c.TrackID = trkID;
                    validDistErr = distErr(validMask);
                    velDiff = tData.Vel(iT, :) - trkData.Vel(iTrk, :);
                    validVelErr = sqrt(sum(velDiff(validMask,:).^2, 2));
                    c.PosRMSE = sqrt(mean(validDistErr.^2));
                    c.VelRMSE = sqrt(mean(validVelErr.^2));
                    c.NumSamples = numValid;
                    c.MatchedDuration = numValid * dt; % Uses Passed DT
                    c.StartTime = trkData.Time(1);
                    c.MatchedTimes = tData.Time(iT(validMask));
                    candidates = [candidates; c]; %#ok<AGROW>
                end
            end
        end
        
        function metrics = calculate_mota_metrics(obj)
            % (Standard MOTA logic - no changes needed)
            totalTruths = 0; totalMisses = 0; totalFalseAlarms = 0; totalIDSwitches = 0;
            lastFrameMatches = containers.Map('KeyType', 'double', 'ValueType', 'double');
            for k = 1:length(obj.History)
                step = obj.History(k);
                truths = step.Truths; tracks = step.Tracks;
                validTracks = [];
                for j=1:length(tracks), if strcmp(tracks(j).Status,'Confirmed'), validTracks=[validTracks; tracks(j)]; end; end
                
                numTruths = length(truths); numTracks = length(validTracks);
                totalTruths = totalTruths + numTruths;
                if numTruths==0 && numTracks==0, continue; end
                
                costMatrix = inf(numTruths, numTracks);
                for t=1:numTruths
                    tPos = truths(t).State([1,3,5]);
                    for r=1:numTracks
                        if norm(tPos - validTracks(r).State([1,3,5])) < obj.Config.AssociationGate
                            costMatrix(t,r) = norm(tPos - validTracks(r).State([1,3,5]));
                        end
                    end
                end
                
                currentFrameMatches = containers.Map('KeyType','double','ValueType','double');
                matches = 0;
                while true
                    [minVal, minIdx] = min(costMatrix(:));
                    if isempty(minVal) || isinf(minVal), break; end
                    [row, col] = ind2sub(size(costMatrix), minIdx);
                    tID = truths(row).ID; rID = validTracks(col).ID;
                    matches = matches + 1;
                    currentFrameMatches(tID) = rID;
                    if isKey(lastFrameMatches, tID) && lastFrameMatches(tID) ~= rID
                        totalIDSwitches = totalIDSwitches + 1;
                    end
                    costMatrix(row,:) = inf; costMatrix(:,col) = inf;
                end
                totalMisses = totalMisses + (numTruths - matches);
                totalFalseAlarms = totalFalseAlarms + (numTracks - matches);
                lastFrameMatches = currentFrameMatches;
            end
            if totalTruths>0, mota = 1 - (totalMisses+totalFalseAlarms+totalIDSwitches)/totalTruths; else, mota=0; end
            metrics.MOTA_Pct = mota*100; metrics.Misses = totalMisses;
            metrics.FalseAlarms = totalFalseAlarms; metrics.ID_Switches = totalIDSwitches;
        end
        
        function [truthMap, trackMap] = consolidate_histories(obj)
            % Helper to organize data by ID rather than by Time
            truthMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            trackMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for k = 1:length(obj.History)
                tStep = obj.History(k); cleanTime = round(tStep.Time, 4);
                for i = 1:length(tStep.Truths)
                    c = tStep.Truths(i);
                    if ~isKey(truthMap, c.ID), s.Time=[]; s.Pos=[]; s.Vel=[]; truthMap(c.ID)=s; end
                    s = truthMap(c.ID); s.Time(end+1,1)=cleanTime; s.Pos(end+1,:)=c.State([1,3,5])'; s.Vel(end+1,:)=c.State([2,4,6])'; truthMap(c.ID)=s;
                end
                for i = 1:length(tStep.Tracks)
                    c = tStep.Tracks(i);
                    if ~strcmp(c.Status, 'Confirmed'), continue; end
                    if ~isKey(trackMap, c.ID), s.Time=[]; s.Pos=[]; s.Vel=[]; s.Status={}; trackMap(c.ID)=s; end
                    s = trackMap(c.ID); s.Time(end+1,1)=cleanTime; s.Pos(end+1,:)=c.State([1,3,5])'; s.Vel(end+1,:)=c.State([2,4,6])'; s.Status{end+1,1}=c.Status; trackMap(c.ID)=s;
                end
            end
        end
    end
end