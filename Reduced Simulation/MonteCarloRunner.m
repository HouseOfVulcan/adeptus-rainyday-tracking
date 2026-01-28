% MONTECARLORUNNER.m
% =========================================================================
% PURPOSE:
%   The "Batch Processor" for the simulation. It executes the simulation 
%   loop N times silently to gather statistical data. It transforms the 
%   project from a visual demo into a rigorous scientific instrument.
%
% HOW IT CONNECTS:
%   - Called by: MainReducedSim.m (when sim.enableMonteCarlo = true).
%   - Orchestrates: Instantiates fresh copies of TruthGenerator, 
%     DetectionGenerator, and GNNTracker for every single run.
%   - Output: Generates the "Executive Summary" with statistical confidence 
%     intervals (Mean +/- Standard Deviation).
%
% KEY LOGIC:
%   - Fresh Initialization: Crucially resets all random number generators 
%     and tracker history between runs to ensure statistical independence.
%   - "Silent Mode": Suppresses the console spam and plotting overhead to 
%     maximize execution speed.
%   - Advanced Ghost Logic: Distinguishes between "Clutter" (raw false dots), 
%     "Fragments" (valid tracks that broke), and true "Ghosts" (confirmed 
%     tracks on empty space).
%   - Aggregator: Computes deep metrics per-target (Survival, RMSE) across 
%     all runs to find specific weak points in the tracker logic.
% =========================================================================
% MONTE CARLO METRICS & INTERPRETATION GUIDE
% =========================================================================
% This section defines the statistical metrics used to grade the system 
% over multiple runs.
%
% 1. SYSTEM HEALTH METRICS
%   - MOTA (Multi-Object Tracking Accuracy): The "Gold Standard" single 
%      score (0-100%). Accounts for Misses, False Alarms, and ID Switches.
%   - OSPA (Optimal Subpattern Assignment): A distance-based metric (meters).
%      Combines position error and cardinality error (wrong number of targets).
%
% 2. FALSE POSITIVE ANALYSIS
%   - False Alarms: Raw radar detections that do not map to a target. 
%      (High in 'Storm' scenarios due to clutter).
%   - Ghost Tracks: A much more severe failure where the tracker "Confirms" 
%      a trajectory that does not exist. Ideally should be < 0.5 per run.
%
% 3. TARGET FIDELITY (Per-ID Breakdown)
%   - Survival %: Percentage of the target's life it was successfully tracked.
%   - Fragmentation: How many times a single target's track broke into 
%      pieces (e.g., during a crossing event).
%   - Purity: Measures identity consistency. 100% means the target held 
%      the same ID number for the entire duration.
% =========================================================================
classdef MonteCarloRunner < handle
    
    properties
        SimParams       % struct(.dt, .duration, .numRuns)
        SensorConfig    % Weather/Sensor parameters
        ScenarioType    % String name for reporting
        Stats           % System-wide stats arrays
        TargetStats     % Map to store per-target history across runs
    end
    
    methods
        function obj = MonteCarloRunner(simParams, sensorConfig, weatherName)
            obj.SimParams = simParams;
            obj.SensorConfig = sensorConfig;
            obj.ScenarioType = weatherName;
            
            % Initialize System Stats
            N = simParams.numRuns;
            obj.Stats.MOTA           = zeros(N, 1);
            obj.Stats.OSPA           = zeros(N, 1);
            obj.Stats.IDSwitches     = zeros(N, 1);
            obj.Stats.FalseAlarms    = zeros(N, 1); % Clutter Points
            obj.Stats.GhostCounts    = zeros(N, 1); % Actual False Tracks
            obj.Stats.RunTimes       = zeros(N, 1);
            obj.Stats.TargetsDetected= zeros(N, 1);
            obj.Stats.AvgSurvival    = zeros(N, 1);
            obj.Stats.AvgQuality     = zeros(N, 1);
            obj.Stats.AvgPosError    = zeros(N, 1);
            obj.Stats.TotalTargets   = zeros(N, 1);
            
            % Initialize Per-Target Storage
            obj.TargetStats = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
        
        function run(obj)
            fprintf('\n============================================================\n');
            fprintf('STARTING MONTE CARLO SIMULATION\n');
            fprintf('Scenario: %s | Runs: %d | dt: %.2fs\n', ...
                obj.ScenarioType, obj.SimParams.numRuns, obj.SimParams.dt);
            fprintf('============================================================\n\n');
            
            totalTimer = tic;
            
            for k = 1:obj.SimParams.numRuns
                obj.execute_single_run(k);
            end
            
            totalTime = toc(totalTimer);
            fprintf('\nTotal Batch Time: %.2f seconds\n', totalTime);
            
            obj.print_aggregated_report();
        end
        
        function execute_single_run(obj, runIdx)
            % 1. INITIALIZATION
            truthGen = TruthGenerator(obj.SensorConfig); 
            detGen   = DetectionGenerator(obj.SensorConfig);
            tracker  = GNNTracker(obj.SimParams.dt, obj.SensorConfig.TrackerTuning); 
            evaluator = PerformanceEvaluator();
            
            ScenarioGenerator.load_scenario3(truthGen);
            
            % 2. EXECUTION
            numSteps = ceil(obj.SimParams.duration / obj.SimParams.dt);
            runStart = tic;
            
            for i = 1:numSteps
                time = (i-1) * obj.SimParams.dt;
                truth = truthGen.get_states(time);       
                dets  = detGen.step(truth, time, obj.SimParams.dt); 
                tracker.update(dets, obj.SimParams.dt); 
                evaluator.record_step(time, truth, tracker.Tracks, dets);
            end
            
            runDuration = toc(runStart);
            
            % 3. SYSTEM METRICS
            [avgOSPA, ~, ~] = evaluator.calculate_ospa_metric();
            motaMetrics = evaluator.calculate_mota_metrics();
            
            % 4. DETAILED ANALYSIS
            [truthMap, trackMap] = evaluator.consolidate_histories();
            truthIDs = cell2mat(keys(truthMap));
            numTargets = length(truthIDs);
            
            claimedTrackIDs = []; % Tracks that matched at least one target
            
            runSurv = zeros(numTargets, 1);
            runQual = zeros(numTargets, 1);
            runDetected = 0;
            totalPosRMSESum = 0; validRMSECount = 0;
            
            for t = 1:numTargets
                tID = truthIDs(t);
                tData = truthMap(tID);
                truthDuration = tData.Time(end) - tData.Time(1);
                
                if ~isKey(obj.TargetStats, tID)
                    s.Surv = []; s.Frag = []; s.TTFT = []; 
                    s.PosRMSE = []; s.Qual = []; s.Purity = []; s.Detected = [];
                    obj.TargetStats(tID) = s;
                end
                ts = obj.TargetStats(tID);
                
                cands = evaluator.find_robust_matches(tData, trackMap, obj.SimParams.dt);
                
                if ~isempty(cands)
                    % Mark ALL candidates as "Claimed" so they don't count as Ghosts
                    for c = 1:length(cands)
                        claimedTrackIDs = [claimedTrackIDs, cands(c).TrackID]; %#ok<AGROW>
                    end
                    
                    [~, bestIdx] = max([cands.MatchedDuration]);
                    bestMatch = cands(bestIdx);
                    
                    matchedTime = bestMatch.MatchedDuration;
                    survPct = min(1.0, matchedTime / max(truthDuration, 1e-6));
                    finalPosRMSE = bestMatch.PosRMSE;
                    accFactor = exp(-(finalPosRMSE^2) / (2 * evaluator.Config.QualityTolerance^2));
                    qualScore = survPct * accFactor * 100;
                    
                    runSurv(t) = survPct * 100;
                    runQual(t) = qualScore;
                    if survPct > 0.10, runDetected = runDetected + 1; end
                    if ~isnan(finalPosRMSE)
                        totalPosRMSESum = totalPosRMSESum + finalPosRMSE;
                        validRMSECount = validRMSECount + 1;
                    end
                    
                    ts.Surv(end+1)     = survPct * 100;
                    ts.Frag(end+1)     = length(cands);
                    ts.TTFT(end+1)     = bestMatch.StartTime - tData.Time(1);
                    ts.PosRMSE(end+1)  = finalPosRMSE;
                    ts.Qual(end+1)     = qualScore;
                    ts.Purity(end+1)   = (bestMatch.MatchedDuration / max(1e-6, (trackMap(bestMatch.TrackID).Time(end) - trackMap(bestMatch.TrackID).Time(1)))) * 100;
                    ts.Detected(end+1) = (survPct > 0.10);
                else
                    runSurv(t) = 0; runQual(t) = 0;
                    ts.Surv(end+1) = 0; ts.Frag(end+1) = 0; ts.TTFT(end+1) = NaN;
                    ts.PosRMSE(end+1) = NaN; ts.Qual(end+1) = 0; ts.Purity(end+1) = 0;
                    ts.Detected(end+1) = false;
                end
                obj.TargetStats(tID) = ts;
            end
            
            % --- CALCULATE GHOSTS ---
            allTrackIDs = cell2mat(keys(trackMap));
            ghostCount = 0;
            for i = 1:length(allTrackIDs)
                trkID = allTrackIDs(i);
                if ~ismember(trkID, claimedTrackIDs)
                    tData = trackMap(trkID);
                    if any(strcmp(tData.Status, 'Confirmed'))
                        ghostCount = ghostCount + 1;
                    end
                end
            end
            
            % 5. STORE RUN STATS
            obj.Stats.MOTA(runIdx)           = motaMetrics.MOTA_Pct;
            obj.Stats.OSPA(runIdx)           = avgOSPA;
            obj.Stats.IDSwitches(runIdx)     = motaMetrics.ID_Switches;
            obj.Stats.FalseAlarms(runIdx)    = motaMetrics.FalseAlarms;
            obj.Stats.GhostCounts(runIdx)    = ghostCount;
            obj.Stats.RunTimes(runIdx)       = runDuration;
            obj.Stats.TotalTargets(runIdx)   = numTargets;
            obj.Stats.TargetsDetected(runIdx)= runDetected;
            obj.Stats.AvgSurvival(runIdx)    = mean(runSurv);
            obj.Stats.AvgQuality(runIdx)     = mean(runQual);
            if validRMSECount > 0
                obj.Stats.AvgPosError(runIdx) = totalPosRMSESum / validRMSECount;
            else
                obj.Stats.AvgPosError(runIdx) = NaN;
            end
            
            fprintf('Run %02d/%d | MOTA: %5.1f%% | OSPA: %5.1fm | Ghosts: %d\n', ...
                runIdx, obj.SimParams.numRuns, motaMetrics.MOTA_Pct, avgOSPA, ghostCount);
        end
        
        function print_aggregated_report(obj)
            fprintf('\n[ MONTE CARLO SUMMARY (N=%d | dt=%.3fs) ]\n', ...
                obj.SimParams.numRuns, obj.SimParams.dt);
            
            obj.print_stat_row('Overall MOTA Score',    obj.Stats.MOTA,        '%',  true,  '(System Health Index)');
            obj.print_stat_row('OSPA Metric (Avg)',     obj.Stats.OSPA,        'm',  false, '(Lower is Better)');
            
            muDet   = mean(obj.Stats.TargetsDetected);
            maxTgt  = max(obj.Stats.TotalTargets); 
            muSurv  = mean(obj.Stats.AvgSurvival);
            fprintf('Targets Detected       :  %4.1f / %d   (Avg Survival: %4.1f%%)\n', muDet, maxTgt, muSurv);
            
            obj.print_stat_row('Average Track Quality', obj.Stats.AvgQuality,  '%',  true,  '(Combined Coverage & Accuracy)');
            obj.print_stat_row('Average Position Error',obj.Stats.AvgPosError, 'm',  false, '(RMSE)');
            obj.print_stat_row('Total ID Switches',     obj.Stats.IDSwitches,  '',   false, '(Stability Count)');
            obj.print_stat_row('Total False Alarms',    obj.Stats.FalseAlarms, '',   false, '(Clutter Count)');
            fprintf('-----------------------------------------------------------------------------------\n');
            
            % --- DETAILED TARGET BREAKDOWN ---
            fprintf('\n[ DETAILED TARGET BREAKDOWN (Averaged over %d runs) ]\n', obj.SimParams.numRuns);
            fprintf('%-4s | %-6s | %-6s | %-6s | %-9s | %-7s | %-7s | %-8s\n', ...
                    'ID', 'Surv %', 'Frags', 'TTFT', 'Pos RMSE', 'Qual %', 'Purity', 'Det %');
            fprintf('-----------------------------------------------------------------------------------\n');
            
            ids = cell2mat(keys(obj.TargetStats));
            ids = sort(ids);
            
            for i = 1:length(ids)
                id = ids(i);
                s = obj.TargetStats(id);
                muSurv = mean(s.Surv);
                muFrag = mean(s.Frag);
                muTTFT = mean(s.TTFT(~isnan(s.TTFT))); 
                muRMSE = mean(s.PosRMSE(~isnan(s.PosRMSE)));
                muQual = mean(s.Qual);
                muPur  = mean(s.Purity);
                detRate= mean(s.Detected) * 100;
                
                fprintf('%-4d | %5.1f%% | %-6.1f | %5.1fs | %8.1fm | %6.1f%% | %6.1f%% | %5.1f%%\n', ...
                        id, muSurv, muFrag, muTTFT, muRMSE, muQual, muPur, detRate);
            end
            
            % --- GHOST ANALYSIS ---
            fprintf('\n[ GHOST TRACK ANALYSIS ]\n');
            avgGhosts = mean(obj.Stats.GhostCounts);
            
            if avgGhosts < 0.1
                fprintf('Negligible Ghost Tracks (Average < 0.1 per run). Clean!\n');
            else
                fprintf('Average Ghost Tracks (Confirmed Fakes) per Run: %.1f\n', avgGhosts);
            end
            fprintf('===================================================================================\n');
        end
        
        function print_stat_row(~, label, data, unit, ~, desc)
            mu = mean(data, 'omitnan');
            sigma = std(data, 'omitnan');
            fprintf('%-23s: %6.1f%s ± %4.1f%s   %s\n', ...
                label, mu, unit, sigma, unit, desc);
        end
    end
end