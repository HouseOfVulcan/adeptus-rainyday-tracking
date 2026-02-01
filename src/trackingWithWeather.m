
function trackingWithWeather(configName)
% trackingWithWeather: Top-level driver for tracking + degradation experiments.
%
% PURPOSE
%   Compare tracking performance across multiple tracker types (GNN, TOMHT, JPDA)
%   and motion models (CV, IMM) under ideal and degraded (rainy) conditions.
%
% WORKFLOW
%   1. Load configuration (scenario parameters)
%   2. Create scenario and generate/load detections
%   3. Configure tracker parameters based on degradation mode
%   4. Run all tracker combinations and collect metrics
%
% BASELINE SOURCE
%   Adapted from MathWorks example:
%   https://www.mathworks.com/help/fusion/ug/tracking-closely-spaced-targets-under-ambiguity.html
% 
% NOTES FOR TEAM
%   - This is the top-level driver script. All "real work" happens in helpers:
%       helperCreateScenario3D   -> scenario definition (truth + radar)
%       helperRunDetections      -> generates detection logs
%       helperRunTracker         -> runs a tracker and produces metrics
%       initCVFilter/initIMMFilter -> defines the filter used per track
%   - If something looks "off", most debugging starts in helperRunDetections
%     (time stamps, measurement noise, clutter injection, etc.).

arguments
    configName (1,1) string = "default"
end

clc; close all;

%% Setup
addProjectPaths();

cfg = load_config("default");
scenarioMode  = cfg.scenario.mode;        % "2D" or "3D" from config
numTargets    = cfg.scenario.num_targets; % from config
sceneDuration = cfg.scenario.duration_s;  % from config

%% Scenario Configurations
enableDegradation = true; % false = IdEAL, true = RAINY

fprintf("\n==============================\n");
fprintf(" RUN START | Scenario = %s\n", ternary(enableDegradation,"RAINY","IDEAL"));
fprintf("==============================\n\n");

%% Generate/Load Detections
% 
% dataLog is the key output:
%   dataLog.Time        -> time stamps per scan/update
%   dataLog.Truth       -> ground truth platform states
%   dataLog.Detections  -> detections per scan (cell array of objectDetections)
%
% IMPORTANT:
%   helperRunDetections is where "RAINY" degradation is injected:
%     - detection dropouts (effective Pd)
%     - inflated measurement noise
%     - added clutter (false alarms)

useSavedDataLog = false;  % true = load, false = generate+save

dataLogFile = fullfile(pwd, "cache", "myRun1.mat");
dataLogDir  = fileparts(dataLogFile);

% Create directory if it doesn't exist
if ~exist(dataLogDir, "dir")
    mkdir(dataLogDir);
end

if useSavedDataLog
    load(dataLogFile, "dataLog");
    fprintf("[INFO] Loaded dataLog from %s\n", dataLogFile);
else
    if scenarioMode == "3D"
        scenario = createScenario3D( ...
            "NumTargets", numTargets, ...
            "SceneDuration", sceneDuration);
    else
        error("Only 3D scenarios are currently supported.");
    end

    dataLog = runDetections(scenario, enableDegradation);

    save(dataLogFile, "dataLog", "-v7.3");
    fprintf("[INFO] Saved dataLog to %s\n", dataLogFile);
end

plotInitialScenario(dataLog);

%% Quick stats on detection count per scan
% Helps confirm degradation is actually happening:
%   IDEAL -> generally stable count close to number of targets
%   RAINY -> increased variability (clutter + dropouts)
nPerScan = cellfun(@numel, dataLog.Detections);
fprintf("Detections/scan stats: min=%g, mean=%.2f, max=%g\n", ...
    min(nPerScan), mean(nPerScan), max(nPerScan));


%% Global tracker knobs (shared baseline assumptions)
% These parameters influence track initiation/maintenance and association.
%
% MaxNumTracks:
%   Upper bound on how many tracks can exist at once (prevents blow-up)
%
% Volume:
%   Assumed surveillance volume (m^3). Used with FAR and Beta.
%   Bigger volume + same FAR implies fewer false alarms per unit volume.
%
% Beta:
%   New track "birth" intensity / expected new target rate (model-dependent).
%   Too low -> tracker refuses to initiate tracks.
%   Too high -> tracker spawns too many tracks (especially under clutter).
%
% pd:
%   DetectionProbability - assumed probability a target is detected each scan.
%   If this is far from reality, track scoring can behave poorly.
numTracks = 500;
volume       = 1e9;
beta      = 1e-14;

if (enableDegradation) 
    pd = 0.7;
else 
    pd = 0.9; 
end

%% Scenario-dependent knobs (gating + clutter assumptions)
% gate:
%   Association threshold. Larger = more permissive; smaller = stricter.
%   Under heavy clutter, a too-large gate increases wrong associations.
%   Under low SNR / noisy measurement, too-small gate loses the target.
%
% farX:
%   FalseAlarmRate assumptions for each tracker type.
%   These are not "truth"; they are the tracker's internal clutter model.
%   If FAR is set too low in a cluttered environment, some trackers can
%   over-trust measurements and diverge or fail initiation logic.


if enableDegradation
    % Under degradation, we assume more false alarms than ideal.
    % (We keep these relatively moderate so association isn't totally starved.)
    gate = 35;
    farGNN  = 1e-2;
    farMHT  = 1e-3;
    farJPDA = 1e-3;

    confirmThresh = 30;
    deleteThresh = -15;
    tomhtThresh = [0.2, 5, 5] * gate;
    maxBranches = 3;
    
    betaJPDA    = 1e-11;        % encourage initiation (tune 1e-12..1e-10)
    gateJPDA    = gate + 10;    % looser gate helps association under clutter
    timeTolJPDA = 0.05;         % detections already snapped to scan times
    numTracksJPDA = 200;        % cap so it doesn't explode

else
    gate = 45;
    farGNN  = 1e-6;
    farMHT  = 1e-6;
    farJPDA = 1e-6;

    confirmThresh = 20;
    deleteThresh = -5;
    tomhtThresh = [0.2, 1, 1] * gate;
    maxBranches = 5;

    betaJPDA    = beta;
    gateJPDA    = gate;
    timeTolJPDA = 0.05;
    numTracksJPDA = numTracks;
end

fprintf("[Config] Scenario=%s | gate=%.1f | farGNN=%.2e | farMHT=%.2e | farJPDA=%.2e | pd=%.2f | volume=%.2e | beta=%.2e\n", ...
    ternary(enableDegradation,"RAINY","IDEAL"), gate, farGNN, farMHT, farJPDA, pd, volume, beta);

%% ============ GNN + CV ============
% GNN: Global Nearest Neighbor (hard assignment)
% CV:  Constant Velocity filter initialization
fprintf("\n============ GNN + CV ============\n");
tracker = trackerGNN( ...
    'FilterInitializationFcn', @initCVFilter, ...
    'MaxNumTracks', numTracks, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', gate, ...
    'TrackLogic', 'Score', ...
    'DetectionProbability', pd, ...
    'FalseAlarmRate', farGNN, ...
    'Volume', volume, ...
    'Beta', beta);

% helperRunTracker does:
%   - step through dataLog detections
%   - call tracker(detections, time)
%   - compute metrics vs truth
[trackSummary, truthSummary, trackMetrics, truthMetrics, timeGNNCV] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

%% ============ GNN + IMM ============
% IMM: Interacting Multiple Model (handles maneuvering better than CV)
fprintf("\n============ GNN + IMM ============\n");
tracker = trackerGNN( ...
    'FilterInitializationFcn', @initIMMFilter, ...
    'MaxNumTracks', numTracks, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', gate, ...
    'TrackLogic', 'Score', ...
    'DetectionProbability', pd, ...
    'FalseAlarmRate', farGNN, ...
    'Volume', volume, ...
    'Beta', beta);

[trackSummary, truthSummary, trackMetrics, truthMetrics, timeGNNIMM] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

%% ============ TOMHT + CV ============
% TOMHT: Track-Oriented MHT (multiple hypotheses; good under ambiguity)
fprintf("\n============ TOMHT + CV ============\n");
tracker = trackerTOMHT( ...
    'FilterInitializationFcn', @initCVFilter, ...
    'MaxNumTracks', numTracks, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', tomhtThresh, ...
    'DetectionProbability', pd, ...
    'FalseAlarmRate', farMHT, ...
    'Volume', volume, ...
    'Beta', beta, ...
    'ConfirmationThreshold', confirmThresh, ...
    'DeletionThreshold', deleteThresh, ...
    'MaxNumHistoryScans', 10, ...
    'MaxNumTrackBranches', maxBranches, ...
    'NScanPruning', 'Hypothesis', ...
    'OutputRepresentation', 'Tracks');

[trackSummary, truthSummary, trackMetrics, truthMetrics, timeTOMHTCV] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

%% ============ TOMHT + IMM ============
fprintf("\n============ TOMHT + IMM ============\n");
tracker = trackerTOMHT( ...
    'FilterInitializationFcn', @initIMMFilter, ...
    'MaxNumTracks', numTracks, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', tomhtThresh, ...
    'DetectionProbability', pd, ...
    'FalseAlarmRate', farMHT, ...
    'Volume', volume, ...
    'Beta', beta, ...
    'ConfirmationThreshold', confirmThresh, ...
    'DeletionThreshold', deleteThresh, ...
    'MaxNumHistoryScans', 10, ...
    'MaxNumTrackBranches', maxBranches, ...
    'NScanPruning', 'Hypothesis', ...
    'OutputRepresentation', 'Tracks');

[trackSummary, truthSummary, trackMetrics, truthMetrics, timeTOMHTIMM] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

%% ============ JPDA + CV ============
% JPDA: Joint Probabilistic Data Association (soft association)
% Integrated logic means it internally manages existence based on association probs.
fprintf("\n============JPDA + CV==================\n");
tracker = trackerJPDA( ...
    'FilterInitializationFcn', @initCVFilter, ...
    'MaxNumTracks', numTracksJPDA, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', gateJPDA, ...
    'TrackLogic', 'Integrated', ...
    'DetectionProbability', pd, ...
    'ClutterDensity', farJPDA/volume, ...
    'NewTargetDensity', betaJPDA, ...
    'TimeTolerance', timeTolJPDA);

[trackSummary, truthSummary, trackMetrics, truthMetrics, timeJPDACV] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

%% ============ JPDA + IMM ============
fprintf("\n============JPDA + IMM==================\n");
tracker = trackerJPDA( ...
    'FilterInitializationFcn', @initIMMFilter, ...
    'MaxNumTracks', numTracksJPDA, ...
    'MaxNumSensors', 1, ...
    'AssignmentThreshold', gateJPDA, ...
    'TrackLogic', 'Integrated', ...
    'DetectionProbability', pd, ...
    'ClutterDensity', farJPDA/volume, ...
    'NewTargetDensity', betaJPDA, ...
    'TimeTolerance', timeTolJPDA);

[trackSummary, truthSummary, trackMetrics, truthMetrics, timeJPDAIMM] = helperRunTracker(dataLog, tracker, false);
disp(trackSummary); disp(truthSummary);
disp(trackMetrics); disp(truthMetrics);

fprintf("\n==============================\n");
fprintf(" RUN END | Scenario = %s\n", ternary(enableDegradation,"RAINY","IDEAL"));
fprintf("==============================\n\n");
end

%% -------- Local helper: ternary --------
% Small utility so we can write:
%   ternary(cond, "A", "B")
% instead of MATLAB's  if/else just for printing.
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function addProjectPaths()
root = pwd;
addpath(genpath(fullfile(root, "src", "helpers")));
addpath(genpath(fullfile(root, "src", "visualization")));
end

