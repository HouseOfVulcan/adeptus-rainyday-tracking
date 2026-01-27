classdef SimulationPlotter < handle
    % =========================================================================
    % VISUALIZATION ENGINE (FIXED LEGENDS)
    % =========================================================================
    % 1. Truth + Detections
    % 2. Tracker + Detections
    % 3. Truth Clean (No Detections)
    % 4. Tracker Clean (No Detections)
    % =========================================================================
    
    properties
        % Noisy Windows
        hFigTruth, hAxTruth, hTruthObj, hDetsTruth
        hFigTrack, hAxTrack, hTrackObj, hDetsTrack
        
        % Clean Windows
        hFigTruthClean, hAxTruthClean, hTruthObjClean
        hFigTrackClean, hAxTrackClean, hTrackObjClean
        
        Config, ColorPalette
        
        % Legend Handles (Arrays to store the proxy lines for the legend)
        LegTruth, LegTrack, LegTruthClean, LegTrackClean
    end
    
    methods
        function obj = SimulationPlotter(config)
            if nargin < 1, config.Limits = [-5000 15000; -10000 10000; 0 10000]; end
            obj.Config = config;
            obj.ColorPalette = lines(12);
            
            % --- WINDOW 1: TRUTH (NOISY) ---
            obj.hFigTruth = figure('Name', 'TRUTH (With Detections)', 'Color', 'w', 'Position', [50, 400, 500, 400]);
            obj.hAxTruth = axes('Parent', obj.hFigTruth, 'FontSize', 9, 'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTruth, 'Truth + Detections');
            obj.hDetsTruth = animatedline('Parent', obj.hAxTruth, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0.7 0], 'MarkerSize', 4, 'MaximumNumPoints', 2000);
            
            % Initial Legend (Detections)
            hold(obj.hAxTruth, 'on');
            hDet = plot3(obj.hAxTruth, nan, nan, nan, '.', 'Color', [0 0.7 0], 'MarkerSize', 10, 'DisplayName', 'Detections');
            obj.LegTruth = hDet;
            legend(obj.hAxTruth, obj.LegTruth, 'Location', 'northeast', 'AutoUpdate', 'off');
            obj.hTruthObj = struct();
            
            % --- WINDOW 2: TRACKER (NOISY) ---
            obj.hFigTrack = figure('Name', 'TRACKER (With Detections)', 'Color', 'w', 'Position', [560, 400, 500, 400]);
            obj.hAxTrack = axes('Parent', obj.hFigTrack, 'FontSize', 9, 'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTrack, 'Tracker + Detections');
            obj.hDetsTrack = animatedline('Parent', obj.hAxTrack, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0.7 0], 'MarkerSize', 4, 'MaximumNumPoints', 2000);
            
            hold(obj.hAxTrack, 'on');
            hDet = plot3(obj.hAxTrack, nan, nan, nan, '.', 'Color', [0 0.7 0], 'MarkerSize', 10, 'DisplayName', 'Detections');
            obj.LegTrack = hDet;
            legend(obj.hAxTrack, obj.LegTrack, 'Location', 'northeast', 'AutoUpdate', 'off');
            obj.hTrackObj = struct();

            % --- WINDOW 3: TRUTH (CLEAN) ---
            obj.hFigTruthClean = figure('Name', 'TRUTH (Clean)', 'Color', 'w', 'Position', [50, 50, 500, 400]);
            obj.hAxTruthClean = axes('Parent', obj.hFigTruthClean, 'FontSize', 9, 'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTruthClean, 'Truth Only (Clean)');
            obj.LegTruthClean = [];
            legend(obj.hAxTruthClean, 'off');
            obj.hTruthObjClean = struct();
            
            % --- WINDOW 4: TRACKER (CLEAN) ---
            obj.hFigTrackClean = figure('Name', 'TRACKER (Clean)', 'Color', 'w', 'Position', [560, 50, 500, 400]);
            obj.hAxTrackClean = axes('Parent', obj.hFigTrackClean, 'FontSize', 9, 'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTrackClean, 'Tracker Only (Clean)');
            obj.LegTrackClean = [];
            legend(obj.hAxTrackClean, 'off');
            obj.hTrackObjClean = struct();
            
            drawnow;
        end
        
        function update(obj, time, truth_list, detections, track_list)
            % 1. UPDATE DETECTIONS
            if ~isempty(detections)
                pos = cellfun(@(d) d.Measurement, detections, 'UniformOutput', false);
                posMat = [pos{:}]; 
                addpoints(obj.hDetsTruth, posMat(1,:), posMat(2,:), posMat(3,:));
                addpoints(obj.hDetsTrack, posMat(1,:), posMat(2,:), posMat(3,:));
            end
            
            % 2. UPDATE TRUTH (Windows 1 & 3)
            for i = 1:length(truth_list)
                t = truth_list(i); fieldID = sprintf('T%d', t.ID);
                pos = t.State([1,3,5]); c = obj.get_color(t.ID);
                
                % Window 1 (Noisy)
                if ~isfield(obj.hTruthObj, fieldID)
                    obj.hTruthObj.(fieldID) = obj.create_target_graphics(obj.hAxTruth, c);
                    % Add to Legend
                    hL = plot3(obj.hAxTruth, nan, nan, nan, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', sprintf('Truth %d', t.ID));
                    obj.LegTruth(end+1) = hL; 
                    legend(obj.hAxTruth, obj.LegTruth, 'Location', 'northeast');
                end
                obj.update_graphics(obj.hTruthObj.(fieldID), pos);
                
                % Window 3 (Clean)
                if ~isfield(obj.hTruthObjClean, fieldID)
                    obj.hTruthObjClean.(fieldID) = obj.create_target_graphics(obj.hAxTruthClean, c);
                    % Add to Legend
                    hL = plot3(obj.hAxTruthClean, nan, nan, nan, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', sprintf('Truth %d', t.ID));
                    obj.LegTruthClean(end+1) = hL;
                    legend(obj.hAxTruthClean, obj.LegTruthClean, 'Location', 'northeast');
                end
                obj.update_graphics(obj.hTruthObjClean.(fieldID), pos);
            end
            
            % 3. UPDATE TRACKS (Windows 2 & 4)
            numConfirmed = 0;
            for i = 1:length(track_list)
                trk = track_list(i);
                if ~strcmp(trk.Status, 'Confirmed'), continue; end
                numConfirmed = numConfirmed + 1;
                
                fieldID = sprintf('Trk%d', trk.ID); pos = trk.State([1,3,5]); c = obj.get_color(trk.ID);
                
                % Window 2 (Noisy)
                if ~isfield(obj.hTrackObj, fieldID)
                    obj.hTrackObj.(fieldID) = obj.create_track_graphics(obj.hAxTrack, c, trk.ID);
                    % Add to Legend
                    hL = plot3(obj.hAxTrack, nan, nan, nan, 's-', 'Color', c, 'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d', trk.ID));
                    obj.LegTrack(end+1) = hL;
                    legend(obj.hAxTrack, obj.LegTrack, 'Location', 'northeast');
                end
                obj.update_track_graphics(obj.hTrackObj.(fieldID), pos, trk.ID);
                
                % Window 4 (Clean)
                if ~isfield(obj.hTrackObjClean, fieldID)
                    obj.hTrackObjClean.(fieldID) = obj.create_track_graphics(obj.hAxTrackClean, c, trk.ID);
                    % Add to Legend
                    hL = plot3(obj.hAxTrackClean, nan, nan, nan, 's-', 'Color', c, 'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d', trk.ID));
                    obj.LegTrackClean(end+1) = hL;
                    legend(obj.hAxTrackClean, obj.LegTrackClean, 'Location', 'northeast');
                end
                obj.update_track_graphics(obj.hTrackObjClean.(fieldID), pos, trk.ID);
            end
            
            % Titles
            title(obj.hAxTruth, sprintf('Truth (Noisy) | T=%.1f', time));
            title(obj.hAxTruthClean, sprintf('Truth (Clean) | T=%.1f', time));
            title(obj.hAxTrack, sprintf('Tracker (Noisy) | Conf=%d', numConfirmed));
            title(obj.hAxTrackClean, sprintf('Tracker (Clean) | Conf=%d', numConfirmed));
            
            drawnow limitrate;
        end
        
        function plot_final_history(obj, history)
            fprintf('Generating 4-View Final Plots...\n');
            
            truthMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            trackMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            allDetPoints = [];
            
            % 1. Unpack
            for k = 1:length(history)
                step = history(k);
                if ~isempty(step.Detections)
                    if isstruct(step.Detections)
                         posMat = [step.Detections.Measurement]'; 
                         if ~isempty(posMat), allDetPoints = [allDetPoints; posMat]; end
                    elseif iscell(step.Detections)
                        posCells = cellfun(@(d) d.Measurement', step.Detections, 'UniformOutput', false);
                        if ~isempty(posCells), allDetPoints = [allDetPoints; vertcat(posCells{:})]; end
                    end
                end
                for i = 1:length(step.Truths)
                    t = step.Truths(i);
                    if ~isKey(truthMap, t.ID), s.Pos=[]; truthMap(t.ID)=s; end
                    s = truthMap(t.ID); s.Pos(end+1,:) = t.State([1,3,5])'; truthMap(t.ID) = s;
                end
                for i = 1:length(step.Tracks)
                    t = step.Tracks(i);
                    if ~strcmp(t.Status, 'Confirmed'), continue; end
                    if ~isKey(trackMap, t.ID), s.Pos=[]; trackMap(t.ID)=s; end
                    s = trackMap(t.ID); s.Pos(end+1,:) = t.State([1,3,5])'; trackMap(t.ID) = s;
                end
            end
            
            % 2. Helper to Plot Targets WITH Legend Handling
            function draw_targets(axH, withDets)
                cla(axH); hold(axH, 'on'); legendHandles = [];
                
                % Plot Detections
                if withDets && ~isempty(allDetPoints)
                    if size(allDetPoints,1) > 10000
                        dsIdx = 1:2:size(allDetPoints,1);
                        hD = plot3(axH, allDetPoints(dsIdx,1), allDetPoints(dsIdx,2), allDetPoints(dsIdx,3), '.', 'Color', [0 0.8 0], 'MarkerSize', 1, 'DisplayName', 'Detections');
                    else
                        hD = plot3(axH, allDetPoints(:,1), allDetPoints(:,2), allDetPoints(:,3), '.', 'Color', [0 0.8 0], 'MarkerSize', 2, 'DisplayName', 'Detections');
                    end
                    legendHandles = [legendHandles, hD];
                end
                
                % Plot Truth Lines
                ids = cell2mat(keys(truthMap));
                for ii = 1:length(ids)
                    id = ids(ii); d = truthMap(id); c = obj.get_color(id);
                    hL = plot3(axH, d.Pos(:,1), d.Pos(:,2), d.Pos(:,3), '-', 'LineWidth', 1.5, 'Color', c, 'DisplayName', sprintf('Truth %d', id));
                    legendHandles = [legendHandles, hL];
                end
                legend(axH, legendHandles, 'Location', 'northeast');
            end
            
            % 3. Helper to Plot Tracks WITH Legend Handling
            function draw_tracks(axH, withDets)
                cla(axH); hold(axH, 'on'); legendHandles = [];
                
                if withDets && ~isempty(allDetPoints)
                     if size(allDetPoints,1) > 10000
                        dsIdx = 1:2:size(allDetPoints,1);
                        hD = plot3(axH, allDetPoints(dsIdx,1), allDetPoints(dsIdx,2), allDetPoints(dsIdx,3), '.', 'Color', [0 0.8 0], 'MarkerSize', 1, 'DisplayName', 'Detections');
                    else
                        hD = plot3(axH, allDetPoints(:,1), allDetPoints(:,2), allDetPoints(:,3), '.', 'Color', [0 0.8 0], 'MarkerSize', 2, 'DisplayName', 'Detections');
                    end
                    legendHandles = [legendHandles, hD];
                end
                
                ids = cell2mat(keys(trackMap));
                for ii = 1:length(ids)
                    id = ids(ii); d = trackMap(id); c = obj.get_color(id);
                    hL = plot3(axH, d.Pos(:,1), d.Pos(:,2), d.Pos(:,3), '.-', 'LineWidth', 1, 'Color', c, 'MarkerSize', 4, 'DisplayName', sprintf('Track %d', id));
                    legendHandles = [legendHandles, hL];
                end
                legend(axH, legendHandles, 'Location', 'northeast');
            end
            
            % 4. Execute
            draw_targets(obj.hAxTruth, true);      title(obj.hAxTruth, 'Truth + Detections');
            draw_targets(obj.hAxTruthClean, false); title(obj.hAxTruthClean, 'Truth (CLEAN)');
            
            draw_tracks(obj.hAxTrack, true);       title(obj.hAxTrack, 'Tracker + Detections');
            draw_tracks(obj.hAxTrackClean, false);  title(obj.hAxTrackClean, 'Tracker (CLEAN)');
            
            drawnow;
        end
        
        % --- GRAPHICS HELPERS ---
        function gfx = create_target_graphics(obj, ax, c)
            gfx.Line = animatedline('Parent', ax, 'Color', c, 'LineWidth', 1.5);
            gfx.Head = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 5, 'MarkerFaceColor', c);
        end
        
        function update_graphics(~, gfx, pos)
            addpoints(gfx.Line, pos(1), pos(2), pos(3));
            set(gfx.Head, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
        end
        
        function gfx = create_track_graphics(obj, ax, c, id)
            gfx.Line = animatedline('Parent', ax, 'Color', c, 'LineWidth', 1.5);
            gfx.Head = plot3(ax, nan, nan, nan, 's', 'MarkerSize', 6, 'MarkerFaceColor', c);
            gfx.Text = text(ax, nan, nan, nan, sprintf('T%d', id), 'Color', c, 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        function update_track_graphics(~, gfx, pos, id)
            addpoints(gfx.Line, pos(1), pos(2), pos(3));
            set(gfx.Head, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
            set(gfx.Text, 'Position', pos + [0;0;500], 'String', sprintf('T%d', id));
        end
        
        function c = get_color(obj, id)
            idx = mod(id - 1, size(obj.ColorPalette, 1)) + 1;
            c = obj.ColorPalette(idx, :);
        end
        
        function setup_axes(obj, axH, titleStr)
            hold(axH, 'on'); grid(axH, 'on'); axis(axH, 'equal'); view(axH, 3);
            xlim(axH, obj.Config.Limits(1,:)); ylim(axH, obj.Config.Limits(2,:)); zlim(axH, obj.Config.Limits(3,:));
            xlabel(axH, 'X'); ylabel(axH, 'Y'); zlabel(axH, 'Z'); title(axH, titleStr);
        end
    end
end