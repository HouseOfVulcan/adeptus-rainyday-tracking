classdef SimulationPlotter < handle
    % =========================================================================
    % VISUALIZATION ENGINE (HIGH PERFORMANCE)
    % =========================================================================
    %
    % PURPOSE:
    %   Handles all graphical outputs for the simulation in real-time.
    %   - Window 1: Ground Truth + Detections
    %   - Window 2: Tracker Output + Detections
    %
    % PERFORMANCE FEATURES:
    %   - Uses 'animatedline' for efficient path rendering.
    %   - Uses 'limitrate' to decouple simulation speed from rendering framerate.
    %     (The sim runs as fast as possible, updating graphics ~30fps).
    %
    % VISUAL FEATURES:
    %   - Dynamic Legends: Automatically detects new Targets/Tracks and assigns
    %     consistent colors (Red for T1, Blue for T2, etc.).
    %   - History Tails: Visualizes the path history to show curvature.
    % =========================================================================
    
    properties
        hFigTruth, hAxTruth, hTruthObj, hDetsTruth
        hFigTrack, hAxTrack, hTrackObj, hDetsTrack
        Config, ColorPalette, LegendHandlesTruth, LegendHandlesTrack
    end
    
    methods
        function obj = SimulationPlotter(config)
            % CONSTRUCTOR: Initialize Windows
            if nargin < 1, config.Limits = [-5000 15000; -10000 10000; 0 10000]; end
            obj.Config = config;
            obj.ColorPalette = lines(12); % Standard high-contrast colors
            
            % --- WINDOW 1: TRUTH ---
            obj.hFigTruth = figure('Name', 'SCENARIO TRUTH', 'Color', 'w', ...
                                   'Position', [100, 300, 600, 500]);
            obj.hAxTruth = axes('Parent', obj.hFigTruth, 'FontSize', 10, ...
                                'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTruth, 'Ground Truth & Detections');
            
            hold(obj.hAxTruth, 'on');
            hDetT = plot3(obj.hAxTruth, nan, nan, nan, '.', 'Color', [0 0.7 0], ...
                          'MarkerSize', 6, 'DisplayName', 'Detections');
            obj.LegendHandlesTruth = [hDetT];
            legend(obj.hAxTruth, obj.LegendHandlesTruth, 'Location', 'northeast', 'AutoUpdate', 'off');
            
            % Data Containers
            obj.hTruthObj = struct();
            obj.hDetsTruth = animatedline('Parent', obj.hAxTruth, 'LineStyle', 'none', ...
                                          'Marker', '.', 'Color', [0 0.7 0], ...
                                          'MarkerSize', 6, 'MaximumNumPoints', Inf);
            
            % --- WINDOW 2: TRACKER ---
            obj.hFigTrack = figure('Name', 'TRACKER OUTPUT', 'Color', 'w', ...
                                   'Position', [710, 300, 600, 500]);
            obj.hAxTrack = axes('Parent', obj.hFigTrack, 'FontSize', 10, ...
                                'GridAlpha', 0.3);
            obj.setup_axes(obj.hAxTrack, 'Tracker Output');
            
            hold(obj.hAxTrack, 'on');
            hDetTrk = plot3(obj.hAxTrack, nan, nan, nan, '.', 'Color', [0 0.7 0], ...
                            'MarkerSize', 6, 'DisplayName', 'Detections');
            obj.LegendHandlesTrack = [hDetTrk];
            legend(obj.hAxTrack, obj.LegendHandlesTrack, 'Location', 'northeast', 'AutoUpdate', 'off');
            
            % Data Containers
            obj.hTrackObj = struct();
            obj.hDetsTrack = animatedline('Parent', obj.hAxTrack, 'LineStyle', 'none', ...
                                          'Marker', '.', 'Color', [0 0.7 0], ...
                                          'MarkerSize', 6, 'MaximumNumPoints', Inf);
        end
        
        function update(obj, time, truth_list, detections, track_list)
            % 1. UPDATE DETECTIONS
            if ~isempty(detections)
                pos = cellfun(@(d) d.Measurement, detections, 'UniformOutput', false);
                posMat = [pos{:}]; 
                addpoints(obj.hDetsTruth, posMat(1,:), posMat(2,:), posMat(3,:));
                addpoints(obj.hDetsTrack, posMat(1,:), posMat(2,:), posMat(3,:));
            end
            
            % 2. UPDATE TRUTH (Dynamic Legend)
            for i = 1:length(truth_list)
                t = truth_list(i); fieldID = sprintf('T%d', t.ID);
                pos = t.State([1,3,5]); c = obj.get_color(t.ID);
                
                if ~isfield(obj.hTruthObj, fieldID)
                    % Initialize Graphics
                    obj.hTruthObj.(fieldID).Line = animatedline('Parent', obj.hAxTruth, ...
                        'Color', c, 'LineWidth', 1.5);
                    obj.hTruthObj.(fieldID).Head = plot3(obj.hAxTruth, nan, nan, nan, 'o', ...
                        'MarkerSize', 5, 'MarkerFaceColor', c);
                    
                    % Update Legend
                    hLegend = plot3(obj.hAxTruth, nan, nan, nan, '-', 'Color', c, ...
                                    'LineWidth', 1.5, 'DisplayName', sprintf('Target %d', t.ID));
                    obj.LegendHandlesTruth(end+1) = hLegend;
                    legend(obj.hAxTruth, obj.LegendHandlesTruth, 'Location', 'northeast');
                end
                
                addpoints(obj.hTruthObj.(fieldID).Line, pos(1), pos(2), pos(3));
                set(obj.hTruthObj.(fieldID).Head, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
            end
            
            % 3. UPDATE TRACKS (Dynamic Legend)
            numConfirmed = 0;
            for i = 1:length(track_list)
                trk = track_list(i);
                if ~strcmp(trk.Status, 'Confirmed'), continue; end
                numConfirmed = numConfirmed + 1;
                
                fieldID = sprintf('Trk%d', trk.ID); pos = trk.State([1,3,5]); c = obj.get_color(trk.ID);
                
                if ~isfield(obj.hTrackObj, fieldID)
                    % Initialize Graphics
                    obj.hTrackObj.(fieldID).Line = animatedline('Parent', obj.hAxTrack, ...
                        'Color', c, 'LineWidth', 1.5);
                    obj.hTrackObj.(fieldID).Head = plot3(obj.hAxTrack, nan, nan, nan, 's', ...
                        'MarkerSize', 6, 'MarkerFaceColor', c);
                    obj.hTrackObj.(fieldID).Text = text(obj.hAxTrack, nan, nan, nan, '', ...
                        'Color', c, 'FontSize', 8, 'FontWeight', 'bold');
                    
                    % Update Legend
                    hLegend = plot3(obj.hAxTrack, nan, nan, nan, 's-', 'Color', c, ...
                                    'MarkerFaceColor', c, 'DisplayName', sprintf('Track %d', trk.ID));
                    obj.LegendHandlesTrack(end+1) = hLegend;
                    legend(obj.hAxTrack, obj.LegendHandlesTrack, 'Location', 'northeast');
                end
                
                addpoints(obj.hTrackObj.(fieldID).Line, pos(1), pos(2), pos(3));
                set(obj.hTrackObj.(fieldID).Head, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
                set(obj.hTrackObj.(fieldID).Text, 'Position', pos + [0;0;500], 'String', sprintf('T%d', trk.ID));
            end
            
            title(obj.hAxTruth, sprintf('Time: %.1f s | Truths: %d', time, length(truth_list)));
            title(obj.hAxTrack, sprintf('Time: %.1f s | Confirmed Tracks: %d', time, numConfirmed));
            
            % PERFORMANCE FIX: Restore limitrate for high-speed simulation
            drawnow limitrate;
        end
        
        function c = get_color(obj, id)
            idx = mod(id - 1, size(obj.ColorPalette, 1)) + 1;
            c = obj.ColorPalette(idx, :);
        end
        
        function setup_axes(obj, axH, titleStr)
            hold(axH, 'on'); grid(axH, 'on'); axis(axH, 'equal'); view(axH, 3);
            xlim(axH, obj.Config.Limits(1,:)); 
            ylim(axH, obj.Config.Limits(2,:)); 
            zlim(axH, obj.Config.Limits(3,:));
            xlabel(axH, 'X'); ylabel(axH, 'Y'); zlabel(axH, 'Z'); 
            title(axH, titleStr);
        end
    end
end