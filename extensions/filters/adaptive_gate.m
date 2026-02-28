function gate = adaptive_gate(trackState, params)
% adaptive_gate  User extension point for adaptive gating.

arguments
    trackState (1,1) struct
    params (1,1) struct
end

gate = NaN;
error("adaptive_gate:notImplemented", "Implement adaptive gating logic in extensions/filters/adaptive_gate.m");
end
