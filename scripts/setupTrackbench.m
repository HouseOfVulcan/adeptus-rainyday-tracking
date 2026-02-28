function ctx = setupTrackbench(rootHint)
%setupTrackbench Initialize Trackbench pathing and runtime root context.
%
% Usage:
%   setupTrackbench
%   setupTrackbench("/abs/path/to/trackbench")
%
% This function is intended to be called once at the top of scripts/tests.

arguments
    rootHint (1,1) string = ""
end

if rootHint == ""
    thisFile = mfilename('fullpath');
    rootDir = fileparts(fileparts(thisFile));
else
    rootDir = char(rootHint);
end

srcDir = fullfile(rootDir, "src");
addpath(srcDir);

% Set shared runtime root so package code avoids pwd-based paths.
trackbench.util.rootDir(rootDir);

ctx = struct();
ctx.root = string(rootDir);
ctx.src = string(srcDir);
ctx.config = string(fullfile(rootDir, "config"));
ctx.cache = string(fullfile(rootDir, "cache"));
ctx.outputs = string(fullfile(rootDir, "outputs"));
end
