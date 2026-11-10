function p = project_paths()
%PROJECT_PATHS Resolve all repository paths from the repository root.
% This file lives in <root>/code; every other directory is derived from it.
% Optional environment overrides support raw data or outputs stored outside
% the checkout; no user-specific paths are versioned.
%   p.root       repository root
%   p.code       analysis and reconstruction code (<root>/code)
%   p.data       derived inputs (<root>/data)
%   p.raw        raw .dat files (VSMSI_RAW_DIR, default p.data); GRE .IMA files stay in p.data
%   p.output     all regenerated outputs (VSMSI_OUTPUT_DIR, default <root>/output, git-ignored)
%   p.figures    regenerated figures (<output>/figures)
%   p.results    regenerated numbers*.tex, statistics_results_final.txt/.mat (<output>/results)
%   p.supporting working files of code/supporting_analyses (<output>/supporting)
%   p.reference  committed reference outputs of the paper (<root>/results); never written by the code
%   p.mapvbvd    Siemens raw-data reader (<root>/third_party/mapVBVD)
p.root = fileparts(fileparts(mfilename('fullpath')));
p.code = fullfile(p.root,'code');
p.data = fullfile(p.root,'data');
p.raw = getenv('VSMSI_RAW_DIR');
if isempty(p.raw), p.raw = p.data; end
p.output = getenv('VSMSI_OUTPUT_DIR');
if isempty(p.output), p.output = fullfile(p.root,'output'); end
p.figures = fullfile(p.output,'figures');
p.results = fullfile(p.output,'results');
p.supporting = fullfile(p.output,'supporting');
p.reference = fullfile(p.root,'results');
p.mapvbvd = fullfile(p.root,'third_party','mapVBVD');
if ~isfolder(p.raw), error('vsmsi:MissingRawDir','Raw-data directory does not exist: %s. Fix: unset VSMSI_RAW_DIR or point it to the folder holding the .dat files.',p.raw); end
d = {p.figures,p.results,p.supporting};
for i = 1:numel(d), if ~isfolder(d{i}), mkdir(d{i}); end, end
end
