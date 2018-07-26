function init_embam()
%INIT_EMBAM   Initialize EMBAM dependencies.
%
%  init_embam()

proj_dir = fileparts(mfilename('fullpath'));

addpath(fullfile(proj_dir, 'helpers', 'masks'))
addpath(fullfile(proj_dir, 'helpers', 'matrixops'))
addpath(fullfile(proj_dir, 'helpers', 'transitions'))
addpath(fullfile(proj_dir, 'paradigms', 'fr'))
addpath(fullfile(proj_dir, 'paradigms', 'fr', 'plot'))
addpath(fullfile(proj_dir, 'utils'))
