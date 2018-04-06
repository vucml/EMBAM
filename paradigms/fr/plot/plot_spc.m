function fig_handle = plot_spc(spcs, params)
%  PLOT_SPC  Plots a matrix of within-subject SPCs
%
%  fig_handle = plot_spc(spcs, params)
% 
%  INPUTS:
%    spcs:  
%      a matrix, or cell array of matrices, whose rows contain an SPC
%      curve for each subject, as output by spc(): a position (i,j)
%      should contain the probability that subject i recalled serial
%      position j.
% 
%    params:  
%      a structure specifying plot options with any of the following
%      fields:
%               
%        cols - vector of column indices to plot
% 
%      Other params will be passed to plot_general; see its docstring
%      for more information.
% 
%  OUTPUTS:
%  fig_handle:  a handle to the figure once it is plotted
% 

% brief sanity checks
if ~exist('params', 'var')
  params = struct();
end
if ~iscell(spcs)
  spcs = {spcs};
end

% size checks
[num_rows, num_cols] = size(spcs{1});
for i = 1:length(spcs)
  assert(size(spcs{i}, 1) == num_rows && ...
	 size(spcs{i}, 2) == num_cols, ...
	 'spcs matrices must all be the same size');
end

% we assume below that params has .xlabel and .ylabel substructs
if ~isfield(params, 'xlabel')
  params.xlabel = struct();
end
if ~isfield(params, 'ylabel')
  params.ylabel = struct();
end

  
% construct plotting default parameters  
def.errorbars = [];
def.cols = [1 : num_cols];
def.xlim = [0 num_cols+1];
def.ylim = [0 1];
def.xlabel = struct('label', 'Serial Position', 'size', 25);
def.ylabel = struct('label', 'Recall Probability', ...
		    'size', 25);
def.xtick = [def.xlim(1)+1 : 2 : def.xlim(end)];

% merge defaults with user-defined options.  Since merge_structs is
% not recursive, we also merge in the xlabel and ylabel sub-structs.
merged_params = merge_structs(params, def);
merged_params.xlabel = merge_structs(params.xlabel, def.xlabel);
merged_params.ylabel = merge_structs(params.ylabel, def.ylabel);

% prepare plot arguments for plot_general
x_vals = merged_params.cols;

% y_vals should have one row for each input spc matrix
y_vals = NaN(length(spcs), size(x_vals, 2));
for i = 1:length(spcs)
  M = spcs{i};
  y_vals(i, :) = nanmean(M(:, merged_params.cols), 1);
end

% hand off to plot_general
fig_handle = plot_general(x_vals, y_vals, merged_params);
