function fig_handle = plot_crp(lag_crps, params)
%  PLOT_CRP  Plots a matrix of within-subject CRPs
%
%  fig_handle = plot_crp(lag_crps, params)
% 
%  INPUTS:
%    lag_crps:  
%       a matrix, or cell array of matrices, whose rows contain a CRP
%       curve for each subject, where the column index is equal to
%       (list_length + lag), as output by crp().
% 
%    params:  
%       a structure specifying plot options with any of the
%       following fields:
%          
%          zerocol - index of column corresponding to lag 0
%          maxlag - maximum lag to plot
%          cols - vector of column indices to plot
% 
%          Note that zerocol is required if lag_crps has an
%          even number of columns.  Otherwise, the middle
%          column is assumed to correspond to lag 0.
%
%          Other params will be passed to plot_general; see
%          its docstring for more information.
% 
%  OUTPUTS:
%    fig_handle:  
%       a handle to the figure once it is plotted.  Each lag_crps
%       matrix will correspond to a single line in the plot.
%       
% 

% brief sanity checks
if ~exist('params', 'var')
  params = struct();
end
if ~iscell(lag_crps)
  % assume that a non-cell input is a single matrix
  lag_crps = {lag_crps};
end

% size checks
[num_rows, num_cols] = size(lag_crps{1});
for i = 1:length(lag_crps)
  assert(size(lag_crps{i}, 1) == num_rows && ...
	 size(lag_crps{i}, 2) == num_cols, ...
	 'lag_crps matrices must all be the same size');
end
if ~isfield(params, 'zerocol') 
  if mod(num_cols, 2) == 0
    error('If lag_crps has no central column, you must specify params.zerocol')
  else
    params.zerocol = (num_cols + 1) / 2;
  end
end


% we assume below that params has .xlabel and .ylabel substructs
if ~isfield(params, 'xlabel')
  params.xlabel = struct();
end
if ~isfield(params, 'ylabel')
  params.ylabel = struct();
end
  
% construct plotting default parameters  
if isfield(params, 'maxlag')
  maxlag = params.maxlag;
elseif isfield(params,'cols')
  maxlag = floor(length(params.cols)/2);
else
  maxlag = 5;
end
def.errorbars = [];
def.cols = (params.zerocol - maxlag : params.zerocol + maxlag);
def.xlim = [-(maxlag+1) (maxlag+1)];
def.ylim = [0 0.6];
def.xlabel = struct('label', 'Lag', 'size', 25);
def.ylabel = struct('label', 'Conditional Response Probability', ...
		    'size', 25);
def.xtick = unique([round(linspace(def.xlim(1)+2,0,3)) ...
                    round(linspace(0,def.xlim(2)-2,3))]);

% merge defaults with user-defined options.  Since merge_structs is
% not recursive, we also merge in the xlabel and ylabel sub-structs.
merged_params = merge_structs(params, def);
merged_params.xlabel = merge_structs(params.xlabel, def.xlabel);
merged_params.ylabel = merge_structs(params.ylabel, def.ylabel);


% prepare plot arguments for plot_general
x_vals = merged_params.cols - merged_params.zerocol;

% y_vals should have one row for each input crp matrix
y_vals = NaN(length(lag_crps), size(x_vals, 2));
for i = 1:length(lag_crps)
  M = lag_crps{i};
  y_vals(i, :) = mean(M(:, merged_params.cols), 1, 'omitnan');
end

if ~isempty(merged_params.errorbars)
  % slice out the columns for the error bars we want
  err = NaN(size(y_vals));
  for i = 1:length(lag_crps)
    err(i,:) = merged_params.errorbars{i}(merged_params.cols);
  end
  merged_params.errorbars = err;
end

% hand off to plot_general
fig_handle = plot_general(x_vals, y_vals, merged_params);
