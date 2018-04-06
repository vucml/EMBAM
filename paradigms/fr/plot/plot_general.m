function fig_handle = plot_general(x_vals, y_vals, params)
% PLOT_GENERAL  Generic interface to plot() and errorbar() for
%               analysis-specific plotting functions
%
% fig_handle = plot_general(x_vals, y_vals, params)
%
% INPUTS:
%   x_vals:  a row vector of X values to plot
%   y_vals:  a row vector or matrix of Y values to plot.
%            If a matrix, each row will be plotted as a separate line.
%
%   params:  a structure with any of the following fields:
%
%            errorbars - vector of error bar values to pass to errorbar()
%            xlim - 1x2 vector of X axis limits
%            ylim - 1x2 vector of Y axis limits
%            xlabel - struct with fields: label, size, extras
%              (label - string;
%               size - scalar for font size;
%               extras - cell array of extra inputs to xlabel)
%            ylabel - struct with fields: label, size, extras
%            xtick - vector of X values to tick
%            xticklabel - cell array of strings to label X ticks with
%            ytick - vector of Y values to tick
%            yticklabel - cell array of strings to label Y ticks with
%            fontsize - scalar for figure fontsize
%            title - string for figure title
%            linecolor - color(s) to use for the lines in the plot;
%              any of: blue, green, red, cyan, magenta, yellow,
%              black, white; or their corresponding one-letter codes
%            linetype -  type(s) of line to use in the plot;
%              any of: solid, dotted, dashed, or dashdot; or their
%              corresponding one-letter codes
%            marker - marker type(s) for points in the plot;
%              any of: point, circle, x-mark, plus, star, square, diamond;
%              or their one-letter codes
%            legend - cell array of strings to pass to legend()
%            extras - cell array of extra values for plot()
%            savefig - path (including filename) to save a .fig file to
%            saveps - path (including filename) to save a .eps file to
%
%            NOTE: linecolor, linetype, and marker parameters may
%            be strings or cell arrays of strings.  If cell arrays,
%            the values in each will be cycled through, in order,
%            for each successive line in the plot.

% sanity checking will be handled by plot() or errorbar() and/or
% higher-up functions, so I'm not doing anything here except
% ensuring that params exists and has .xlabel and .ylabel
% substructs, which we assume below.
if ~exist('params', 'var')
    params = struct();
end
if ~isfield(params, 'xlabel')
    params.xlabel = struct();
end
if ~isfield(params, 'ylabel')
    params.ylabel = struct();
end

% generic plotting defaults for the behavioral toolbox:
def.errorbars = [];
def.xlim = []; % allow matlab to handle choosing x limits
def.ylim = []; % allow matlab to handle choosing y limits
def.xlabel = struct('label', '', 'size', 25);
def.xlabel.extras = {};
def.ylabel = struct('label', '', 'size', 25);
def.ylabel.extras = {};
def.xtick = []; % allow matlab to handle choosing x ticks
def.ytick = []; % allow matlab to handle choosing y ticks
def.xticklabel = {};
def.yticklabel = {};
def.fontsize = 25;
def.title = '';
% line options: solid black line, white circle markers, width 3 and
% marker size 10
def.linecolor = 'k';
def.linetype = {'solid', 'dashed', 'dotted', 'dashdot'};
def.marker = {'circle', 'diamond', 'x-mark', 'plus'};

%def.linetype = {'solid', 'solid', 'solid', 'solid'};
%def.marker = {'none', 'none', 'none', 'none'};

def.extras = {'LineWidth', 3, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 10};
def.legend = {};
% save options
def.savefig = false;
def.saveps = false;

% merge default options with user-defined options, allowing users
% to override defaults.  Since merge_structs is not recursive, we
% also merge in the xlabel and ylabel structs.
merged_params = merge_structs(params, def);
merged_params.xlabel = merge_structs(params.xlabel, def.xlabel);
merged_params.ylabel = merge_structs(params.ylabel, def.ylabel);

% create plot()-friendly line format string
% colors and linetypes: now conveniently in English!
colors = {'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', ...
    'black', 'white'};
color_codes = {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'};
lines = {'solid', 'dotted', 'dashdot', 'dashed'};
line_codes = {'-', ':', '-.', '--'};
markers = {'point', 'circle', 'x-mark', 'plus', 'star', 'square', 'diamond'};
marker_codes = {'.', 'o', 'x', '+', '*', 's', 'd'};

% plot, using errorbar() if requested
num_lines = size(y_vals, 1);
fig_handle = NaN(1, num_lines);

for i = 1:num_lines
    if i > 1
        % prevent user from needing to clear the figure just to get a new plot
        hold on
    end
    linecolor = get_plot_code(i, merged_params.linecolor, colors, color_codes);
    linetype = get_plot_code(i, merged_params.linetype, lines, line_codes);
    
    if strcmp(merged_params.marker,'none')
        line_spec = [linecolor, linetype];
    else
        marker = get_plot_code(i, merged_params.marker, markers, marker_codes);
        
        line_spec = [linecolor, linetype, marker];
    end    
        
    if ~isempty(merged_params.errorbars)
        fig_handle(i) = errorbar(x_vals, y_vals(i, :), ...
            merged_params.errorbars(i, :), ...
            line_spec, merged_params.extras{:});
    else
        fig_handle(i) = plot(x_vals, y_vals(i, :), ...
            line_spec, merged_params.extras{:});
    end
    
end

% adjust axes if requested, as well as font sizes and ticks
if ~isempty(merged_params.xlim)
    xlim(merged_params.xlim);
end
if ~isempty(merged_params.ylim)
    ylim(merged_params.ylim)
end
if ~isempty(merged_params.xtick)
    set(gca, 'XTick', merged_params.xtick)
end
if ~isempty(merged_params.ytick)
    set(gca, 'YTick', merged_params.ytick)
end
if ~isempty(merged_params.xticklabel)
    set(gca, 'XTickLabel', merged_params.xticklabel);
end
if ~isempty(merged_params.yticklabel)
    set(gca, 'YTickLabel', merged_params.yticklabel);
end
set(gca, 'FontSize', merged_params.fontsize)

% title and labels
if ~isempty(merged_params.xlabel.label)
    xlabel(merged_params.xlabel.label, 'FontSize', merged_params.xlabel.size, ...
        merged_params.xlabel.extras{:})
end
if ~isempty(merged_params.ylabel.label)
    ylabel(merged_params.ylabel.label, 'FontSize', merged_params.ylabel.size, ...
        merged_params.ylabel.extras{:})
end
if ~isempty(merged_params.title)
    title(merged_params.title)
end

% legend
if ~isempty(merged_params.legend)
    legend(fig_handle, merged_params.legend{:});
end

% save .fig and .eps files if requested
if isstr(merged_params.savefig)
    hgsave(fig_handle, merged_params.savefig);
end
if isstr(merged_params.saveps)
    print('-deps', merged_params.saveps);
end

hold off

%endfunction

function code = get_plot_code(line_number, param, english, codes)
% Helper to translate English words for plot params into
% plot()-friendly codes.  Will cycle through a cell array of
% values in param based on line_number
if iscell(param)
    % if the user passed a cell array of values for linetype,
    % marker, etc., cycle through them based on the line number.
    % Because N modN == 0, use N instead when we get to the last value.
    num_values = length(param);
    ind = mod(line_number, num_values);
    if ind == 0,  ind = num_values; end
    item = param{ind};
else
    item = param;
end

if ismember(item, codes)
    % directly-specified plotting codes take precedence
    code = item;
    return
else
    [tf, loc] = ismember(item, english);
    if loc
        code = codes{loc};
    else
        % assume the user knew what they were doing, passing in some
        % string we weren't aware of; let plot() handle any errors
        code = item;
    end
end
%endfunction
