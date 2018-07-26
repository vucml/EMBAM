function [crps, act_bin, poss_bin] = dist_item_crp(actual, possible, pair_sim, varargin)
%DIST_ITEM_CRP   Conditional response probabilty as a function of distance.
%
%  After calculating possible and actual transitions for individual
%  item pairs, this function will calculate conditional response
%  probabilities for individual distance bins (e.g. based on semantic
%  similarity vales).
%
%  [crps, act_bin, poss_bin] = dist_item_crp(actual, possible, 
%      pair_sim, ...)
%
%  INPUTS:
%  actual - [items x items x groups] numeric array
%      Count of actual transitions that occurred between each pair
%      of items in the pool.
%
%  possible - [items x items x groups] numeric array
%      Count of times that a given transition between item pairs
%      was possible (i.e., the first item was recalled, and the
%      second item was still available for recall).
%
%  pair_sim - [items x items] numeric array
%      Some meaure of similarity between each pair of items (e.g.,
%      semantic similarity).
%
%  OUTPUTS:
%  crps - [groups x bins] numeric array
%      Conditional response probability for each bin.
%
%  act_bin - [groups x bins] numeric array
%      Count of actual transitions in each bin.
%
%  poss_bin - [groups x bins] numeric array
%      Count of possible transitions in each bin.
%
%  OPTIONS:
%  edges - [bins + 1] numeric array - []
%      Edges of distance bins to use. A distance x is included in
%      bin i if edges(i) <= x < edges(i+1).
%
%  percentiles - [bins + 1] numeric array - 0:10:100
%      If edges is not specified, bins will be defined by
%      percentiles of the pairwise similarities.
%
%  mask - [items x items] logical array
%      Pairs of items to include in the analysis.

% options
def.percentiles = 0:10:100;
def.edges = [];
def.mask = true(size(pair_sim));
opt = propval(varargin, def);

% determine edges to use to define bins
if ~isempty(opt.edges)
  edges = opt.edges;
else
  edges = make_dist_bins(pair_sim, opt.percentiles);
end

n_bin = length(edges);
n_group = size(actual, 3);

% find the bin for each pair of items
[~, bin] = histc(squareform(pair_sim), edges);
bin_mat = squareform(bin);

% calculate the CRP for each bin
act_bin = NaN(n_group, n_bin);
poss_bin = NaN(n_group, n_bin);
crps = NaN(n_group, n_bin);
for i = 1:n_group
  sub_actual = actual(:,:,i);
  sub_possible = possible(:,:,i);
  for j = 1:length(edges)
    act_bin(i,j) = nansum(sub_actual(bin_mat == j & opt.mask));
    poss_bin(i,j) = nansum(sub_possible(bin_mat == j & opt.mask));
    
    crps(i,j) = act_bin(i,j) ./ poss_bin(i,j);
  end
end
