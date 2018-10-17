function [edges, centers] = make_dist_bins_adapt(s, start_edges, min_n)
%MAKE_DIST_BINS_ADAPT   Make distance bins based on minimum bin count.
%
%  [edges, centers] = make_dist_bins_adapt(s, start_edges, min_n)
%
%  INPUTS:
%            s:  vector of similarity values.
%
%  start_edges:  initial edges to start from. Will grow bins
%                starting from this to meet the minimum bin count
%                constraint.
%
%        min_n:  minimum count for each bin.
%
%  OUTPUTS:
%    edges:  edges of all bins.
%
%  centers:  centers, based on the mean of values in each bin.

edges = start_edges(1);
start_ind = 1;
finish_ind = 2;
final = false;
while start_ind < length(start_edges) - 1 && ~final
  if finish_ind > length(start_edges)
    finish_ind = length(start_edges);
    final = true;
  end
  
  count = histc(s, [start_edges(start_ind) start_edges(finish_ind)]);

  if count(1) > min_n
    edges = [edges start_edges(finish_ind)];
    start_ind = finish_ind;
    finish_ind = start_ind + 1;
  else
    finish_ind = finish_ind + 1;
  end
end

centers = NaN(1, length(edges) - 1);
for i = 1:length(edges) - 1
  centers(i) = mean(s(s >= edges(i) & s < edges(i+1)));
end

