function [lag_crps] = crp(recalls_matrix, subjects, list_length, ...
                          from_mask_rec, to_mask_rec, from_mask_pres, ...
                          to_mask_pres)
%CRP   Conditional response probability as a function of lag (lag-CRP).
%
%  lag_crps = crp(recalls_matrix, subjects, list_length, ...
%                 from_mask_rec, to_mask_rec, from_mask_pres, to_mask_pres)
%
%  INPUTS:
%  recalls_matrix:  a matrix whose elements are serial positions of recalled
%                   items.  The rows of this matrix should represent recalls
%                   made by a single subject on a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls_matrix
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_matrix(find(subjects==S), :)
%
%     list_length:  a scalar indicating the number of serial positions in the
%                   presented lists.  serial positions are assumed to run 
%                   from 1:list_length.
%
%   from_mask_rec:  if given, a logical matrix of the same shape as 
%                   recalls_matrix, which is false at positions (i, j) where
%                   the transition FROM recalls_matrix(i, j) to
%                   recalls_matrix(i, j+1) should be excluded from
%                   the calculation of the CRP.  
%
%     to_mask_rec:  if given, a logical matrix of the same shape as
%                   recalls_matrix, which is false at positions (i,
%                   j) where the transition from
%                   recalls_matrix(i, j-1) TO recalls_matrix(i, j)
%                   should be excluded from the calculation of the CRP.
% 
%                   If neither from_mask_rec nor to_mask_rec is given, a
%                   standard clean recalls mask is used, which
%                   excludes repeats, intrusions and empty cells. If
%                   from_mask_rec is given but to_mask_rec is not,
%                   from_mask_rec will be used for both masks (i.e.,
%                   transitions both to and from masked out elements
%                   will be excluded).
%
%  from_mask_pres:  [trials X list_length] logical matrix. False at
%                   serial positions/items where transitions from that
%                   item should be excluded. If not given, all items are
%                   included.
%
%    to_mask_pres:  [trials X list_length] logical matrix. False at
%                   serial positions/items where transitions to that
%                   item should be excluded. If not given, all items are
%                   included.
%
%  OUTPUTS:
%        lag_crps:  a matrix of lag-CRP values.  Each row contains the values
%                   for one subject.  It has as many columns as there are
%                   possible transitions (i.e., the length of
%                   (-list_length + 1) : (list_length - 1) ).
%                   The center column, corresponding to the "transition of
%                   length 0," is guaranteed to be filled with NaNs.
%
%                   For example, if list_length == 4, a row in lag_crps
%                   has 7 columns, corresponding to the transitions from
%                   -3 to +3:
%                   lag-CRPs:     [ 0.1  0.2  0.3  NaN  0.3  0.1  0.0 ]
%                   transitions:    -3   -2    -1   0    +1   +2   +3

% sanity is always important
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls_matrix, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
elseif ~exist('from_mask_rec', 'var')
  % create standard clean recalls mask if none was given
  from_mask_rec = make_clean_recalls_mask2d(recalls_matrix);
end

if ~exist('to_mask_rec', 'var')
  % assume to_mask_rec should be the same as from_mask_rec (i.e.,
  % transitions to and from the same points should be excluded)
  to_mask_rec = from_mask_rec;
end

if size(from_mask_rec) ~= size(recalls_matrix) | ...
   size(to_mask_rec) ~= size(recalls_matrix)
  error('recalls_matrix and from and to masks must have the same shape.')
end

pres_mask_size = [size(recalls_matrix, 1), list_length];
if ~exist('from_mask_pres', 'var')
  % assume pres_mask should allow all presented items
  from_mask_pres = true(pres_mask_size);
elseif size(from_mask_pres) ~= pres_mask_size
  error(['pres_mask must have the same number of rows as' ...
	 ' recalls_matrix and list_length columns'])
end
if ~exist('to_mask_pres', 'var')
  % for backwards compatibility with previous versions, where there was just
  % a pres_mask input, and this was used for both to and from mask
  to_mask_pres = from_mask_pres;
elseif size(to_mask_pres) ~= pres_mask_size
  error(['pres_mask must have the same number of rows as' ...
	 ' recalls_matrix and list_length columns'])
end

% call on the more generalized form of crp
[lag_crps] = crp_core(recalls_matrix, subjects, list_length, ...
                      from_mask_rec, to_mask_rec, ...
                      from_mask_pres, to_mask_pres);