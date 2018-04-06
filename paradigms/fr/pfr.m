function [p_recalls] = pfr(recalls_matrix, subjects, list_length, ...
			   rec_mask, pres_mask)
%PFR   Probability of first recall.
%
% Computes probability of recall by serial position for the
% first output position.  Wraps spc.m.
%
%  [p_recalls] = pfr(recalls_matrix, subjects, list_length, ...
%                    mask, condition_mask)
%
%  INPUTS:
%  recalls_matrix:
%
%        subjects:
%
%     list_length:
%
%            mask: ???
%
%  condition_mask: ???
%
%  OUTPUTS:
%       p_recalls:

% sanity checks:
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls_matrix, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
end

if ~exist('rec_mask', 'var')
  % create standard clean recalls mask if none was given
  rec_mask = make_clean_recalls_mask2d(recalls_matrix);
elseif size(rec_mask) ~= size(recalls_matrix)
  error('recalls_matrix and mask must have the same shape.')
end

if ~exist('pres_mask', 'var')
  % create a mask for the "standard" condition: all items count.
  pres_mask = true(size(recalls_matrix, 1), list_length);
elseif size(pres_mask, 1) ~= size(recalls_matrix, 1)
  error(['recalls_matrix and pres_mask must have the same number' ...
	 ' of rows'])
elseif size(pres_mask, 2) ~= list_length
  error('pres_mask must have a column for each serial position')
end

% call on the more generalized pfr function
p_recalls = pfr_core(recalls_matrix, subjects, list_length, ...
		rec_mask, pres_mask);