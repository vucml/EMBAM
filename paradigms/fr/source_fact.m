function source_facts = source_fact(source_matrix, subjects, from_mask, ...
                                    to_mask)
%SOURCE_FACT  Source-based clustering factor.
%
%  source_facts = source_fact(source_matrix, subjects, from_mask, to_mask)
%
%  INPUTS:
%  source_matrix:  a matrix indexing the source of each recalled item.
%
%       subjects:  a column vector which indexes the rows of
%                  recalls_matrix with a subject number (or other
%                  identifier).  That is, the recall trials of subject S
%                  should be located in
%                  recalls_matrix(find(subjects==S), :)
%
%      from_mask:  a logical matrix of the same shape as
%                  recalls_itemnos, which is false at a position (i, j)
%                  if transitions from that position should be
%                  excluded.
%
%        to_mask:  a logical matrix of the same shape as
%                  recalls_itemnos, which is false at a position (i, j)
%                  if transitions to that position should be
%                  excluded.
%  Note: if only one mask is given, it will be used for both the
%  from_mask and to_mask.
%
%  OUTPUTS:
%  source_facts:  a vector of proportion of transitions made to same
%                 source items, one for each subject.

% sanity checks
if ~exist('source_matrix', 'var')
  error('You must pass a source matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('from_mask', 'var')
  error('You must pass a from_mask.')
elseif ~exist('to_mask', 'var')
  %error('You must pass a to_mask.')
  to_mask = from_mask;
end

if ~isequal(size(from_mask), size(source_matrix)) || ...
   ~isequal(size(to_mask), size(source_matrix))
  error('source matrix and masks must have the same shape.')
end
if size(source_matrix, 1) ~= length(subjects)
  error('source matrix must have the same number of rows as subjects.')
end

% source_fact_for_subj will do the work:
source_facts = apply_by_index(@source_fact_for_subj, ...
			    subjects, ...
			    1, ...
			    {from_mask, to_mask, source_matrix});


function subj_source_fact = source_fact_for_subj(from_mask, to_mask, ...
                                                 source_matrix)
% helper for source_fact

step = 1;
params = [];

% call transitions for each trial
num_trials = size(source_matrix, 1);
trial_facts = zeros(num_trials, 1);

for i = 1:num_trials
  % only pass through the transitions for which items were recalled.
  op_facts = transitions(source_matrix(i,:), from_mask(i,:), ...
                         to_mask(i, :), ...
                         @is_same_source, ...
                         step, params); 

  % the main point: source clustering defined as the proportion of
  % transitions to same source out of all possible transitions
  trial_facts(i) = nanmean(op_facts);
end

% aggregate the trial factors to get the subject source factor
subj_source_fact = nanmean(trial_facts);
%endfunction


function d = is_same_source(item1, item2, params)
  % Helper to return whether the transition is to an item of the same or
  % different source. This is the transition function passed to the
  % transitions function.
  d = isequal(item1, item2);
%endfunction