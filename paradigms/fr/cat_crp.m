function [lag_crps, numer, denom] = ...
                    cat_crp(recalls_matrix, pres_catlabels, ...
                            rec_catlabels, subjects, cat_type, ...
                            from_mask_rec, to_mask_rec, from_mask_pres, ...
                            to_mask_pres, lag_type)
%CAT_CRP   Lag-CRP, conditional on category structure.
%
%  Computes (lag) conditional response probabilities from a matrix of
%  recalled serial positions, for either within- or between-category
%  transitions.
%
%  [lag_crps, numer, denom] = cat_crp(recalls_matrix, ...
%                                     pres_catlabels, rec_catlabels, ...
%			              subjects, cat_type, ...
%                                     from_mask_rec, to_mask_rec, ...
%                                     from_mask_pres, to_mask_pres,
%                                     lag_type);
%
%  INPUTS:
%  recalls_matrix:  a matrix whose elements are serial positions of
%                   recalled items.  The rows of this matrix should
%                   represent recalls made by a single subject on a
%                   single trial.
%
%  pres_catlabels:  the category labels of the presented items.
%
%   rec_catlabels:  the category labels of the recalled items.
%
%        subjects:  a column vector which indexes the rows of
%                   recalls_matrix with a subject number (or other
%                   identifier). Recall trials of subject S should be
%                   located in recalls_matrix(find(subjects==S),:)
%
%        cat_type:  if 1, examines same-category transitions; if
%                   anything else (2), examines between-category
%                   transitions.
%
%   from_mask_rec:  if given, a logical matrix of the same shape as 
%                   recalls_matrix, which is false at positions (i, j)
%                   where the transition FROM recalls_matrix(i, j) to
%                   recalls_matrix(i, j+1) should be excluded from the
%                   calculation of the CRP.  
%
%     to_mask_rec:  if given, a logical matrix of the same shape as
%                   recalls_matrix, which is false at positions (i,j)
%                   where the transition from recalls_matrix(i, j-1) TO
%                   recalls_matrix(i, j) should be excluded from the
%                   calculation of the CRP.
% 
%                   If neither mask is given, a standard clean recalls
%                   mask is used, which excludes repeats, intrusions
%                   and empty cells.  If from_mask_rec is given but
%                   to_mask_rec is not, from_mask_rec will be used for
%                   both masks (i.e., transitions both to and from
%                   masked out elements will be excluded).
%
%  from_mask_pres:  a logical matrix of allowable serial positions that
%                   transitions can come FROM.
%
%    to_mask_pres:  a logical matrix of allowable serial positions that
%                   transitions can go TO.
%
%        lag_type:  'cat_pos' (default) or 'serial_pos'
%
%  OUTPUTS:
%  lag_crps:  a matrix of lag-CRP values.  Each row contains the values
%             for one subject.  It has as many columns as there are
%             possible transitions (i.e., the length of
%             (-list_length + 1) : (list_length - 1) ).
%             The center column, corresponding to the "transition of
%             length 0," is guaranteed to be filled with NaNs.
%
%             For example, if list_length == 4, a row in lag_crps
%             has 7 columns, corresponding to the transitions from
%             -3 to +3:
%              lag-CRPs:     [ 0.1  0.2  0.3  NaN  0.3  0.1  0.0 ]
%              transitions:    -3   -2    -1   0    +1   +2   +3

% sanity is always important
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
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
  error('recalls_matrix and masks must have the same shape.')
end

if ~exist('pres_catlabels', 'var')
  error('You must pass in a pres_catlabels matrix');
end

list_length = size(pres_catlabels,2);

if ~exist('to_mask_pres', 'var')
  to_mask_pres = ones(length(subjects),list_length);
end

if ~exist('from_mask_pres', 'var')
  from_mask_pres = ones(length(subjects),list_length);
end

if ~exist('lag_type', 'var')
  lag_type = 'cat_pos';
end

% crp_for_subj is going to do all the real work:
packed_info = apply_by_index(@cat_crp_for_subj, ...
			  subjects, ...
			  1, ...
			  {recalls_matrix, pres_catlabels, ...
		           from_mask_rec, to_mask_rec, ...
                           from_mask_pres, to_mask_pres}, ...
			  cat_type, lag_type);

crp_len = (list_length*2)-1;
lag_crps = packed_info(:,1:crp_len);
numer = packed_info(:,crp_len+1:2*crp_len);
denom = packed_info(:,2*crp_len+1:3*crp_len);

%endfunction

function packed_info = cat_crp_for_subj(recalls, pres_catlabels, ...
             from_mask_rec, to_mask_rec, from_mask_pres, to_mask_pres, ...
             cat_type, lag_type)
  % Helper for crp:
  % calculates the lag-CRP for each possible transition for one subject's
  % recall trials; returns a row vector of CRPs, indexed by
  % (-list_length + index)
  
  % the golden goose: these will store actual and possible transitions for
  % this subject's trials
  actual_transitions = [];
  poss_transitions = [];

  % arguments for conditional_transitions:
  step = 1;   % for now, step of 1 is hard-coded; this may change   
  params =  struct('list_length', size(pres_catlabels,2), ...
		   'cat_type', cat_type);

  switch lag_type
    case 'cat_pos'
      transit_func = @cat_lag;
      condition_func = @cat_transitions;
    case 'serial_pos'
      transit_func = @lag;
      condition_func = @cat_transitions_lag;
    case 'cat_nocond'
      % a transition is possible if that item hasn't been recalled; the
      % category condition is not taken into account when determining
      % possible transitions. Transitions that don't match the
      % category type are not counted. Conditional on some transition
      % having been made, and conditional on a transition of this lag
      % and category type being available, what was the probability of
      % making that transition? Transitions of other category types
      % may be available, and if they made a transition to some
      % other category type when this type and lag were available,
      % that type and lag will be "punished" for not having
      % happened when it could have happened.
      % As opposed to: conditional on having made a transition of
      % this category type, and conditional on this lag being
      % available, what is the probability of making that transition?
      transit_func = @cat_lag;
      condition_func = @possible_transitions;
  end
  
  % conditional_transitions returns the actual and possible transitions for
  % each trial; we simply concatenate them all together and use collect() to
  % gather them up, since we don't care about distinctions between individual
  % trials within a subject
  num_trials = size(recalls, 1);
  for i = 1:num_trials
    params.pres_catlabels = pres_catlabels(i,:);
    params.from_mask_pres = from_mask_pres(i,:);
    params.to_mask_pres = to_mask_pres(i,:);
    [trial_actuals, trial_possibles] = conditional_transitions(...
					recalls(i, :), ...
    	                                from_mask_rec(i, :), to_mask_rec(i, :), ...
					transit_func, condition_func, ...
			                step, params);
    actual_transitions = [actual_transitions, trial_actuals];
    % only grab the cells from trial_possibles if the corresponding
    % element in actual_transitions is a number (not a NaN)
    poss_transitions = [poss_transitions, catcell(trial_possibles(~isnan(trial_actuals)))];
    % poss_transitions = [poss_transitions, catcell(trial_possibles)];
  end

  % the range of all possible transitions depends only on the list length;
  % we want to gather all the subject's actual and possible transitions in
  % this range
  all_possible_transitions = [-params.list_length + 1 : params.list_length - 1];
  actuals_counts = collect(actual_transitions, all_possible_transitions);
  possibles_counts = collect(poss_transitions, all_possible_transitions);

  % and that's it!
  subj_crp = actuals_counts ./ possibles_counts;

  packed_info = [subj_crp, actuals_counts, possibles_counts];
  
%endfunction

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
%endfunction
  

