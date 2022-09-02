function lag_crps = bin_crp(recalls_matrix, subjects, list_length, ...
			    bins, from_mask_rec, to_mask_rec, ...
                            from_mask_pres, to_mask_pres)
%BIN_CRP  Conditional response probability as a function of lag bin.
%
%  lag_crps = bin_crp(recalls_matrix, subjects, list_length, ...
%                     bins, from_mask_rec, to_mask_rec, ...
%                     from_mask_pres, to_mask_pres);
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
%            bins:  a set of bin edges, as you would pass to histc,
%                   in which to bin actual and possible transitions
%                   before calculating conditional response
%                   probabilities.
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
%                   If neither mask is given, a standard clean recalls
%                   mask is used, which excludes repeats, intrusions
%                   and empty cells.  If from_mask is given but
%                   to_mask is not, from_mask will be used for both
%                   masks (i.e., transitions both to and from
%                   masked out elements will be excluded).       
%
%  from_mask_pres:  if given, a logical matrix of the same shape as
%                   the presented items matrix, i.e., 
%                   (num_trials x list_length).  This mask should
%                   be true at positions (t, sp) where an item in
%                   the condition of interest was presented at
%                   serial position sp on trial t; and false
%                   everywhere else.  If NOT given, a blank mask of
%                   this shape is used.
% 
%   to_mask_pres:   if given, logical (complement of from_mask_pres)
%
%  OUTPUTS:
%        lag_crps:  a matrix of lag-CRP values.  Each row contains the values
%                   for one subject.  It has as many columns as there are
%                   bin edges; the value in column i is the CRP for
%                   transitions in the bin bins(i) to bins(i+1).
%                   The last bin contains the CRP for transitions
%                   equal to bins(end).

% sanity is always important
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls_matrix, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
elseif ~exist('bins', 'var')
  error('You must pass some bin edges compatible with histc.')
elseif ~exist('from_mask_rec', 'var')
  % create standard clean recalls mask if none was given
  from_mask_rec = make_clean_recalls_mask2d(recalls_matrix);
end

if ~exist('to_mask_rec', 'var')
  % assume to_mask should be the same as from_mask (i.e.,
  % transitions to and from the same points should be excluded)
  to_mask_rec = from_mask_rec;
end

if size(from_mask_rec) ~= size(recalls_matrix) | ...
   size(to_mask_rec) ~= size(recalls_matrix)
  error('recalls_matrix and masks must have the same shape.')
end

pres_mask_size = [size(recalls_matrix, 1), list_length];
pres_mask = true(pres_mask_size);
if ~exist('from_mask_pres', 'var')
  % if missing, assume pres_mask should allow all presented items
  from_mask_pres = pres_mask;
elseif any(size(from_mask_pres) ~= pres_mask_size)
  error(['from_mask_pres must have the same number of rows as' ...
	 ' recalls_matrix and list_length columns'])
elseif any(size(to_mask_pres) ~= pres_mask_size)
  error(['to_mask_pres must have the same number of rows as' ...
	 ' recalls_matrix and list_length columns'])
end
if ~exist('to_mask_pres', 'var')
    to_mask_pres = pres_mask;
end

% crp_for_subj is going to do all the real work:
lag_crps = apply_by_index(@crp_for_subj, ...
			  subjects, ...
			  1, ...
			  {recalls_matrix, ...
                    from_mask_rec, to_mask_rec, ...
                    from_mask_pres, to_mask_pres}, ...
			  list_length, bins);
%endfunction

function subj_crp = crp_for_subj(recalls, from_mask_rec, to_mask_rec, ...
                                 from_mask_pres, to_mask_pres, ...
				 list_length, bins)
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
  params =  struct('list_length', list_length);

  % conditional_transitions returns the actual and possible transitions for
  % each trial; we simply concatenate them all together and use histc() to
  % gather them up, since we don't care about distinctions between individual
  % trials within a subject
  num_trials = size(recalls, 1);
  for i = 1:num_trials
    % possible_transitions will exclude serial positions masked out by
    % params.to_mask_pres from trial_possibles
    params.to_mask_pres = to_mask_pres(i,:);
    params.from_mask_pres = from_mask_pres(i,:);
    [trial_actuals, trial_possibles] = conditional_transitions(...
					recalls(i, :), ...
    	                                from_mask_rec(i, :), to_mask_rec(i, :), ...
					@lag, @possible_transitions, ...
			                step, params);
    actual_transitions = [actual_transitions, trial_actuals];
    poss_transitions = [poss_transitions, catcell(trial_possibles)];
  end

  % the range of all possible transitions depends only on the list length;
  % we want to gather all the subject's actual and possible transitions in
  % this range
  %all_possible_transitions = [-list_length + 1 : list_length - 1];
  %actuals_counts = collect(actual_transitions, all_possible_transitions);
  %possibles_counts = collect(poss_transitions, all_possible_transitions);

  % the mod: bin actuals and possibles before dividing.
  if length(poss_transitions) == 0
    % unlike collect(), histc() doesn't handle empty vectors very
    % nicely, so return a vector of zeros if there were no possibles
    subj_crp = zeros(size(bins));
  else
    binned_actuals = histc(actual_transitions, bins);
    binned_possibles = histc(poss_transitions, bins);
    subj_crp = binned_actuals ./ binned_possibles;
  end

%endfunction

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
%endfunction
  
% function binned_lags = sumbins(lag_counts, all_lags, bins)
%   binned_lags = zeros(1,length(bins));
%   for i = 1:length(bins) - 1
%     left_edge = bins(i); right_edge = bins(i+1);
%     binned_lags(i) = sum(lag_counts(find(all_lags >= left_edge & ...
% 					 all_lags < right_edge)));
%   end
%   % correct for last bin, emulating histc
%   binned_lags(end) = sum(lag_counts(find(all_lags == bins(end))));
%endfunction

  

