function sem_crps = semantic_crp(recalls_itemnos, pres_itemnos, subjects, ...
				 sem_sims, bins, from_mask, to_mask, pres_mask)
%SEMANTIC_CRP   Conditional response probability by semantic relatedness.
%
% Computes conditional response probabilities as a
% function of semantic similarity from a matrix of recalled serial positions.
%
%  sem_crps = semantic_crp(recalls_itemnos, pres_itemnos, subjects,...
%                          sem_sims, bins, mask)
%
%  INPUTS:
% recalls_itemnos:  a matrix whose elements are identifying item
%                   numbers of recalled items.  The rows of this
%                   matrix should represent recalls made by a
%                   single subject on a single trial.
% 
%    pres_itemnos:  a matrix whose elements are identifying item
%                   numbers of presented items.  The rows of this
%                   matrix should represent presentations to a
%                   single subject in a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls_itemnos
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_itemnos(find(subjects==S), :)
% 
%        sem_sims:  a matrix whose elements are semantic similarity
%                   values for pairs of presented items.  For two
%                   item numbers i1 and i2, sem_sims(i1, i2) should
%                   contain the semantic similarity of those items.
%
%            bins:  bin edges, as they would be passed to histc(), for
%                   bins in which to place (possible and actual)
%                   semantic similarity values
%
%       from_mask:  if given, a logical matrix of the same shape as 
%                   recalls_itemnos, which is false at positions (i, j) where
%                   the value at recalls_itemnos(i, j) should be excluded from
%                   the calculation of the probability of recall.  If NOT
%                   given, a standard clean recalls mask is used, which 
%                   excludes repeats, intrusions and empty cells
%
%         to_mask:  if given, a logical matrix of the same shape as
%                   rec_itemnos, which is false at positions (i,
%                   j) where the transition from
%                   recalls_itemnos(i, j-1) TO recalls_itemnos(i, j)
%                   should be excluded from the calculation of the
%                   semantic CRP.
% 
%                   If neither from_mask nor to_mask is given, a
%                   standard clean recalls mask is used, which
%                   excludes repeats, intrusions and empty cells.  If
%                   from_mask is given but to_mask is not, from_mask
%                   will be used for both masks (i.e., transitions
%                   both to and from masked out elements will be
%                   excluded).
% 
%       pres_mask:  if given, a logical matrix of the same shape as
%                   the pres_itemnos matrix, i.e., 
%                   (num_trials x list_length).  This mask should
%                   be true at positions (t, sp) where an item in
%                   the condition of interest was presented at
%                   serial position sp on trial t; and false
%                   everywhere else.  If NOT given, a blank mask of
%                   this shape is used.
%
%  OUTPUTS:
%        sem_crps:  a matrix of semantic CRP values.  Each row contains
%                   values for one subject, with the value in each
%                   column j corresponding to the subject's probability of
%                   transitioning to items with semantic similarity values
%                   such that bins(j) <= similarity value < bins(j+1).
%
%

% sanity is always important
if ~exist('recalls_itemnos', 'var')
  error('You must pass a recalls-by-item-numbers matrix.')
elseif ~exist('pres_itemnos', 'var')
  error('You must pass a presentations-by-item-numbers matrix.')
elseif ~exist('sem_sims', 'var')
  error('You must pass a semantic-similarities matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('bins', 'var')
  error('You must pass bins')
elseif size(recalls_itemnos, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
elseif ~exist('from_mask', 'var')
  % create standard clean recalls from_mask if none was given
  from_mask = make_clean_recalls_mask2d(recalls_itemnos);
end

if ~exist('to_mask', 'var')
  % assume to_mask should be the same as from_mask, i.e.,
  % transitions to and from the same points should be excluded
  to_mask = from_mask;
end

if (size(from_mask) ~= size(recalls_itemnos)) | ...
   (size(to_mask) ~= size(recalls_itemnos))
  error('recalls_itemnos and masks must have the same shape.')
end

if ~exist('pres_mask', 'var')
  % assume pres_mask should allow all presented items
  pres_mask = make_blank_mask(pres_itemnos);
end

if max(unique(pres_itemnos)) > size(sem_sims, 1) || ...
       size(sem_sims, 1) ~= size(sem_sims, 2)
  % sem_sims must be a square matrix with as many rows/cols as item numbers
  % to guarantee sem_sims(i, j) contains semantic similarity of all <i,j>
  error(['sem_sims must have the same number of rows and columns as ' ...
	 'there are presented item numbers.'])
end

% sem_crp_for_subj is going to do all the real work:
sem_crps = apply_by_index(@sem_crp_for_subj, ...
			  subjects, ...
			  1, ...
			  {recalls_itemnos, pres_itemnos, ...
		           from_mask, to_mask, pres_mask}, ...
			  bins, sem_sims);
%endfunction

function subj_crp = sem_crp_for_subj(recalls_itemnos, pres_itemnos, ...
				     from_mask, to_mask, pres_mask, ...
				     bins, sem_sims)
  % Helper for semantic_crp:
  % calculates the sem-CRP for each possible transition for one subject's
  % recall trial
  
  % arguments for conditional_transitions:
  step = 1;   % for now, step of 1 is hard-coded; this may change   
  params =  struct('sem_sims', sem_sims);

  % state variables:
  num_trials = size(recalls_itemnos, 1);
  num_bins = size(bins, 2);
  actuals_counts = zeros(1, num_bins);  % number of actual transits in each bin
  possibles_counts = zeros(1, num_bins); % number of possible
                                         % transits in each bin

  for i = 1:num_trials
    params.trial_pres_itemnos = pres_itemnos(i, :);
    params.trial_pres_mask = pres_mask(i, :);
    [trial_actuals, trial_possibles] = conditional_transitions(...
					recalls_itemnos(i, :), ...
	                                from_mask(i, :), ...
	                                to_mask(i, :), ...
					@semantic_similarity, ...
                                        @possible_sem_transitions, ...
			                step, params);

    % NOTES: the right and proper way to do this, apparently, is to
    % "collapse" after every transition, and add 1 to each actual
    % and possible bin where an actual or possible value occurred
    % (regardless of how many occurred in that bin on that transition).
    % So if a given transition had an actual relatedness of 0.5 and possible
    % relatedness values [0.2 0.3 0.5 0.8], and the bins were 
    % [0->0.4, 0.4->0.7, 0.7->0.9, 0.9->1.0], then we would add
    % [0 1 0 0] to the actuals count and [1 1 1 0] to the possibles count

    collapsed_actuals = histc(trial_actuals, bins);
    collapsed_possibles = zeros(1, num_bins);
    for j = 1:size(trial_actuals, 2)
      % collapse the possible transitions from each actual
      % transition; this ensures we don't count more than one
      % possible transition in each bin for each actual
      % transition in a trial
      if isempty(trial_possibles{j})
	continue
      end
      trial_possibles_counts = histc(trial_possibles{j}, bins);

      collapsed_possibles = collapsed_possibles + ...
	                    logical(trial_possibles_counts);
    end

    actuals_counts = actuals_counts + collapsed_actuals;
    possibles_counts = possibles_counts + collapsed_possibles;
  end

  % and that's it!
  subj_crp = (actuals_counts ./ possibles_counts);

%endfunction

function ps = possible_sem_transitions(item, prior_items, item_sem_sim, params)
  % Helper to return the possible semantic similarity values for each
  % transition.  The possible items for a given transition are those
  % which were presented but have not been recalled, and which are not
  % masked out by params.trial_pres_mask.  This is the condition
  % function passed to conditional_transitions.
  pres_itemnos = params.trial_pres_itemnos;
  unmasked_pres_itemnos = pres_itemnos(params.trial_pres_mask);
  poss_items = setdiff(unmasked_pres_itemnos, prior_items);
  ps = params.sem_sims(item, poss_items);
%endfunction
  

function s = semantic_similarity(item1, item2, params)
  % Helper to return the actual semantic similarity of items
  % in a transition.  This is the transition function passed to
  % conditional_transitions.
  s = params.sem_sims(item1, item2);
%endfunction
  
