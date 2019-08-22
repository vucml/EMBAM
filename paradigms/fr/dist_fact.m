function dist_facts = dist_fact(recalls_itemnos, pres_itemnos, subjects, ...
                                dist_mat, rec_mask, pres_mask)
%DIST_FACT   Distance-based clustering factor.
%
%  dist_facts = dist_fact(recalls_itemnos, pres_itemnos, subjects, ...
%                         dist_mat, rec_mask, pres_mask)
%
%  INPUTS:
%  recalls_itemnos: a matrix whose elements are identifying item numbers of
%                   recalled items.  The rows of this matrix should represent
%                   recalls made by a single subject on a single trial.
%
%     pres_itemnos: a matrix whose elements are identifying item numbers of
%                   presented items.  The rows of this matrix should represent
%                   presentations to a single subject on a single trial.
%
%         dist_mat: a matrix whose elements are similarity values for pairs of 
%                   presented items.  For two item numbers i1 and i2, 
%                   dist_mat(i1, i2) contains the similarity of those items.
%
%        subjects:  a column vector which indexes the rows of recalls_matrix
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_matrix(find(subjects==S), :)
%
%        rec_mask:  if given, a logical matrix of the same shape as 
%                   recalls_itemnos, which is false at positions (i, j) where
%                   the value at recalls_itemnos(i, j) should be excluded from
%                   the analyses.  If NOT given, a standard clean recalls mask
%                   is used, which excludes repeats, intrusions and empty cells
%
%       pres_mask:  if given, a logical matrix of the same shape as the
%                   item presentation matrix (i.e., num_trials x
%                   list_length).  This mask should be true at
%                   positions (t, sp) where an item in the condition
%                   of interest was presented at serial position sp on
%                   trial t; and false everywhere else.  If NOT given,
%                   a blank mask of this shape is used.
%  OUTPUTS:
%      dist_facts:  a vector of distance clustering factors, one
%                   for each subject.


% sanity checks
if ~exist('recalls_itemnos', 'var')
  error('You must pass a recalls item number matrix.')
elseif ~exist('pres_itemnos', 'var')
  error('You must pass a presented item number matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif size(recalls_itemnos, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
end
if ~exist('rec_mask', 'var')
  % create standard clean recalls mask if none was given
  rec_mask = make_clean_recalls_mask2d(recalls_itemnos);
  rec_mask = remove_intrusions(rec_mask, recalls_itemnos, pres_itemnos);
end
if ~exist('pres_mask', 'var')
  % if no mask given, include all item presentations
  pres_mask = true(size(pres_itemnos));
end
if size(rec_mask) ~= size(recalls_itemnos)
  error('recalls_matrix and mask must have the same shape.')
end


% dist_fact_for_subj will do the work:
dist_facts = apply_by_index(@dist_fact_for_subj, subjects, 1, ...
                            {recalls_itemnos, pres_itemnos, rec_mask, ...
                             pres_mask}, dist_mat);
%endfunction

function subj_dist_fact = dist_fact_for_subj(recalls_itemnos, pres_itemnos, ...
					     rec_mask, pres_mask, dist_mat)
% helper for dist_fact

step = 1;
params = struct('dist_mat',dist_mat);

% call conditional_transitions for each trial
num_trials = size(recalls_itemnos, 1);
trial_facts = [];

for i = 1:num_trials
  % set parameters for this trial to be sent to conditional_transitions
  params.trial_pres_itemnos = pres_itemnos(i,:);
  params.trial_pres_mask = pres_mask(i,:);
  fact = [];
    
  [act_dists, pos_dists] = conditional_transitions(recalls_itemnos(i,:), ...
      rec_mask(i,:), rec_mask(i,:), @similarity, @possible_sim_transitions, ...
      step, params);
            
  % determine actual and possible distances that could be transitioned to
  valid = ~isnan(act_dists);
  act_dists = act_dists(valid);
  pos_dists = pos_dists(valid);

  if isempty(act_dists)
    trial_facts(i) = NaN;
    continue
  end

  % using the possible and actual distances, rank the transition
  for j = 1:length(act_dists)
    possible = pos_dists{j};
    actual = act_dists(j);
    fact(j) = percentile_rank(actual,possible);
  end
  % aggregate the transition ranks to get the trial level rank
  trial_facts(i) = mean(fact);
end

% aggregate the trial factors to get the subject temporal factor
subj_dist_fact = nanmean(trial_facts);
%endfunction

function pd = possible_sim_transitions(item, prior_items, item_sim, params)
  % Helper to return the possible similarity values for each
  % transition.  This is the condition function passed
  % to conditional_transitions.
  unmasked_pres_items = params.trial_pres_itemnos(params.trial_pres_mask);
  poss_items = setdiff(unmasked_pres_items, [item prior_items]);
  pd = params.dist_mat(item, poss_items);
%endfunction


function d = similarity(item1, item2, params)
  % Helper to return the actual similarity of items
  % in a transition.  This is the transition function passed to
  % conditional_transitions.
  d = params.dist_mat(item1, item2);
%endfunction


function mask_out = remove_intrusions(mask_in, rec_itemnos, pres_itemnos)
  % Helper to remove intrusions based on item numbers
  mask_out = mask_in;
  for i=1:size(mask_out,1)
    mask_out(i,:) = mask_out(i,:) & ...
                    ismember(rec_itemnos(i,:), pres_itemnos(i,:));
  end
%endfunction  