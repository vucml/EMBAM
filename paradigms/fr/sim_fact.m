function sim_facts = sim_fact(recalls_itemnos, pres_itemnos, ...
                              sim_mat, subjects, rec_mask, pres_mask)
%SIM_FACT   Similarity-based clustering factor.
%
%  Determine whether recall is organized by some measure that is
%  continuous and varies between item pairs, for example an
%  estimate of semantic similarity.
%
%  sim_facts = sim_fact(recalls_itemnos, pres_itemnos, sim_mat, 
%                       subjects, rec_mask, pres_mask)
%
%  INPUTS
%  recalls_itemnos - [lists x recalls] numeric array
%      Item numbers of recalled items.
%
%  pres_itemnos - [lists x serial position] numeric array
%      Item numbers of presented items.
%
%  sim_mat - [items x items] numeric array
%      sim_mat(i,j) gives the similarity between item numbers i and j.
%
%  subjects - [lists x 1] numeric array
%      Labels indicating which lists belong to each subject. Each
%      subject label can be any numeric value, as long as it is
%      unique to that subject.
%
%  rec_mask - [lists x recalls] numeric array
%      Logical matrix that is true only for recalls that should be
%      included in the analysis. If not specified, will exclude
%      repeats and intrusions.
%
%  pres_mask - [lists x serial position] numeric array
%      Logical matrix that is true only for items that should be
%      included in the analysis. If not specified, all items will
%      be included.
%
%  OUTPUTS
%  sim_facts - [subjects x 1] numeric array
%      Similarity clustering factor. The mean percentile rank of
%      similarity of actual transitions, relative to possible
%      transitions. Values above 0.5 indicate clustering by
%      similarity. Output is in order of sorted subject number.
%
%  See also: sim_fact_cat, temp_fact.

% sanity checks
if size(recalls_itemnos, 1) ~= length(subjects)
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
if any(size(rec_mask) ~= size(recalls_itemnos))
  error('recalls_matrix and mask must have the same shape.')
end

% run analysis on each subject
sim_facts = apply_by_index(@sim_fact_for_subj, subjects, 1, ...
                            {recalls_itemnos, pres_itemnos, rec_mask, ...
                             pres_mask}, sim_mat);
%endfunction


function subj_sim_fact = sim_fact_for_subj(recalls_itemnos, pres_itemnos, ...
                                           rec_mask, pres_mask, sim_mat)

% calculate similarity percentile for each trial
step = 1;
params = struct('sim_mat', sim_mat);
num_trials = size(recalls_itemnos, 1);
trial_facts = NaN(num_trials, 1);
for i = 1:num_trials
  % set parameters for this trial to be sent to conditional_transitions
  params.trial_pres_itemnos = pres_itemnos(i,:);
  params.trial_pres_mask = pres_mask(i,:);
  
  % actual and possible similarity for included transitions only
  [act_sims, pos_sims] = conditional_transitions(recalls_itemnos(i,:), ...
      rec_mask(i,:), rec_mask(i,:), @similarity, @possible_sim_transitions, ...
      step, params);
  
  % if no included transitions, percentile for this trial is undefined
  valid = ~isnan(act_sims);
  if nnz(valid) == 0
    continue
  end
  
  % rank the actual transitions relative to possible transitions
  act_sims = act_sims(valid);
  pos_sims = pos_sims(valid);
  fact = NaN(length(act_sims), 1);
  for j = 1:length(act_sims)
    fact(j) = percentile_rank(act_sims(j), pos_sims{j});
  end

  % average percentile for this trial
  trial_facts(i) = mean(fact);
end

% average percentile for this subject
subj_sim_fact = nanmean(trial_facts);
%endfunction

function ps = possible_sim_transitions(item, prior_items, item_sim, params)
  % Helper to return the possible similarity values for each
  % transition.  This is the condition function passed
  % to conditional_transitions.
  unmasked_pres_items = params.trial_pres_itemnos(params.trial_pres_mask);
  poss_items = setdiff(unmasked_pres_items, [item prior_items]);
  ps = params.sim_mat(item, poss_items);
%endfunction


function s = similarity(item1, item2, params)
  % Helper to return the actual similarity of items
  % in a transition.  This is the transition function passed to
  % conditional_transitions.
  s = params.sim_mat(item1, item2);
%endfunction


function mask_out = remove_intrusions(mask_in, rec_itemnos, pres_itemnos)
  % Helper to remove intrusions based on item numbers
  mask_out = mask_in;
  for i=1:size(mask_out,1)
    mask_out(i,:) = mask_out(i,:) & ...
                    ismember(rec_itemnos(i,:), pres_itemnos(i,:));
  end
%endfunction
