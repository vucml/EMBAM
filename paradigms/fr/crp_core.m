function [lag_crps] = crp_core(recalls_matrix, subjects, list_length, ...
                                        from_mask_rec, to_mask_rec, ...
                               from_mask_pres, to_mask_pres)
%CRP_CORE   Conditional response probability as a function of lag (lag-CRP).
%
%  Unlike CRP, this function does not have any error-checking. It
%  assumes all of the inputs listed below are given and formatted
%  correctly.
%
%  lag_crps = crp_core(recalls_matrix, subjects, list_length, ...
%                      from_mask_rec, to_mask_rec, ...
%                      from_mask_pres, to_mask_pres)
%
%  See crp for details.

% crp_for_subj is going to do all the real work:
lag_crps = apply_by_index(@crp_for_subj, subjects, 1, ...
                          {recalls_matrix, from_mask_rec, to_mask_rec, ...
                           from_mask_pres, to_mask_pres}, ...
                          list_length);


function subj_crp = crp_for_subj(recalls, from_mask_rec, to_mask_rec, ...
                                 from_mask_pres, to_mask_pres, list_length)
  % helper for crp: calculates the lag-CRP for each possible transition
  % for one subject's recall trials; returns a row vector of CRPs,
  % indexed by (-list_length + index)
  
  % the golden goose: these will store actual and possible transitions
  % for this subject's trials
  actual_transitions = [];
  poss_transitions = [];

  % arguments for conditional_transitions:
  step = 1;   % for now, step of 1 is hard-coded; this may change   
  params = struct('list_length', list_length);

  % conditional_transitions returns the actual and possible transitions
  % for each trial; we simply concatenate them all together and use
  % collect() to gather them up, since we don't care about
  % distinctions between individual trials within a subject
  num_trials = size(recalls, 1);
  for i = 1:num_trials
    % possible_transitions will exclude serial positions masked out by
    % params.to_mask_pres from trial_possibles
    params.to_mask_pres = to_mask_pres(i,:);
    params.from_mask_pres = from_mask_pres(i,:);
    [trial_actuals, trial_possibles] = ...
        conditional_transitions(recalls(i,:), from_mask_rec(i,:), ...
                                to_mask_rec(i,:), ...
        			@lag, @possible_transitions, ...
                                step, params);
    actual_transitions = [actual_transitions, trial_actuals];
    poss_transitions = [poss_transitions, catcell(trial_possibles)];
  end

  % the range of all possible transitions depends only on the list
  % length; we want to gather all the subject's actual and possible
  % transitions in this range
  all_possible_transitions = [-list_length + 1 : list_length - 1];
  actuals_counts = collect(actual_transitions, all_possible_transitions);
  possibles_counts = collect(poss_transitions, all_possible_transitions);
  
  % and that's it!
  subj_crp = actuals_counts ./ possibles_counts;

%endfunction

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
%endfunction