function temp_facts = general_temp_fact(recalls, pres_catlabels, ...
					subjects, cat_type, from_mask_rec, ...
					to_mask_rec, from_mask_pres, ...
					to_mask_pres, signed)
%GENERAL_TEMP_FACT   Temporal clustering conditional on category structure.
%
%  temp_facts = general_temp_fact(recalls, pres_catlabels, ...
%                                 subjects, cat_type, from_mask_rec, ...
%                                 to_mask_rec, from_mask_pres, ...
%                                 to_mask_pres, signed);
%
%  INPUTS:
%         recalls:  a matrix whose elements are serial positions of
%                   recalled items.  The rows of this matrix should
%                   represent recalls made by a single subject on a
%                   single trial.
%
%  pres_catlabels:  a matrix whose elements are the category label of
%                   the item at each serial position.  The rows of
%                   this matrix represent the category presentation
%                   schedule for a single subject on a single trial.
%                   The number of columns is used to set list_length.
%
%        subjects:  a column vector which indexes the rows of
%                   recalls with a subject number (or other
%                   identifier).  That is, the recall trials of
%                   subject S should be located in
%                   recalls(find(subjects==S), :)
%
%        cat_type:  a flag.  If set to 1, calculates the temporal
%                   factor for same category transitions.  If set to
%                   2, calculates the temporal factor for between
%                   category transitions.
%
%   from_mask_rec:  a logical matrix of the same shape as
%                   recalls, which is false at positions (i, j)
%                   where the value at recalls(i, j) should be
%                   excluded from the calculation of the probability
%                   of recall. from_mask excludes transitions
%                   that initiate with item j.
%
%     to_mask_rec:  a logical matrix, as described above.  to_mask
%                   excludes transitions that terminate on item j+1.
% 
%  from_mask_pres:  a [list X item] logical matrix, which is false at
%                   positions (i,j) where transitions from that item
%                   should be excluded.
%
%    to_mask_pres:  a [list X item] logical matrix, which is false at
%                   positions (i,j) where transitions to that item
%                   should be excluded.
%
%          signed:  a logical value that indicates whether to
%                   preserve the sign of each lag in the set of
%                   percentile ranks.  If true, a lag of -2 will
%                   have temporal factor:
%                     (-1 * percentile_rank(-2, set_of_possible_transitions))
%                   while a lag of +2 will have temporal factor:
%                     percentile_rank(-2, set_of_possible_transitions)
%
%  OUTPUTS:
%      temp_facts:  a vector of temporal clustering factors, one
%                   for each subject.
%
%  NOTES:
%   This version aggregates by trial, then takes the set of
%   trial means.  Prior version aggregated all transitions per
%   subject and then took the mean.
%
%  See also temp_fact.

% sanity checks
if ~exist('recalls', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('pres_catlabels', 'var')
  error('You must pass in pres_catlabels.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('cat_type', 'var')
  error('You must specify type of category transitions.')
elseif ~exist('from_mask_rec', 'var')
  error('You must provide from-item recalls mask.')
elseif ~exist('to_mask_rec', 'var')
  error('You must provide to-item recalls mask.')
elseif ~exist('from_mask_pres', 'var')
  error('You must provide from_item presentation mask.')
elseif ~exist('to_mask_pres', 'var')
  error('You must provide to_item presentation mask.')
elseif size(recalls, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
elseif size(pres_catlabels, 1) ~= length(subjects)
  error('pres_catlabels must have the same number of rows as subjects.')
elseif size(from_mask_rec, 1) ~= length(subjects)
  error('from_mask_rec must have the same number of rows as subjects.')
elseif size(to_mask_rec, 1) ~= length(subjects)
  error('to_mask_rec must have the same number of rows as subjects.')
elseif size(from_mask_pres, 1) ~= length(subjects)
  error('from_mask_pres must have the same number of rows as subjects.')
elseif size(to_mask_pres, 1) ~= length(subjects)
  error('to_mask_pres must have the same number of rows as subjects.')
elseif ~exist('signed', 'var')
  error('You must indicate whether you want to preserve lag signedness')
elseif ~islogical(signed) || length(signed) > 1
  error('signed must be a logical scalar')
end

% from_mask_rec and to_mask_rec are used to determine which
% transitions are examined.  Here, we transfer additional
% conditions from the presentation masks onto these recall masks.
to_mask_rec = update_recalls_mask(recalls, to_mask_rec, to_mask_pres);
from_mask_rec = update_recalls_mask(recalls, from_mask_rec, from_mask_pres);

% temp_fact_for_subj will do the work:
temp_facts = apply_by_index(@general_temp_fact_for_subj, ...
			    subjects, ...
			    1, ...
			    {recalls, pres_catlabels, ... 
		             from_mask_rec, to_mask_rec, ...
		             from_mask_pres, to_mask_pres}, ...
			    cat_type, signed);
%endfunction

function subj_temp_fact = general_temp_fact_for_subj(recalls, ...
				pres_catlabels, from_mask_rec, to_mask_rec, ...
		                from_mask_pres, to_mask_pres, cat_type, signed)
  % helper for general_temp_fact:
  % calculates the percentile score for each recall transition,
  % aggregates these by trial into a mean score.  Aggregates across
  % trials within subject to a mean score, which is returned.
  %
  %

  step = 1;
  params = struct('list_length', size(pres_catlabels,2), ...
		  'cat_type', cat_type);

  % call conditional_transitions for each trial
  num_trials = size(recalls, 1);
  % trial_facts = [];
  all_facts = [];
  for i = 1:num_trials
    fact = [];
    params.pres_catlabels = pres_catlabels(i,:);
    params.from_mask_pres = from_mask_pres(i,:);
    params.to_mask_pres = to_mask_pres(i,:);
    [act_trans, poss_trans] = conditional_transitions(...
	recalls(i, :), from_mask_rec(i, :), to_mask_rec(i,:), ...
	@cat_lag, @cat_transitions, ...
	step, params);
    % for each valid transition call percentile_rank
    valid = ~isnan(act_trans);
    act_trans = act_trans(valid);
    poss_trans = poss_trans(valid);
    for j = 1:length(act_trans)
      possible = poss_trans{j};
      possible = -1 * abs(possible);
      actual = -1 * abs(act_trans(j));
      if signed
	fact(j) = sign(act_trans(j)) * percentile_rank(actual, possible);
      else
	fact(j) = percentile_rank(actual,possible);
      end
    end
    % aggregate the transition ranks to get the trial level rank
    % trial_facts(i) = mean(fact);
    % aggregate the transition ranks into an across-trial vector 
    all_facts = [all_facts fact];
  end
  % aggregate the trial factors to get the subject temporal factor
  % subj_temp_fact = mean(trial_facts);
  % aggregate all factors to get the subject temporal factor
  % nanmean: can have an undefined percentile rank if all items of
  % a particular category are recalled
  subj_temp_fact = nanmean(all_facts);

