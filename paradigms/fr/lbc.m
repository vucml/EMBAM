function lbc_scores = lbc(study_cats, recall_cats, subjects, varargin)
%LBC   List-Based semantic Clustering index.
%
%  Calculate the List-Based semantic Clustering index (LBC-SEM)
%  developed by Stricker et al. (2002). Unlike ARC, LBC takes the list
%  length and number of presented categories into account.
%
%  lbc_scores = lbc(study_cats, recall_cats, subjects, ...)
%
%  INPUTS:
%   study_cats:  [trials X items] matrix of category labels of studied
%                items.
%
%  recall_cats:  [trials X recalls] matrix of category labels of
%                recalled items.
%
%     subjects:  numeric vector of length [trials] indicating the
%                subject identifier for each trial.
%
%  OUTPUTS:
%   lbc_scores:  [subjects X 1] vector of LBC scores, averaged across
%                trials within each subject.
%
%  PARAMS:
%  These options may be specified using parameter, value pairs or by
%  passing a structure.
%   study_mask  - [trials X items] logical array indicating items to
%                 include for calculating LBC. Probably the only use
%                 for this is to exclude whole trials from the
%                 calculation. Default is to include all items.
%   recall_mask - [trials X recalls] logical array indicating recalls
%                 to include. Default is to exclude (remove entirely
%                 from the recall sequence) NaN elements of recall_cats.

% conditions tested in Stricker et al. 2002
% c = {[1 1 2 2] ...
%      [1 1 1 1] ...
%      [1 1 1 1 2 2 2 2] ...
%      [1 1 2 2 3 3 4 4] ...
%      [1 1 1 1 2 2 2 2 3 3 3 3] ...
%      [1 1 1 2 2 2 3 3 3 4 4 4] ...
%      [1 1 2 2 3 4 3 4] ...
%      [1 1 1 1 3 4 3 4] ...
%      [1 1 2 2 1 2 1 2] ...
%      [1 1 1 2 2 2 3 3 3 4 4 4 1 2 3 4]};

% sanity checks
if ~exist('study_cats', 'var')
  error('You must pass a study category matrix.')
elseif ~exist('recall_cats', 'var')
  error('You must pass a recall category matrix.')
elseif size(study_cats, 1) ~= length(subjects)
  error('Study category matrix must have the same number of rows as subjects.')
elseif size(recall_cats, 1) ~= length(subjects)
  error('Recall category matrix must have the same number of rows as subjects.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector')
end

% options
defaults.study_mask = true(size(study_cats));
defaults.recall_mask = ~isnan(recall_cats);
params = propval(varargin, defaults);

lbc_scores = apply_by_index(@lbc_for_subj, subjects, 1, ...
                            {study_cats, recall_cats, ...
                             params.study_mask, params.recall_mask});

function subj_lbc_score = lbc_for_subj(study_cats, recall_cats, ...
                                       study_mask, recall_mask)
  %SUBJ_LBC_SCORE   Mean LBC for one subject.

  n_trials = size(study_cats, 1);
  trial_lbc_scores = [];
  for i=1:n_trials
    trial_study_cats = study_cats(i, study_mask(i,:));
    trial_recall_cats = recall_cats(i, recall_mask(i,:));
    
    if isempty(trial_study_cats) || isempty(trial_recall_cats)
      % LBC is undefined for this trial
      continue
    end
    
    % list length
    nl = length(trial_study_cats);
    
    % check study category
    cats = unique(trial_study_cats);
    n_cats = length(cats);
    m = nl / n_cats;
    if ~all(collect(trial_study_cats, cats) == m)
      error('Assumption of equal number of items per category is violated.')
    end
    
    % number of recalls
    r = length(trial_recall_cats);
    
    % clustering expected by chance
    exp = ((r - 1) * (m - 1)) / (nl - 1);
    
    % actual clustering (number of same-category pairs)
    obs = nnz(diff(trial_recall_cats) == 0);
    
    trial_lbc_scores = [trial_lbc_scores obs - exp];
  end
  
  % sanity check
  if any(isnan(trial_lbc_scores))
    error('LBC score is NaN for at least one trial.')
  end
  
  subj_lbc_score = mean(trial_lbc_scores);
  
