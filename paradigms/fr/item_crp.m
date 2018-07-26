function [actual, possible] = item_crp(recalls, pres_itemnos, index, n_wordpool, varargin)
%ITEM_CRP   Conditional response probability by item pair.
%
%  [actual, possible] = item_crp(recalls, pres_itemnos, index, n_wordpool, ...)
%
%  INPUTS:
%  recalls - [lists x items] numeric array
%      Serial positions of recalled items.
%
%  pres_itemnos - [lists x items] numeric array
%      Position of each presented item in the wordpool.
%
%  index - [lists] numeric array
%      Unique label for each group of lists that should be separately
%      analyzed.
%
%  n_wordpool - scalar
%      Number of items in the wordpool.
%
%  OUTPUTS:
%  actual - [items x items x subjects] numeric array
%      Count of the number of times a transition between each pair
%      of items occurred, for a given group of lists/subject.
%
%  possible - [items x items x subjects] numeric array
%      Count of the number of times a transition could have
%      occurred, given recall of the first item when the second
%      item was still available for recall.
%
%  OPTIONS:
%  from_mask_pres - [lists x items] logical array - true(size(pres_itemnos))
%      Include transitions from these items.
%
%  to_mask_pres - [lists x items] logical array - true(size(pres_itemnos))
%      Include transitions to these items.
%
%  from_mask_rec - [lists x items] logical array - 
%    make_clean_recalls_mask2d(recalls)
%      Include transitions from these recalls.
%
%  to_mask_rec - [lists x items] logical array - 
%    make_clean_recalls_mask2d(recalls)
%      Include transitions to these recalls.

% options
def.from_mask_pres = [];
def.to_mask_pres = [];
def.from_mask_rec = [];
def.to_mask_rec = [];
opt = propval(varargin, def);

[n_trials, n_items] = size(pres_itemnos);

% default presentation masks
if isempty(opt.from_mask_pres)
  opt.from_mask_pres = true(n_trials, n_items);
end
if isempty(opt.to_mask_pres)
  opt.to_mask_pres = true(n_trials, n_items);
end

% default recall masks
if isempty(opt.from_mask_rec)
  opt.from_mask_rec = make_clean_recalls_mask2d(recalls);
end
if isempty(opt.to_mask_rec)
  opt.to_mask_rec = make_clean_recalls_mask2d(recalls);
end

crps = apply_by_index(@subj_item_crp, index, 3, ...
                      {recalls, pres_itemnos, ...
                       opt.from_mask_pres, opt.to_mask_pres, ...
                       opt.from_mask_rec, opt.to_mask_rec}, ...
                      n_wordpool);

% unpack
actual = crps(:,:,:,1);
possible = crps(:,:,:,2);


function subj_crp = subj_item_crp(recalls, pres_itemnos, ...
             from_mask_pres, to_mask_pres, from_mask_rec, to_mask_rec, ...
             n_wordpool)

  [n_trials, n_items] = size(pres_itemnos);
  n_recalls = size(recalls, 2);
  
  actual = zeros(n_wordpool, n_wordpool, 'uint8');
  possible = zeros(n_wordpool, n_wordpool, 'uint8');
  n_trans = 0;
  for i = 1:n_trials
    priors = [];
    for j = 1:(n_recalls - 1)
      from_pt = recalls(i,j);
      to_pt = recalls(i,j+1);
      
      % if this isn't a valid transition, nothing to do
      if ~from_mask_rec(i,j) || ~to_mask_rec(i,j+1)
        priors = [priors from_pt];
        continue
      end
      n_trans = n_trans + 1;
      
      % get valid lags for recall j+1
      act_lag = to_pt - from_pt;
      params.from_mask_pres = from_mask_pres(i,:);
      params.to_mask_pres = to_mask_pres(i,:);
      poss_lag = possible_transitions(from_pt, priors, [], params);
      
      % add this recall to the record of prior recalls
      priors = [priors from_pt];
      
      poss_pt = from_pt + poss_lag;
      from_item = pres_itemnos(i,from_pt);
      to_item = pres_itemnos(i,to_pt);
      poss_item = pres_itemnos(i,poss_pt);
      
      actual(from_item, to_item) = actual(from_item, to_item) + 1;
      possible(from_item, poss_item) = possible(from_item, poss_item) + 1;
    end
  end
  
  if any(possible(:) == 255)
    error('Hit limit of uint8.')
  end
  
  % pack outputs
  subj_crp = cat(4, actual, possible);
  