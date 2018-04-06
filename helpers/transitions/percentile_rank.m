function rank = percentile_rank(actual,possible)
%  PERCENTILE_RANK  Computes the rank of actual from among the set
%  of possible values.  
%
%  rank = percentile_rank(actual,possible)
%
%  INPUTS:
%      actual:  A scalar specifying the *actual* value.
%
%    possible:  A vector of all *possible* values.
%
%  OUTPUTS:
%        rank:  A scalar, the percentile ranking of the actual value from
%               among the possible values.
%

% sanity check: is the actual transition one of the possible transitions?
if ~ismember(actual,possible)
    error('The actual value must be one of the possible values.')
end

poss_ranks = tiedrank(possible);
act_rank = poss_ranks(actual==possible);
length_p_ranks = sum(~isnan(poss_ranks)) -1;

rank = unique((act_rank - 1)  / length_p_ranks);