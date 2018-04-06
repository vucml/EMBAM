function temp_facts = temp_fact(recalls, subjects, list_length,...
    from_mask_rec, from_mask_pres)
%TEMP_FACT   Lag-based temporal clustering factor.
%
%  temp_facts = temp_fact(recalls, subjects, list_length, ...
%                         rec_mask, pres_mask)
%
%  INPUTS:
%         recalls:  a matrix whose elements are serial positions of recalled
%                   items.  The rows of this matrix should represent recalls
%                   made by a single subject on a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls(find(subjects==S), :)
%
%     list_length:  a scalar indicating the number of serial positions in the
%                   presented lists.  serial positions are assumed to run 
%                   from 1:list_length.
%     
%  from_mask_pres:  ?
%
%  OUTPUTS:
%      temp_facts:  a vector of temporal clustering factors, one
%                   for each subject.
%
%  NOTES:
%      This version aggregates by trial, then takes the set of
%      trial means.  Prior version aggregated all transitions per
%      subject and then took the mean.

% sanity checks
if ~exist('recalls', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
end

% default values
pres_catlabels = ones(length(subjects),list_length);
cat_type = 1;

% create default masks
if ~exist('from_mask_rec', 'var')
  from_mask_rec = make_clean_recalls_mask2d(recalls);
end
to_mask_rec = from_mask_rec;

if ~exist('from_mask_pres', 'var')
  from_mask_pres = ones(length(subjects),list_length);
end
%to_mask_pres = ones(length(subjects),list_length);
to_mask_pres = from_mask_pres;


% general_temp_fact will do the work:
temp_facts = general_temp_fact(recalls, pres_catlabels, ...
			       subjects, cat_type, ...
			       from_mask_rec, to_mask_rec, ...
			       from_mask_pres, to_mask_pres, ...
                               false);
