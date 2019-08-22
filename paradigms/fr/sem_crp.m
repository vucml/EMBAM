function [bin_mean,sem_freq_binned_all,slope_subj,slope] = sem_crp(data,sem,numbins,pIncrease,lim_vec,collapse_neg)
%SEM_CRP   Conditional response probability by semantic relatedness.
%
% FUNCTION:
%   [bin_mean,sem_freq_binned,slope_subj,slope,sem_freq,all_vals] = sem_crp(data,sem,numbins,pIncrease,lags)
%
% INPUT ARGS:
%   data = data structure with the following required fields:
%          data.subject
%          data.listLength
%          data.recalls
%          data.rec_itemnos
%          data.pres_itemnos
%   sem = matrix of semantic similarity values
%   numbins = number of bins into which semantic similarity values are
%             seperated and plotted
%   pIncrease = number of times larger the last bin is than the
%               first (in order to spread out the points more evenly on the figure)
%   lags = absolute value of lags you want to analyze, default is
%          all (1:listlength-1)
%
% OUTPUT ARGS:
%   bin_mean = average semantic value in each bin
%   sem_freq_binned = conditional probability of transitioning to
%   each bin
%   slope_subj = slope of each subjs semantic CRP
%   slope = slope of CRP across subjs
%   sem_freq = conditional probability of transitioning to every possible
%               semantic value
%   all_vals = every possible semantic value, same length as
%              sem_freq
%
%
%  Example:
%  [bin_mean,sem_freq_binned,slope_subj,slope,sem_freq,all_vals] = sem_crp(data,WAS,15,4,1:15)
%


% Get info from data structure
iRECNOS=data.rec_itemnos;
iRECALLS=data.recalls;
iPRESNOS=data.pres_itemnos;
nSUBJ = unique(data.subject);
SUBJ = data.subject;
LL = data.listLength;

% set defaults if var is missing
if ~exist('lags','var') || isempty(lags)
  lags = 1:LL-1;
end

if ~exist('sem','var') || isempty(sem)
  fprintf('Semantic Matrix Needed!\n')
end

if ~exist('numbins','var') || isempty(numbins)
  numbins = 15;
end

if ~exist('pIncrease','var') || isempty(pIncrease)
  pIncrease = 4;
end

if ~exist('collapse_neg','var') || isempty(collapse_neg)
  collapse_neg = 0;
end

% find all values in sementic matrix
all_vals = sort(unique(sem(sem<1)));

% find starting range of values to place in each bin
if ~collapse_neg
  binsize = range(all_vals)/numbins;
else
  binsize = max(all_vals)/(numbins-1);
end

% calculate new bin ranges based on pIncrease
if ~collapse_neg
  d = ((pIncrease*binsize)-binsize)/(pIncrease + 1);
  m = (2*d)/(numbins - 1);
  range_vec = m*(1:numbins)+binsize-d-m;
else
  range_vec(1) = range([min(all_vals) 0]);
  d = ((pIncrease*binsize)-binsize)/(pIncrease + 1);
  m = (2*d)/(numbins - 2);
  range_vec(2:numbins) = m*(1:numbins-1)+binsize-d-m;
end

if ~exist('lim_vec','var') || isempty(lim_vec)
  % upper semantic similarity bounds for each bin
  lim_vec = zeros(1,numbins);
  
  % upper indices into all_vals for each bin, bases on lim_vec
  bin_vec = zeros(1,numbins);
  start = min(all_vals);
  for i = 1:numbins;
    lim_vec(i) = start + range_vec(i);
    if i == 1
      new_bin = find(all_vals >= start & all_vals <= lim_vec(i),1,'last');
    else
      new_bin = find(all_vals > start & all_vals <= lim_vec(i),1,'last');
    end
    if isempty(new_bin)
      new_bin = NaN;
    end
    bin_vec(i) = new_bin;
    
    % save the starting limit
    start = lim_vec(i);
  end
  
else
  
  % upper indices into all_vals for each bin, bases on lim_vec
  bin_vec = zeros(1,numbins);
  start = min(all_vals);
  for i = 1:numbins;
    %lim_vec(i) = start + range_vec(i);
    if i == 1
      new_bin = find(all_vals >= start & all_vals <= lim_vec(i),1,'last');
    else
      new_bin = find(all_vals > start & all_vals <= lim_vec(i),1,'last');
    end
    if isempty(new_bin)
      new_bin = NaN;
    end
    bin_vec(i) = new_bin;
    
    % save the starting limit
    start = lim_vec(i);
  end
end

% remove the bins that have no values
lim_vec = lim_vec(~isnan(bin_vec));
bin_vec = bin_vec(~isnan(bin_vec));
new_numbins = length(bin_vec);
if new_numbins ~= numbins
  fprintf('Warning: Some bins did not have any values.\n');
  fprintf('The new numbins is %d.\n',new_numbins);
  numbins = new_numbins;
end

% calculate the mean bin values
bin_mean = zeros(numbins,1);
for b = 1:numbins
  if b > 1
    startInd = bin_vec(b-1);
  else
    startInd = 1;
  end
  
  endInd = bin_vec(b);
  bin_mean(b) = mean(all_vals(startInd:endInd));
end

% will hold probability for each bin for each subj
%sem_freq_all = zeros(length(all_vals),LL-1,'single');

% will hold probability for new, larger bins
sem_freq_binned_all = zeros(length(nSUBJ),numbins,LL-1);

% will hold the slopes for each subj
slope_subj = NaN(length(nSUBJ),LL-1);

slope = NaN(length(nSUBJ),1);

% s goes through each subject
for s = 1:length(nSUBJ)
  
  % create numerator and denominator
  numerator_all = zeros(length(all_vals),LL-1,'int16');
  denominator_all = zeros(length(all_vals),LL-1,'int16');

  % get info for each subject
  iRECNOS_subj = iRECNOS(SUBJ == nSUBJ(s),:);
  iRECALLS_subj = iRECALLS(SUBJ == nSUBJ(s),:);
  iPRESNOS_subj = iPRESNOS(SUBJ == nSUBJ(s),:);
  
   
  nLists = size(iRECNOS_subj,1);
  nPossRec = size(iRECNOS_subj,2);

  % i says what trial this is
  for i = 1:nLists
    
    % j says what output position this is
    for j = 2:nPossRec
      
      % only do the following if prev is an item from the wordpool and is
      % a valid recall
      
      if iRECALLS_subj(i,j-1) > 0 && ... % prev. item is correct recall
            iRECALLS_subj(i,j) > 0 && ... % cur. item is correct recall
            ~isnan(sem(iRECNOS_subj(i,j),iRECNOS_subj(i,j-1))) && ... % has sim
            iRECNOS_subj(i,j-1) ~= iRECNOS_subj(i,j) % is not same word (or intrusion)
        
        % determine the semantic similarity val
        this_sem = sem(iRECNOS_subj(i,j-1),iRECNOS_subj(i,j));
	
	% determine the lag for the transition
	this_lag = abs(iRECALLS_subj(i,j-1) - iRECALLS_subj(i,j));
        
        % increment the numerator for that sim val and lag
        binIndN = this_sem == all_vals;
        numerator_all(binIndN,this_lag) = numerator_all(binIndN,this_lag) + 1;
        
        % find all possible transitions
        pos_trans = iPRESNOS_subj(i,~ismember(iPRESNOS_subj(i,:), ...
                                              iRECNOS_subj(i,1:j-1)));
	
	% find semantic values for all possible transitions
        pos_trans_val = sem(iRECNOS_subj(i,j-1),pos_trans);
	
	% NaNs?
	NaNLoc = find(isnan(pos_trans_val));
        pos_trans_val(NaNLoc) = [];
        pos_trans(NaNLoc) = [];
        
        % find semantic bins for all possible transitions
        binIndD = zeros(1,length(pos_trans));
        for k = 1:length(pos_trans)
          binIndD(k) = find(pos_trans_val(k) == all_vals);
        end
        	
	% determine the abs(lag) for all possible transitions
	pos_lags = zeros(1,length(pos_trans));
	for k = 1:length(pos_trans)
	  pos_lags(k) = abs(iRECALLS_subj(i,j-1) - find(iPRESNOS_subj(i,:) == pos_trans(k)));
	end
	
	% increment the denominator
	for k = 1:length(binIndD)
	  denominator_all(binIndD(k),pos_lags(k)) = denominator_all(binIndD(k),pos_lags(k)) + 1;
	end
	
      end
    end
  end
        
    
  % calculate probability for each bin

  denominator_all = single(denominator_all);
  denominator_all(denominator_all == 0) = NaN;
  numerator_all = single(numerator_all);
  sem_freq_all = numerator_all./denominator_all;
  if length(find(denominator_all == 0)) > 0
    keyboard
  end

  
 for lag = 1:LL-1
    % find prob for new bins for each subj
    start = 1;
    for i = 1:numbins
      
      % number of all possible transitions within new, larger bin
      if nansum(denominator_all(start:bin_vec(i),lag)) == 0;
	sumBin_denom = NaN;
      else
	sumBin_denom = nansum(denominator_all(start:bin_vec(i),lag));
      end
      
      % find weight each semantic value gives to the bin
      weights = denominator_all(start:bin_vec(i),lag)/sumBin_denom;
      
      if isempty(weights(weights > 0))
	sem_freq_binned_all(s,i,lag) = NaN;
      else
	%find probability for each new bin
	sem_freq_binned_all(s,i,lag) = nansum(weights.*sem_freq_all(start:bin_vec(i),lag));
      end
      
      % update starting index
      start = bin_vec(i) + 1;
      
    end
    
    % regression to find slope
    if length(find(~isnan(sem_freq_binned_all(s,:,lag)))) > 2
      %B = robustfit(all_vals,sem_freq_all(:,lag));
%      B = regress(sem_freq_all(:,lag),[all_vals ones(length(all_vals),1)]);
      B = robustfit(bin_mean,sem_freq_binned_all(s,:,lag));
      slope_subj(s,lag) = B(2);
    else
      slope_subj(s,lag) = NaN;
    end
    if length(find(~isnan(nanmean(sem_freq_binned_all(s,:,:),3)))) > 2
        B = robustfit(bin_mean,nanmean(sem_freq_binned_all(s,:,:),3));
        slope(s) = B(2);
    else
        slope(s) = NaN;
    end
    
  end

  fprintf('%f\t',s/length(nSUBJ)*100)
end % subj
fprintf('\n')

% create denom and num for all subjs combined
%denominator = double(squeeze(nansum(denominator_all,1)));
%numerator = double(squeeze(nansum(numerator_all,1)));

% change zeros to NaNs
%denominator(denominator == 0) = NaN;

% finally, find overall semCRP by meaning across lag
%sem_freq = nanmean(numerator./denominator,2);

% regression to find slope (unbinned)
%B = regress(sem_freq,[all_vals ones(length(all_vals),1)]);
%slope = B(1);

