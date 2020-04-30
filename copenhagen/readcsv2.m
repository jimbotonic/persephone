% readcsv2.m
% created 4-22-2020 by Jeffrey W. Herrmann

% This reads data but then creates matrix only 
% for certain range of times.

% skip first row

% nrows = 5474289; % number of rows to read
%  
% M = csvread('bt_symmetric.csv',1,0,[1 0 nrows 2]);
%   
% 
% maxtime = max(M(:,1));  % max time (in seconds)
% maxid = max(M(:,2))+1;  % max variable number
% maxvarn = maxid; 
% 
% M = M(M(:,3)>=0,:);  % keep rows with non-negative entry 3 

% RANGES FOR WEEKDAY START HERE
% nranges = 20;
% timerange = zeros(nranges,2);   % beginning and end of each time range in hours
% timerange(1,1:2) = [24+7,24+17]; % weekdays go from 700 to 1700
% for rj=2:5
%     timerange(rj,1:2) = timerange(rj-1,1:2)+24;
% end
% timerange(5,2) = timerange(5,1) + 5;  % Friday ends at 1200
% 
% % add three more weeks
%  for wn=1:3
%      timerange(wn*5+1:wn*5+5,1:2) = timerange(wn*5-4:wn*5,1:2) + 168;  % add one week (168 hours)
%  end 

% RANGES FOR WEEKDAY NIGHTS START HERE
% nranges = 20;
% timerange = zeros(nranges,2);   % beginning and end of each time range in hours
% timerange(1,1:2) = [24,24+6]; % weekdays go from 0000 to 0600
% for rj=2:5
%     timerange(rj,1:2) = timerange(rj-1,1:2)+24;
% end
% 
% % add three more weeks of five days
% ndays = 5; % days per week
%  for wn=1:3
%      timerange(wn*ndays+1:wn*ndays+ndays,1:2) = timerange(wn*ndays-ndays+1:wn*ndays,1:2) + 168;  % add one week (168 hours)
%  end 

% RANGES FOR WEEKEND NIGHTS START HERE
nranges = 8;
timerange = zeros(nranges,2);   % beginning and end of each time range in hours
timerange(1,1:2) = [5*24+16,6*24]; % Friday 4:00 to midnight
timerange(2,1:2) = [6*24+16,7*24]; % Saturday 4:00 to midnight

% add three more weeks (2 days per week)
 for wn=1:3
     timerange(wn*2+1:wn*2+2,1:2) = timerange(wn*2-1:wn*2,1:2) + 168;  % add one week (168 hours)
 end 

 
 % create sparse matrix
 S = sparse(maxvarn,maxvarn);

 for rj=1:nranges % loop over ranges
     tstart = timerange(rj,1)*3600;  % start in seconds
     tend =   timerange(rj,2)*3600;  % end in seconds
     datalist = find((M(:,1)>=tstart) .* (M(:,1)<=tend)); % rows in M that are in this range

     mv2 = max(M(datalist,2))+1;  % max variable in column 2
     mv3 = max(M(datalist,3))+1;  % max variable in column 3
     
     S(1:mv2,1:mv3) = S(1:mv2,1:mv3) + sparse(M(datalist,2)+1,M(datalist,3)+1,ones(length(datalist),1));  % add matrix for this range

     % create sparse adjacency matrix; add one to every entry
     % note that this accumulates, so it also counts the number of times
     % that the pair appears.

 end % loop over ranges
 

