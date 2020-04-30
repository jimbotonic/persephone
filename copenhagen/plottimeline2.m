% plottimeline2.m

% Created 4-21-2020 by Jeffrey W. Herrmann.
% NOTE: this re-uses and adapts some code written by 
% Erica L. Gralla (GWU).

% This includes clustering from clustering2.m
% and selecting strong, maximal clusters from pickclusters.m

% This should do the following for the clusters with 1's in picklist:
% output timeline of cluster showing number of persons paired in that
% time bin.
% This is for strong and maximal clusters.

% This creates one timeline for each cluster and takes maximum over the
% hour as the value for the hour.

% FROM clustering2.m

% WEEKNIGHTS
% ccol = 1; 
%    k = 28;  % should be number of positive eigenvalues from clustering.m
% WEEKDAYS
% ccol = 2; 
% k = 36;  % should be number of positive eigenvalues from clustering.m
% WEEKEND NIGHTS 
ccol = 3; 
k = 33;  % should be number of positive eigenvalues from clustering.m
    
    U = V(:,I(1:k));  % eigenvectors corresponding to top eigenvalues
    
    % multiply by eigenvalues
    newU = U*Dx(I(1:k),I(1:k));
    
    
    % HIERARCHICAL Clustering
    % See http://www.mathworks.com/help/stats/hierarchical-clustering.html
    
    Y = pdist(newU);  % each element contains the distance between a pair of objects.
    Z = linkage(Y); % each row is a link between objects or clusters. The first two columns identify the objects that have been linked. The third column contains the distance between these objects.
    figure()
    dendrogram(Z,nkv);  % create the dendrogram for the data with all of the kept variables
    [H,T,outperm] = dendrogram(Z,nkv);  % outperm = list of nodes in dendrogram  (kept variables)
    
    treelist = pcols(outperm); % list of nodes in dendrogram, original variables
    
    keepperm = treelist;
    
    Clustercount = smallS(outperm,outperm);  % reorganize count matrix around clusters
    
    
    
    
    
        % START creating ridgeline
        
        % cluster = in which cluster is each variable?
        cluster = 1:nkv;  % initially, every variable is in different cluster
        
        nrows = size(Z,1);  % number of rows in Z
        
        clevel = zeros(nrows,1);  % "level" from 1 (combines elements) to nrows (top)
        csize = zeros(nrows,1);  % size of cluster (number of variables)
        cstart = zeros(nrows,1);  % first element in cluster (position in outperm)
        cavgsim = zeros(nrows,1);  % average relative count (not including self)
        newheight = zeros(nrows,1);  % height in ridgline (1 = bottom)
        
        for r=1:nrows % loop over rows
            % nvk + r is the cluster number used inside Z
            elevel = [0 0]; % levels of elements
            esize = [1 1]; % sizes of elements
            epos = [0 0];  % start positions of elements
            for ei=1:2 % loop over the elements of this clusters
                if Z(r,ei) > nkv  % element is a cluster, not a variable
                    epos(ei) = cstart(Z(r,ei) - nkv);  % start position of cluster
                    esize(ei) = csize(Z(r,ei) - nkv);  % size of cluster
                    elevel(ei) = clevel(Z(r,ei) - nkv);  % level of cluster
                else
                    epos(ei) = find(outperm == Z(r,ei));  % position of variable in outperm
                    % this uses the value from Z in outperm, which has only nkv
                    % variables.
                end
            end % loop over ei
            cstart(r) = min(epos);  % first position
            csize(r) = sum(esize);  % total size of new cluster
            clevel(r) = max(elevel) + 1;  % level of new cluster
            % NOTE: must use Clustercount, which is ordered by outperm
            Subcount = Clustercount(cstart(r):cstart(r)+csize(r)-1,cstart(r):cstart(r)+csize(r)-1);  % submatrix of Clustercount for this cluster
            cavgsim(r) = sum(sum(Subcount))/(csize(r)*csize(r)-csize(r));  % average relative count
            
        end % loop over rows
        
% FROM pickclusters.m

% pickclusters.m
% pick the biggest good clusters

smin = maxtime / (300*48);  % min strength to keep (1/24 of all times)

% The key is to see if any larger, strong cluster also includes this
% cluster.

picklist = ones(1,nrows); % one entry for each cluster
% 1 = keep 
% 0 = lose

for r=nrows:-1:1  % loop over rows (clusters) from top down
    
    if picklist(r) == 0 % already lost
        continue
    end
    
    if cavgsim(r) < smin  % cluster is not strong
        picklist(r) = 0; % lose it
        continue % go to next 
    end
    
    % have a strong cluster
    
    % cstart refers to position in Clustercount 
    % csize is number of elements

    for rj=1:r-1  % loop over lower rows (smaller clusters)
        if (cstart(rj) >= cstart(r)) && (cstart(rj) < cstart(r)+csize(r))  % this cluster is inside
            picklist(rj) = 0;  % lose it
            fprintf('Losing smaller cluster %3.0f inside cluster %3.0f.\n',rj,r);
        end
    end % loop over lower rows (smaller clusters)
    
    
    
end % loop over rows

% PRINT strong, maximal clusters

fprintf('List of clusters with strength >= %10.2f in %8.0f time bins.\n',smin,maxtime/300);

for r=1:nrows  % loop over rows (clusters) from top down
    
    if picklist(r) == 0 % lost
        continue
    end
    
    fprintf('\n Cluster #%2.0f has cluster strength  = %10.4f. \n',r,cavgsim(r));
    
    for vj=0:csize(r)-1  % loop over variables
        % outperm has only values 1 to nkv
        % treelist has original variable numbers
        % fprintf('%2.0f : %2.0f (%2.0f) ',cstart(r)+vj,treelist(cstart(r)+vj),outperm(cstart(r)+vj));
        % fprintf('\n');
        fprintf('Variable %3.0f : ID %3.0f \n',treelist(cstart(r)+vj),treelist(cstart(r)+vj)-1);
        
    end
    

    
    
end % loop over rows


% NEW

% each time bin = 300 seconds;  first bin = 0;
nbins = maxtime / 300 + 1;  % number of time bins
nhours = floor(nbins / 12);  % number of complete hours
ndays = floor(nhours/24);  % number of complete days

% Note: nkv = number of non-zero variables in smallS
clusterlist = find(picklist > 0); % list of picked clusters 
nclusters = length(clusterlist);  % number of picked clusters 


% Create a matrix to go from variable number to cluster number
var2cluster = zeros(max(treelist),1);  % cluster for this variable 
for r = 1:nclusters % loop over clusters in clusterlist
    cid = clusterlist(r);  % original cluster number
    var2cluster(treelist(cstart(cid):cstart(cid)+csize(cid)-1)) = cid;
end % loop over clusters


% NOW CREATE TIMELINE

% makematrix to keep track of which variables are paired in which bins
varTimeline = zeros(max(csize(clusterlist)), nbins);

for counterSubproblem = 1:nclusters % loop over clusters
    cid = clusterlist(counterSubproblem);  % cluster number

    % list of variables in this cluster:
    cvarlist = treelist(cstart(cid):cstart(cid)+csize(cid)-1);
    
    varTimeline = zeros(max(csize(clusterlist)), nbins);    % reset to 0
    
    % loop over items in M by ID in this cluster
    for vj=1:csize(cid) 
        IDinM = cvarlist(vj) - 1;  % subtract one to get back to M
        IDlist = find(M(:,2)==IDinM);  % list of rows in M with this ID
        if isempty(IDlist)  % no rows in M with this ID
            continue   % go to next variable
        end
        
        ID2list = unique(M(IDlist,3));  % list of unique IDs with this ID
        
        for ik=1:length(ID2list)  % loop over unique IDs
            % note: ID2list(ik) + 1 is the variable number
            if var2cluster(ID2list(ik) + 1)==cid  % if this variable is in this cluster!
                % now find rows in IDlist with this ID2
                tlist = find(M(IDlist,3)==ID2list(ik));  
                % NOTE: values in tlist are pointers to IDlist (not M)
                % convert these times to bins
                binlist = M(IDlist(tlist),1)/300 + 1;  % list of bins

                varTimeline(vj,binlist) = 1;  % this variable is used in these bins
                % which variable in this cluster is ID2
                wv = find(cvarlist == ID2list(ik)+1);  
                varTimeline(wv,binlist) = 1;  % this variable is used in these bins
            end 
        end % loop over rows
    end % loop over items
    
    % now take max over hour (12 bins) for each variable

    % timeline by hour
    hourTimeline = zeros(csize(cid),nhours);
    for vj=1:csize(cid) % loop over variables in cluster
        % transpose timeline for this variable and then chop into
        % 12-hour blocks
        thisline = reshape(varTimeline(vj,1:12*nhours)',12,nhours);
        % thisline has one column for each hour (and one row for each bin)
        hourTimeline(vj,:) = max(thisline);  % max over each column
    end % loop over variables in cluster
    
    % look at the first 24*ndays hours and sum over variables
    % and then reshape the transpose
    dayTimeline = reshape(sum(hourTimeline(:,1:24*ndays))',24,ndays);
    % dayTimeline has one column for each day
    
  
    % make the heatmap with one row for each day
    if csize(cid) > 3
        figure();
        h2=heatmap(dayTimeline');
        h2.Title = strcat('Cluster #',num2str(cid));
        h2.XLabel = 'Hour';
        h2.YLabel = 'Days';
    end

end % loop over clusters


% save to appropriate column
var2clustersall(1:max(treelist),ccol) = var2cluster;

fprintf('Number of clusters: %3.0f.\n',nclusters);
fprintf('Number of persons : %3.0f.\n',sum(csize(clusterlist)));
fprintf('Largest cluster   : %3.0f.\n',max(csize(clusterlist)));

