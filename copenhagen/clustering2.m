% clustering2.m
% created 4-21-2020 by Jeffrey W. Herrmann
% this runs the spectral clustering with given value of k
% determined by output of clustering.m

    k = 17;  % should be number of positive eigenvalues from clustering.m
    
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
            
            % fprintf('Row %2.0f: start = %2.0f, size = %2.0f, level = %2.0f.\n',r,cstart(r),csize(r),clevel(r));
            
            if cavgsim(r) < 10  % average is too low
                continue  % go to next row
            end
            
            fprintf('\n Subcluster #%2.0f has cluster strength  = %10.4f. \n',r,cavgsim(r));
            
            for vj=0:csize(r)-1  % loop over variables
                % outperm has only values 1 to nkv
                % treelist has original variable numbers
                fprintf('%2.0f : %2.0f (%2.0f) ',cstart(r)+vj,treelist(cstart(r)+vj),outperm(cstart(r)+vj));
                fprintf('\n');
            end
            
        end % loop over rows
        
            