%calcDistribDist will compute the pairwise distance matrix for
%distrubutions, in which the distance matrix reveals how divergent two
%distributions are from each other.
%
%  DistMat = calcDivergeDist(Ydata)
%
%  DistMat = calcDivergeDist(Ydata, 'Method', Method)
%
%  INPUT
%    Ydata: MxN matrix containing M value entries for N distributions
%      Example, Ydata = [ 0.1 0.2;  %2 distributions, 4 bins.
%                         0.2 0.4;
%                         0.3 0.2;
%                         0.4 0.4]
%    Method: method for calculating divergence between two distributions
%      'overl': (default) Overlapping distance between two distrubtions
%        Distance = 1 -  Overlap / (Total - Overlap), where Overlap
%        is the minimum overlap value between two entry values. 
%      'bhatt': Bhattacharyya distance between two distributions
%        Distance = -ln( sum( (Px * Qx)^0.5 ) )
%  
%  OUTPUT
%    DistMat: NxN pairwise distance for the N distributions

function DistMat = calcDistribDist(Ydata, varargin)
P = inputParser;
addParameter(P, 'Method', 'overl', @(x) ischar(x) && ismember(lower(x), {'overl', 'bhatt'}));
parse(P, varargin{:})
Method = P.Results.Method;

switch Method(1:5)
    case 'overl' %Overlap
        DistMat = zeros(size(Ydata, 2)); %Similarity score
        for j = 1:size(Ydata, 2) - 1
            for k = j+1:size(Ydata, 2)
                %Overlapping area calculation
                AUC = sum(min(Ydata(:, [j, k]), [], 2));
                DC = 1 - AUC / (sum(sum(Ydata(:, [j, k]))) - AUC); %Normalize by total template count - overlap area.
                DistMat(j, k) = DC;
                DistMat(k, j) = DC;
            end
        end
    case 'bhatt' %Bhattacharyya
        DistMat = zeros(size(Ydata, 2)); %Similarity score
        for j = 1:size(Ydata, 2)-1
            Y1 = Ydata(:, j)/sum(Ydata(:, j));
            for k = j+1:size(Ydata, 2)
                Y2 = Ydata(:, k)/sum(Ydata(:, k));
                DC = - log ( sum((Y1.*Y2).^0.5) );
                DistMat(j, k) = DC;
                DistMat(k, j) = DC;
            end
        end
end
