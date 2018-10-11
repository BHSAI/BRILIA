function Combo = combvec2(varargin)
Dim = cellfun('length', varargin);
Tot = prod(Dim(:), 1);
Combo = zeros(Tot, nargin);
for j = 1:nargin
    Val = varargin{j};
    RepCt = prod(Dim(j+1:end));
    RepVec = ones(length(Val), 1) * RepCt;
    RepVal = repmat(repelem(Val(:), RepVec), Tot/sum(RepVec), 1);
    Combo(:, j) = RepVal;
end
Combo = Combo'; %TO match with combvec output