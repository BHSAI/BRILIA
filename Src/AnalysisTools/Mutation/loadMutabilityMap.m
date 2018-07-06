%loadMutabilityMap will return the mutability map as an Java object. Each
%key is a tri-nucleotide, and the value is a 1x3 mutability index for the
%1st, 2nd, and 3rd posiiton.
%
%  REF - H and L of Mouse) Shaprio et al, 2002
%
%  mutMap = loadMutabilityMap;
%
%  OUTPUT:
%    mutMap: key-value mapping of tri-nucleotides to the 1st, 2nd, and 3rd
%      nt's mutability index.
%
%  See also calcSeqMutability
function mutMap = loadMutabilityMap(Chain)
%Get the map
MutData = readDlmFile('MutabilityIndex_Mouse_Shapiro2002.csv', 'Delimiter', ',');
MutData(1, :) = [];
for c = 2:size(MutData, 2)
    for r = 1:size(MutData, 1)
        try
            MutData{r, c} = convStr2NumMEX(MutData{r, c});
        catch
        end
    end
end

%Store the key-value pair in a java map object
mutMap = containers.Map; %Java object so using Java naming convention lowercase 
if strcmpi(Chain, 'H')
    for k = 1:size(MutData, 1)
        mutMap(MutData{k, 1}) = cell2mat(MutData(k, 2:4));
    end
elseif strcmpi(Chain, 'L')
    for k = 1:size(MutData, 1)
        mutMap(MutData{k, 1}) = cell2mat(MutData(k, 5:7));
    end
end
