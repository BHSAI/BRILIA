%normalizeMotifData will normalize the frequency values within the
%MotifData output file such that the sum within each cell is 1. 
%
%  NormMotifData = normalizeMotifData(MotifData)
%
%  INPUT
%    MotfiData: output from collectMotifData
%    
%  OUTPUT
%    NormMotifData: normalized MotifData

function MotifData = normalizeMotifData(MotifData)
for r = 1:size(MotifData.Data, 1)
    for c = 1:size(MotifData.Data, 2)
        CellSum = sum(MotifData.Data{r, c});
        if CellSum > 0
            MotifData.Data{r, c} = MotifData.Data{r, c} / CellSum;
        end
    end
end
