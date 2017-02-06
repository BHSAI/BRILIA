%removeNAN will get rid of NaN in a cell variable and replace them with ''.
%
%  Data = removeNAN(Data)

function Data = removeNAN(Data)
for r = 1:size(Data,1)
    for c = 1:size(Data,2)
        if isnan(Data{r,c})
            Data{r,c} = '';
        end
    end
end
