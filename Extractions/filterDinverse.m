%filterDinverse will get an VDJdata excel file, and take out just the D
%inverse results, save it into a new file.

function filterDinverse
[FileName, FilePath] = uigetfile('*.xlsx','Select the VDJ data file');

[~,~,SampleData] = xlsread([FilePath FileName]);
[SampleData, VDJheader] = filterHeader(SampleData);
H.FamNumLoc = findHeader(VDJheader,{'vMapNum','dMapNum','jMapNum'});


%Extracting location of inverse D's
ExtractThis = zeros(size(SampleData,1),1)==1;
for j = 1:size(SampleData,1)
    Dmap = SampleData{j,H.FamNumLoc(2)};
    if ischar(Dmap)
        Dmap = eval(Dmap);
    end
    %look for ANY even numbers
    for k = 1:length(Dmap)
        if mod(Dmap(k),2) == 0
            ExtractThis(j) = 1;
        end
    end
end

[SaveName, SavePath] = uiputfile('*.xlsx','Save the D inverses as',[FileName(1:end-5)]);
xlswrite([SavePath SaveName '.Dinv.xlsx'],cat(1,VDJheader,SampleData(ExtractThis,:)));
xlswrite([SavePath SaveName '.Dfwd.xlsx'],cat(1,VDJheader,SampleData(ExtractThis==0,:)));
