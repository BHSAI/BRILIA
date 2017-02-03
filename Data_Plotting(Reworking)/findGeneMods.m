%findGeneMods will analyzes Adap-formated files, and look for the deletion
%and insertion frequencies. If Cluster = On, it will only report the
%cluster's results.

function Ndata = findGeneMods(SampleData,Header,HeaderEval,ClusterOn)
%Cluster the data if needed
if ClusterOn == 1
    GrpCol = findHeader(Header,'GroupNum');
    [~,GrpLoc,~] = unique(cell2mat(SampleData(:,GrpCol)));
    SampleData = SampleData(GrpLoc,:);
end

%Look for the column and data
Ndata = cell(1,length(HeaderEval)); %{Ncount Edges}
for j = 1:length(HeaderEval)
    ColNum = findHeader(Header,HeaderEval{j});
    Ndata{j} = cell2mat(SampleData(:,ColNum));
end


