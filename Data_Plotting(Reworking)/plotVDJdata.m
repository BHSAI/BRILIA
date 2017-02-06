 %plotVDJdata will plot data unique to VDJdata format, such as N1 N2 comp,
%mutation frequency. Will make a new folder, named after the orig final,
%saving all the necessary plots.

%FIGURE GUIDELINES
%Most journals accept 3.3" for single column, 6.9" for double column. 

[FileNames, FilePath] = uigetfile('*.xlsx','Select all files to plot data for','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

%Determine if you want to find frequency per cluster or single sequences
ClustThis = input('Do you want to cluster the data? 0 for no, 1 for yes (default)');

[Vmap,Dmap,Jmap] = getCurrentDatabase;

for f = 1:length(FileNames)
    %Open the datafile
    [VDJdata, VDJheader, FileName, ~] = openSeqData([FilePath FileNames{f}]);
    DotLoc = find(FileName == '.');
    SaveFilePre = FileName(1:DotLoc(end)-1);
    SavePath = FilePath;
    % [SaveFilePre, SavePath] = uiputfile('*.*','Set the savefile Prefix and Folder',SaveFilePre);

    % VDJdata = buildTree(VDJdata,VDJheader,2);

    GrpLoc = findHeader(VDJheader,'GroupNum');
    if isempty(ClustThis); ClustThis = 1; end
    if ClustThis ~= 1
        VDJdata(:,GrpLoc) = num2cell([1:size(VDJdata,1)]'); %Gets rid of group info.
    end

    
    %==========================================================================
    %Plotting the VDJ gene usage scatterplot
    %VDJdata = renumberVDJmap(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    [VDJcombo, AxisLabels] = findGeneFamilyUsage(VDJdata,VDJheader,1,Vmap,Dmap,Jmap);
    
    %Select the top 10 V, top 10 D, top 4 J.
    Jsum = zeros(size(VDJcombo,3),1);
    for j = 1:length(Jsum)
        Jsum(j) = sum(sum(VDJcombo(:,:,j)));
    end
    [Jsum, SortIdx] = sort(Jsum,'descend');
    Jdel = SortIdx(5:end);
    VDJcombo(:,:,Jdel) = [];
    AxisLabels{3}(Jdel) = [];
    
    Dsum = zeros(size(VDJcombo,2),1);
    for j = 1:length(Dsum)
        Dsum(j) = sum(sum(VDJcombo(:,j,:)));
    end
    [Dsum, SortIdx] = sort(Dsum,'descend');
    Ddel = SortIdx(11:end);
    VDJcombo(:,Ddel,:) = [];
    AxisLabels{2}(Ddel) = [];
    
    Vsum = zeros(size(VDJcombo,1),1);
    for j = 1:length(Vsum)
        Vsum(j) = sum(sum(VDJcombo(j,:,:)));
    end
    [Vsum, SortIdx] = sort(Vsum,'descend');
    Vdel = SortIdx(11:end);
    VDJcombo(Vdel,:,:) = [];
    AxisLabels{1}(Vdel) = [];
    
    
    Gx = plotGeneFamilyCombo2(VDJcombo, AxisLabels);    
    set(Gx,'PaperPosition',[0 0 1*size(VDJcombo,2) 1*size(VDJcombo,3)+size(VDJcombo,1)])
    saveas(Gx,[SavePath SaveFilePre '_VD_DJcombo.png']);
    
    
    continue
    
%     saveas(Gx{1},[SavePath SaveFilePre '_DJcombo.png']);
%     saveas(Gx{2},[SavePath SaveFilePre '_VDcombo.png']);
%     saveas(Gx{3},[SavePath SaveFilePre '_VJcombo.png']);

    % %==========================================================================
    % %Plotting the VDJ deletion/insertoin stats
    % AllHeaderEval = {'vDeletion', 'n2Insertion', 'd5Deletion', 'd3Deletion', 'n1Insertion', 'jDeletion'};%, 'dMatchCt'};
    % 
    % for k = 1:length(AllHeaderEval)
    %     HeaderEval = AllHeaderEval(k);
    %     ColData = cell(1,length(FileName));
    %     for j = 1:length(FileName)
    %         ColData(j) = findGeneMods(FileData{j,1},FileData{j,2},HeaderEval,ClusterThis);
    %     end
    %     [~, Ax1] = plotGeneMods(ColData,[0:25],'','norm');
    %     modifyPlotLabels(Ax1,'','',HeaderEval,FileLegend,[10 10 14 10]); %Use HeaderEval as Title
    %     
    %     %Print figure as png
    %     set(gcf,'PaperPosition',[0 0 5 5]);
    %     HeaderNospace = HeaderEval{1};
    %     HeaderNospace(HeaderNospace == ' ') = [];
    %     saveas(gcf,[OutputPath OutputFile '_' HeaderNospace '.png']);
    % end


    %==========================================================================
    %Comparing V,D,J mutation frequencies + N1/N2 compositions
    Option = 'single';
    [Vmat, Dmat, Jmat, VDJlen] = findMutationFreq(VDJdata,VDJheader,Option);
    Gx = plotMutRateCorr(Vmat,Dmat,Jmat,'a');

    %Gx = plotMutationFreq(Vmat, Dmat, Jmat, Vmat+Dmat+Jmat,'Legend',{'V','D','J','V+D+J'});
    set(Gx,'PaperPosition',[0 0 10 10])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VDJ' Option '.png']);

    Gx = plotMutationFreq(Vleftmat,CDR3mat,'Legend',{'V before 104C','CDR3'});
    set(Gx,'PaperPosition',[0 0 10 5])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VCDR3' Option '.png']);

    %==========================================================================
    %Comparing V,D,J mutation frequencies + N1/N2 compositions
    [TDTmatrix, TDTcomp, Left2Right] = findTDTmatrix(VDJdata,VDJheader,'single'); %Doesn't divide the N regions
    Gx2 = plotTDTmatrix(TDTmatrix, TDTcomp);

    set(Gx2,'PaperPosition',[0 0 10 5])
    saveas(Gx2,[SavePath SaveFilePre '_TDTsingle.png']);

    [TDTmatrix, TDTcomp, Left2Right] = findTDTmatrix(VDJdata,VDJheader,'flip'); %divide the N regions
    Gx2 = plotTDTmatrix(TDTmatrix, TDTcomp, Left2Right);

    set(Gx2,'PaperPosition',[0 0 10 5])
    saveas(Gx2,[SavePath SaveFilePre '_TDTdivide.png']);

    % CleavageType = {'AG-CT'; 'AC-GT'; 'AT-CG'; 'CT-AG'; 'GT-AC'; 'GC-AT'};
    % for j = 1:size(CleavageType,1)
    %     CleaveNT = regexp(CleavageType{j},'-','split');
    %     LeftNTs = CleaveNT{1};
    %     RightNTs = CleaveNT{2};
    %         
    %     [TDTmatrix, TDTcomp] = findTDTmatrix(VDJdata,VDJheader,'flip',LeftNTs,RightNTs); %divide the N regions
    %     Gx2 = plotTDTmatrix(TDTmatrix, TDTcomp);
    % 
    %     set(Gx2,'PaperPosition',[0 0 10 5])
    %     saveas(Gx2,[OutputPath OutputFile '_TDTdivide' CleavageType{j} '.png']);
    % end

    % 
    % if ClustThis == 1 && length(Gx)>1
    %     set(Gx{2},'PaperPosition',[0 0 10 10])
    %     saveas(Gx{2},[OutputPath OutputFile '_MutFreq-ClustGT1.png']);
    % end
    % 

    % 
    % 
    % %==========================================================================
    % %Comparing ADAR vs AID mutation correlation
    % Output = findMutationCorr(VDJdata,VDJheader); %Returns the matlab default AA order: ARNDCQEGHILKMFPSTWYV
    % Gx3 = plotMutationCorr(Output);
    % 
    % set(Gx3,'PaperPosition',[0 0 5 5])
    % saveas(Gx3,[SavePath SaveFilePre '_MutationCorr.png']);
    % % 
    % 
    % %==========================================================================
    % %Comparing hot spot stats
    % [Amat, Cmat, Gmat, Tmat] = findHotSpots(VDJdata,VDJheader); %Returns the matlab default AA order: ARNDCQEGHILKMFPSTWYV
    % Gx4 = plotHotSpots(Amat, Cmat, Gmat, Tmat);
    % 
    % set(Gx4{1},'PaperPosition',[0 0 20 5])
    % saveas(Gx4{1},[SavePath SaveFilePre '_HotSpots.png']);
    % 
    % set(Gx4{2},'PaperPosition',[0 0 20 5])
    % saveas(Gx4{2},[SavePath SaveFilePre '_AllSpots.png']);
    % 
    % 
    % %Check to see if resolve pur/pyr combo...
    % Apur = sum(Amat{1}([1 3],:),1);
    % Apyr = sum(Amat{1}([2 4],:),1)
    % Cpur = sum(Cmat{1}([1 3],:),1);
    % Cpyr = sum(Cmat{1}([2 4],:),1)
    % Gpur = sum(Gmat{1}([1 3],:),1);
    % Gpyr = sum(Gmat{1}([2 4],:),1)
    % Tpur = sum(Tmat{1}([1 3],:),1);
    % Tpyr = sum(Tmat{1}([2 4],:),1)
    % 
    % figure
    % subplot(2,2,1)
    % bar(Apur)
    % subplot(2,2,2)
    % bar(Cpur)
    % subplot(2,2,3)
    % bar(Gpur)
    % subplot(2,2,4)
    % bar(Tpur)


    %==========================================================================
    %Comparing V,D,J mutation frequencies + N1/N2 compositions
    FreqMat = findCDR3MutFreq(VDJdata,VDJheader); %Returns the matlab default AA order: ARNDCQEGHILKMFPSTWYV
    % 
    % %Based on amino acid VDV properties
    % %http://www.proteinsandproteomics.org/content/free/tables_1/table08.pdf
    % %Based on HPI
    % %Kyte-Doolittle paper
    % 
    [~, ~, AAprop] = xlsread('AminoAcidProp.xlsx');
    [AAprop, AAheader] = filterHeader(AAprop);

    %Order by VDW
    AApropVDW = sortrows(AAprop,[3 2]);
    AAorderVDW = cell2mat(AApropVDW(:,2))';

    %Order by HPI
    AApropHPI = sortrows(AAprop,[4 2]);
    AAorderHPI = cell2mat(AApropHPI(:,2))';

    Gx1 = plotCDR3MutFreq(FreqMat,AAorderVDW);
    set(Gx1,'PaperPosition',[0 0 10 10])
    saveas(Gx1,[SavePath SaveFilePre '_AAMutFreq_VDW.png']);
    % 
    % Gx2 = plotCDR3MutFreq(FreqMat,AAorderHPI);
    % set(Gx2,'PaperPosition',[0 0 10 10])
    % saveas(Gx2,[OutputPath OutputFile '_AAMutFreq_HPI.png']);
    
    close all
end
