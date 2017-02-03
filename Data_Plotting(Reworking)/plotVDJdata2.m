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
if isempty(ClustThis); 
    ClustThis = 1; 
end

[Vmap,Dmap,Jmap] = getCurrentDatabase;

for f = 1:length(FileNames)
    %Open the datafile
    [VDJdata, NewHeader, FileName, ~] = openSeqData([FilePath FileNames{f}]);
    DotLoc = find(FileName == '.');
    SaveFilePre = FileName(1:DotLoc(end)-1);
    SavePath = FilePath;
    
    getHeaderVar;

    if ClustThis ~= 1
        VDJdata(:,GrpLoc) = num2cell([1:size(VDJdata,1)]'); %Gets rid of group info.
    end
    
    %==========================================================================
    %Plotting the VDJ gene usage scatterplot
    %VDJdata = renumberVDJmap(VDJdata,NewHeader,Vmap,Dmap,Jmap);
%    [VDJcombo, AxisLabels] = findGeneFamilyUsage(VDJdata,NewHeader,2,Vmap,Dmap,Jmap);
%    
%     %Select the top 10 V, top 10 D, top 4 J.
%     Jsum = zeros(size(VDJcombo,3),1);
%     for j = 1:length(Jsum)
%         Jsum(j) = sum(sum(VDJcombo(:,:,j)));
%     end
%     [Jsum, SortIdx] = sort(Jsum,'descend');
%     Jdel = SortIdx(5:end);
%     VDJcombo(:,:,Jdel) = [];
%     AxisLabels{3}(Jdel) = [];
%     
%     Dsum = zeros(size(VDJcombo,2),1);
%     for j = 1:length(Dsum)
%         Dsum(j) = sum(sum(VDJcombo(:,j,:)));
%     end
%     [Dsum, SortIdx] = sort(Dsum,'descend');
%     Ddel = SortIdx(11:end);
%     VDJcombo(:,Ddel,:) = [];
%     AxisLabels{2}(Ddel) = [];
%     
%     Vsum = zeros(size(VDJcombo,1),1);
%     for j = 1:length(Vsum)
%         Vsum(j) = sum(sum(VDJcombo(j,:,:)));
%     end
%     [Vsum, SortIdx] = sort(Vsum,'descend');
%     Vdel = SortIdx(11:end);
%     VDJcombo(Vdel,:,:) = [];
%     AxisLabels{1}(Vdel) = [];
%     
%     %Removing "IGH" from the AxisLabels
%     for q = 1:length(AxisLabels)
%         AxLabels = AxisLabels{q};
%         for qq = 1:length(AxLabels)
%             AxLabels{qq} = strrep(AxLabels{qq},'IGH','');
%         end
%         AxisLabels{q} = AxLabels;
%     end
%     
%     Gx = plotGeneFamilyCombo(VDJcombo, AxisLabels);    
%     saveas(Gx,[SavePath SaveFilePre '_VD-DJcombo.png']);

    % %==========================================================================
    % %Plotting the VDJ deletion/insertoin stats
    % AllHeaderEval = {'vDeletion', 'n2Insertion', 'd5Deletion', 'd3Deletion', 'n1Insertion', 'jDeletion'};%, 'dMatchCt'};
    % 
%     AllHeaderEval = {'vMutCt_Germline', 'dMutCt_Germline', 'jMutCt_Germline'};
%     for k = 1:length(AllHeaderEval)
%         HeaderEval = AllHeaderEval(k);
%         ColData = findGeneMods(VDJdata,NewHeader,HeaderEval,ClustThis);
%         [~, Ax1] = plotGeneMods(ColData,[0:25],'','norm');
%         modifyPlotLabels(Ax1,'','',HeaderEval,AllHeaderEval{k},[10 10 14 10]); %Use HeaderEval as Title
%         
%         %Print figure as png
%         set(gcf,'PaperPosition',[0 0 5 5]);
%         HeaderNospace = HeaderEval{1};
%         HeaderNospace(HeaderNospace == ' ') = [];
%         saveas(gcf,[SavePath SaveFilePre '_' HeaderNospace '.png']);
%     end

    %==========================================================================
    %Comparing V,D,J mutation frequencies with respect to germline
    CorrOpt = 'a';

    Option = 'none';
    [Vmat, Dmat, Jmat] = findMutationFreq(VDJdata,NewHeader);
    Gx = plotMutRateCorr(Vmat,Dmat,Jmat,CorrOpt);
    set(Gx,'PaperPosition',[0 0 7 7])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VDJ' Option '.png']);

    Option = 'ham';
    VDJdataH = buildTree(VDJdata,NewHeader,'saveoff','plotoff',Option);
    [Vmat, Dmat, Jmat] = findMutationFreq(VDJdataH,NewHeader);
    Gx = plotMutRateCorr(Vmat,Dmat,Jmat,CorrOpt);
    set(Gx,'PaperPosition',[0 0 7 7])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VDJ' Option '.png']);
    
    Option = 'shmham';
    VDJdataS = buildTree(VDJdata,NewHeader,'saveoff','plotoff',Option);
    [Vmat, Dmat, Jmat] = findMutationFreq(VDJdataS,NewHeader);
    Gx = plotMutRateCorr(Vmat,Dmat,Jmat,CorrOpt);
    set(Gx,'PaperPosition',[0 0 7 7])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VDJ' Option '.png']);

    Gx = plotMutationFreq(Vmat+Dmat+Jmat);
    set(Gx,'PaperPosition',[0 0 7 7])
    saveas(Gx,[SavePath SaveFilePre '_MutFreq-VDJ' '.png']);

%     %==========================================================================
%     %Comparing V,D,J mutation frequencies + N1/N2 compositions
%     [TDTmatrix, TDTcomp, Left2Right] = findTDTmatrix(VDJdata,NewHeader,'single'); %Doesn't divide the N regions
%     Gx2 = plotTDTmatrix(TDTmatrix, TDTcomp);
% 
%     set(Gx2,'PaperPosition',[0 0 10 5])
%     saveas(Gx2,[SavePath SaveFilePre '_TDTsingle.png']);
% 
%     [TDTmatrix, TDTcomp, Left2Right] = findTDTmatrix(VDJdata,NewHeader,'flip'); %divide the N regions
%     Gx2 = plotTDTmatrix(TDTmatrix, TDTcomp, Left2Right);
% 
%     set(Gx2,'PaperPosition',[0 0 10 5])
%     saveas(Gx2,[SavePath SaveFilePre '_TDTdivide.png']);
%     
    %==========================================================================
    %Comparing ADAR vs AID mutation correlation
    Output = findMutationCorr(VDJdataS,NewHeader); %Returns the matlab default AA order: ARNDCQEGHILKMFPSTWYV
    Gx3 = plotMutationCorr(Output);
%     
%     set(Gx3,'PaperPosition',[0 0 5 5])
    saveas(Gx3,[SavePath SaveFilePre '_MutationCorr.png']);
%     
%     %==========================================================================
%     %Comparing hot spot stats
%     [Amat, Cmat, Gmat, Tmat] = findHotSpots(VDJdata,NewHeader);
%     Gx4 = plotHotSpots(Amat, Cmat, Gmat, Tmat);
%     
%     set(Gx4{1},'PaperPosition',[0 0 20 5])
%     saveas(Gx4{1},[SavePath SaveFilePre '_HotSpots.png']);
%     
%     set(Gx4{2},'PaperPosition',[0 0 20 5])
%     saveas(Gx4{2},[SavePath SaveFilePre '_AllSpots.png']);
%     
%     
%     %Check to see if resolve pur/pyr combo...
%     Apur = sum(Amat{1}([1 3],:),1);
%     Apyr = sum(Amat{1}([2 4],:),1)
%     Cpur = sum(Cmat{1}([1 3],:),1);
%     Cpyr = sum(Cmat{1}([2 4],:),1)
%     Gpur = sum(Gmat{1}([1 3],:),1);
%     Gpyr = sum(Gmat{1}([2 4],:),1)
%     Tpur = sum(Tmat{1}([1 3],:),1);
%     Tpyr = sum(Tmat{1}([2 4],:),1)
%     
%     figure
%     subplot(2,2,1)
%     bar(Apur)
%     subplot(2,2,2)
%     bar(Cpur)
%     subplot(2,2,3)
%     bar(Gpur)
%     subplot(2,2,4)
%     bar(Tpur)
% 
% 
    close all
end