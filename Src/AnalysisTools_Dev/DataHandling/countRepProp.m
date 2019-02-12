%countRepProp will cound a property of a repertoire and return a Mx2 cell
%array of frequencies either at the clonotype or individual cell level.
%Valid properties are:
%     'TemplateCount'
%     'TreeChildCount'
%     'H-CDR1_AminoAcid'
%     'H-CDR2_AminoAcid'
%     'H-CDR3_AminoAcid'
%     'H-V_GeneName'
%     'H-D_GeneName'
%     'H-J_GeneName'
%     'H-CDR3_Length'
%     'H-Length_V'
%     'H-Length_Nvd'
%     'H-Length_D'
%     'H-Length_Ndj'
%     'H-Length_J'
%     'H-SHM_V'
%     'H-SHM_Nvd'
%     'H-SHM_D'
%     'H-SHM_Ndj'
%     'H-SHM_J'
%     'H-V_Deletion3'
%     'H-D_Deletion5'
%     'H-D_Deletion3'
%     'H-J_Deletion5'
%     'H-V_Align'
%     'H-D_Align'
%     'H-J_Align'
%
%  Freq = countRepProp(S, SizeFilter, Prop1, Prop2, ...)
%
%  INPUT
%    S: structure of VDJdata's from multiple files (see openMultSeqData.m)
%    SizeFilter ['AC' 'BC' 'BCN' 'TOPN' 'BOTN']: filter for the group based
%         on clonotype sizes
%      AC - all clonotypes
%      BC - branched clonotypes with >= 2 unique sequences per clonotype
%      BCN - branched clonotypes with >= N unique sequences per clonotype
%      TOPN - top N clonotypes based on total Template count per clonotype
%      BOTN - top N clonotpyes based on total Template count per clonotype
%  
%  OUTPUT
%    Freq: Mx2 cell where the 2nd column is the frequency and 1st column is
%      the property. The 1st column format is: [PropertyName: BinName]
%      EX: {'H-SHM_D:3' [30]}

function S = countRepProp(S, SizeFilter, varargin)
for f = 1:length(S)
    G = getGrpIdx(S(f).VDJdata, S(f).VDJheader, SizeFilter);
    Idx = arrayfun(@(x) x.Idx(1), G); 
    Weight = ones(length(Idx), 1);
    
    PropIdx = findCell(S(f).VDJheader, varargin, 'MatchCase', 'any');
    if PropIdx(1) == 0
        error('%s: Invalid or non-existent property', mfilename);
    end
    
    for w = 1:length(PropIdx)
        DataName = sprintf('%s_%s', SizeFilter, strrep(S(f).VDJheader{PropIdx(w)}, '-', ''));
        S(f).(DataName) = countData(S(f).VDJdata(Idx, PropIdx(w)), Weight);
    end
end