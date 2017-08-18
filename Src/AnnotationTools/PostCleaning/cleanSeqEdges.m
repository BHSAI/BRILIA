%cleanSeqEdges will trim the sequence data in the BRILIA output file
%sequences according to the results from analyzeSeqEdges, adjusting all
%necessary alignment information.
function VDJdata = cleanSeqEdges(VDJdata, VDJheader, EdgeAnalysis)

