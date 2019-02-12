function mainMutRateCorr
[VDJdata, VDJheader] = openSeqData;
[Vmat, Dmat, Jmat] = findMutationFreq(VDJdata, VDJheader, 'single');
 plotMutRateCorr(Vmat, Dmat, Jmat)