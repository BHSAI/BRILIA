classdef SimulationParam
    properties
        MeanV3Del;
        MeanD5Del;
        MeanD3Del;
        MeanJ5Del;
        MeanMLen;
        MeanNLen;
        DelCap
        DB;
        Freq;
    end
    
    methods (Access = 'public')
        function Obj = SimulationParam(Species)
            
            switch lower(Species)
                case 'mouse_c57bl6'
                    %These were collected from the C57BL6 mice data set from A Collins,
                    %2015, processed with old version of BRILIA, v1.9.0.
                    Obj.MeanV3Del = 1;
                    Obj.MeanD5Del = 4.5;
                    Obj.MeanD3Del = 3.4;
                    Obj.MeanJ5Del = 3.9;
                    Obj.MeanMLen = 3.8;
                    Obj.MeanNLen = 2.9;
                    Obj.Freq = openFreqData('MOUSE_C57BL6_VDJ_FREQ.csv');
                case 'mouse_balbc'
                    %These were collected from the C57BL6 mice data set from A Collins,
                    %2015, processed with old version of BRILIA, v1.9.0.
                    Obj.MeanV3Del = 1;
                    Obj.MeanD5Del = 4.5;
                    Obj.MeanD3Del = 3.4;
                    Obj.MeanJ5Del = 3.9;
                    Obj.MeanMLen = 3.8;
                    Obj.MeanNLen = 2.9;
                    Obj.Freq = openFreqData('MOUSE_BALBC_VDJ_FREQ.csv');
                case 'human'
                    %Souto-Carneiro, M.M., et al., Characterization of the Human Ig Heavy
                    %P.Chain Antigen Binding Complementarity Determining Region 3 Using a
                    %Newly Developed Software Algorithm, JOINSOLVER. The Journal of
                    %Immunology, 2004. 172(11): p. 6790-6802.
                    Obj.MeanV3Del = 2.1;
                    Obj.MeanD5Del = 4.5;
                    Obj.MeanD3Del = 5.0;
                    Obj.MeanJ5Del = 6.8;
                    Obj.MeanMLen = 7;
                    Obj.MeanNLen = 7;
                    Obj.Freq = openFreqData('HUMAN_VDJ_FREQ.csv');
                otherwise
                    error('%s: No parameter for this species "%s".', mfilename, Species);
            end
            Obj.DB = getGeneDatabase(Species);
            Obj.DelCap = 13; %Maximum deletion set for all species, for now.
        end
        
        function openFreqData(FileName)
            FreqData = readDlmFile(FileName);
            
        end
    end
end