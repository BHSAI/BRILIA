%MultiPanelPlot is an object designed to help plot multi-panel figures by
%intaking axes, and reassembling these into subpanels of unique dimensions
%
%
classdef MultiPanelPlot < handle
    properties
        Panel struct %stores all the information about the axes and desired location, etc
          %.FileName  %location to .fig or .png (any image file)
          %.Position  %Absolution Position you want this to be
          %.Coordinate %Subpanel Row and Col (r, c) coordinate
    end
    
    methods
        function O = MultiPanel(varargin)
            
        end
    
    end
    
    
    
    
end