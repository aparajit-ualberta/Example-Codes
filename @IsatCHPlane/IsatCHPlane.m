classdef IsatCHPlane
    properties
        data = struct('t',[],'r',[],'Isat_diff',[], 'Isat_DC',[]);
    end
    
    methods
        function obj = IsatCHPlane(val1, val2, val3, val4)
            if nargin == 3
                obj.data.t = val1;
                obj.data.r = val2;
                obj.data.Isat_diff = val3;
                obj.data.Isat_DC = val4;
            else
                load('IsatR_diff.mat')
                load('RCH.mat')
                load('TCHR.mat')
                load('IsatR_DC')
                obj.data.t = t;
                obj.data.r = r;
                obj.data.Isat_diff = Isat_diff;
                obj.data.Isat_DC = Isat_DC;
            end
        end
    end
    
    methods
        [Isat_diff_temp,t_temp] = subsample(obj,varargin)
        [H,C,CHFig] = CHPlane_Posn(obj,r,varargin)
        [H,C,CH_Shot_Fig] = CHPlane_PosnRange_Shot(obj,varargin)
        [H,C,CH_Mean_Fig] = CHPlane_PosnRange_Mean(obj,varargin)
        PeakReduction(obj,rpos,varargin)
        [fig1,ax,ax1,ax2] = PeakReduction_Subfigure(obj,rpos,varargin)
        [PF] = PosnFreq_Shot(obj,r,st,varargin)
    end
end
        
        