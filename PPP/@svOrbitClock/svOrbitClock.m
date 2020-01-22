classdef svOrbitClock < handle


    properties

        orbMode % whether to use precise or broadcast orbit
        clkMode % whether to use precise or broadcast clock
        
        PClock % precise clock data
        PEph   % precise orbit data
        
        BEph   % broadcast orbit and clock data
        
        iono   % ionospheric data- could be TEC map
        
        atx    % Antenna phase center file from IGS- parased
    end


    methods
        function obj = svOrbitClock()
            if nargin < 1
                % if no values are provided (this is the only option right
                % now)
                obj.orbMode = 'PRECISE';
                obj.clkMode = 'PRECISE';
            else
                % if values are provided
%                 obj.QVals = QVals;
%                 obj.QVals = QVals;
            end
        end
    end
    
    % function signatures
    methods
        % Load precise orbit and clock
        initPEph(obj,varargin)
        % Interpolate precise orbit
        [svPos,svVel,iod] = propagate(obj,prns,constInds,epochs,settings,varargin)
        % Load precise clock
        initPClock(obj,year,doy,settings,varargin)
        % Interpolate precise clock
        cbias = clock(obj,prns,constInds,epochs);
        % Load IGS TEC map
        initIonoData(obj,year,doy,settings,varargin);
        % Iono delay from TEC map
        [tecs,delays,tecSlant] = ionoDelay(obj,epoch,llh,varargin)    
        % Load atx file
        initAtxData(obj,filenameAtx);
    end


end















