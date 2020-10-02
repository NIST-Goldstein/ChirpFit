classdef testChirpFit < matlab.unittest.TestCase
% Class-based unit testing for chirp fitter
%
% If you are not familiar with Matlsb's Unit Test Framework you are
% missing out on a powerful feature.  See: https://www.mathworks.com/help/matlab/class-based-unit-tests.html
%
%   run these tests with two command-line commands:
%   >> testCase = testChirpFit;
%   >> res = run(testCase);
%
%   You only need to run the testCase = testChirpFit once until you clear
%   all, after editing and saving any file, you only need to run
%   res=run(testCase) to run the unit tests.
%   () 
    
    properties
        Name    % Test name
        numPhases % number of phases (waveforms) in the time series
        SignalParams    % Parameters input to PmuTestSignals
        F0
        t0              % signal start time
        AnalysisCycles
        SettlingTime
        sizeMax         % signal size
        Fs      % sample rate
        TS      % Generated Time series to be analysed
        exp     %expected values    
        fig = 0 % index used for referencing figures
    end
    
%% Signal params.  Note that the labeling convention comes mostly from the
%  standard
%     Xm = signalparams(1,:)*sqrt(2);     % phase amplitude (given by the user in RMS
%     Fin = signalparams(2,:);    % frequency (must be the same for all 6 channels or an error will be thrown
%     Ps = signalparams(3,:);     % phase 
%     Fh = signalparams(4,:);     % Frequency of the interfering signal
%     Ph = signalparams(5,:);     % Phase of the interfering signal
%     Kh = signalparams(6,:);     % index of the interfering signal    
%     Fa = signalparams(7,:);     % phase (angle) moduation frequency
%     Ka = signalparams(8,:);     % phase (angle) moduation index
%     Fx = signalparams(9,:);     % amplitude moduation frequency
%     Kx = signalparams(10,:);    % amplitude moduation index
%     Rf = signalparams(11,:);    % ROCOF (chirp rate)
%     KaS = signalparams(12,:);   % phase (angle) step index
%     KxS = signalparams(13,:);   % magnitude step index
    
% Static method to get indexes into the signalParams matrix    
    methods (Static)
        function [Xm,Fin,Ps,Fh,Ph,Kh,Fa,Ka,Fx,Kx,Rf,KaS,KxS] = getParamIndex()
            Xm=1;Fin=2;Ps=3;Fh=4;Ph=5;Kh=6;Fa=7;Ka=8;Fx=9;Kx=10;Rf=11;KaS=12;KxS=13;
        end
    end
          
%% Test Methods
% Any uncommented test functions will be run in series. These functions are
% found below in the Public Methods in the order they appear in the list.
    methods (Test)
        function regressionTests (testCase) 
            testCase.fig=0;
            %run_OAs_example (testCase);
            testTDCH (testCase)
        end
    end
    
 %% Public Methods
 methods (Access = public)
     
     function run_OAs_example (testCase)
         [~] =chirp_Parameter_Estimation (testCase);
     end
     
     function testTDCH (testCase)      
         %--------------------------------
         % the default chirp signal will be 1 second of 45 Hz steady state, 10
         % seconds of 1 Hz/second chirp, then 1 second of 55 Hz steady state.
         % We need to know the instantanious frequency and chirp rate each
         % 1/50th of a second. this includes periods when the "chirp step" is in the
         % analysis window.  The reason we use 1/50th is because we are using
         % 50 Hz as our nominal power line cycle and we will be comparing these
         % parameter estimates to the output of a PMU which is reporting at 50
         % reports per second.
         setDefaults (testCase)
         
         % Analysiscycles is the number of samples that will be in the
         % analysis window.  I chose 20 arbitrarily but we will need to
         % experimentally choose a number that gives us the kind of resolution we need.
         testCase.AnalysisCycles = 20;
         
         testCase.TS = PmuWaveforms (...
             testCase.t0,...
             testCase.SettlingTime, ...
             testCase.sizeMax, ...
             testCase.Fs,...
             testCase.SignalParams...
             );
         %------------- Debug: Visualize the waveform --------------
         %testCase.fig = testCase.fig+1
         %figure(fig)
         %t = testCase.t0-testCase.SettlingTime:1/testCase.Fs:((testCase.sizeMax-1)/testCase.Fs)+testCase.t0+testCase.SettlingTime;
         %plot(t,testCase.TS(1,:))
         %----------------------------------------------------------------
         
         % Loop through the time series one nominal cycle at a time 
         %    we need to skip the first and last half window in the
         N = length(testCase.TS);
         dur = N/testCase.Fs;
         winSize = ceil((testCase.AnalysisCycles / testCase.F0) * testCase.Fs);    %window size in samples
         sampPerCyc = testCase.Fs/testCase.F0;      % number of samples per nominal cycle
         numReports = dur*testCase.F0 - winSize/sampPerCyc;
         
         
         % storage for our results from each window of data
         act = cell(numReports+1,1);
         act(1,:) = {'XT'};
         
         % Input parameters to the TDCH
         % nn is the sample number vector squared 
         n = 0:1:winSize-1;nn=n.^2;
         rr = 2;  % rr is the maximum chirp rate that can be analyized
         P = testCase.Fs/2; % is the maximum frequency (it seems to me that this should be the Nyquist Frequency)
         step = 2*rr/P; beta = -rr:step:rr-step;
         
         %---- Comment the below if not using the plot on the loop--------
         testCase.fig = testCase.fig+1; figure(testCase.fig)
         %----------------------------------------------------------------
         idx = winSize/2;
         pause on
         for i = 1:numReports
             Y =  testCase.TS(:,idx:idx+winSize-1);
             
             % ------------- Visualize the analysis window ---------------
             plot(Y(1,:));drawnow; pause(1/25);
             % -----------------------------------------------------------

             % ----------process each window with the TDCH--------------------------
             for j = 1 : testCase.numPhases
                 [XT(:,:,i)] = TDCH(Y(1,:),beta,nn,winSize,P);
                 
                 % later on, we need to find out what the instantanious
                 % frequency and chirp rate at the center of each window
                 % is. Ideally we woudl also be able to find out the
                 % magnitude and phase of the signal at the window center.
                 % If TDCH does not provide that, we can use curve fitting
                 % techniques to estimate those.
             end
             
              % For now store all the info in a cell array to analyse later.
              act(i+1,:) = {XT};
              idx = idx + sampPerCyc;    % 1 nominal cycle is added to the end of the window, 1 cycle is removed from the end
             
         end
         pause off
         
         
            
         
         
     end
          
 end
    
%% Private Methods
    methods (Access = private)
        
        %Defaults for a 10 second linear chirp, 45 Hz to 55 Hz with 1 second steady state at each end
        function setDefaults(testCase)
            testCase.F0 = 50;
            testCase.t0 = 0; % beginning of the time series
            testCase.Fs =4800;
            testCase.SettlingTime = 1.0;
            testCase.sizeMax = testCase.Fs * 10;
            testCase.AnalysisCycles = 6;
            
            testCase.numPhases = 1;
            testCase.SignalParams = zeros (13,testCase.numPhases);
            [Xm,Fin,Ps,Fh,Ph,Kh,Fa,Ka,Fx,Kx,Rf,KaS,KxS] = testCase.getParamIndex();
            
            testCase.SignalParams(Xm,:) = 1;
            testCase.SignalParams(Fin,:) = 45;
            testCase.SignalParams(Ps,:) = 0;
            testCase.SignalParams(Rf,:) = 1;
        end    
 
    end
end