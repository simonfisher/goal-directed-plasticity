%
% PP class: defines a paired-pulse recording and computes the ratio
%
% Simon Fisher, 2018
%

classdef PP < handle
    %
    %   
    
    properties
        data                % PSC waveform data
        dataPP              % data limited to 100ms around PSC event, and filtered
        file                % number of file that includes this PSC
        interval            % e.g. during baseline, or while on a drug...
        onsets              % time in ms from start of stim to start of PSC
        amplitudes          % measures of pulse amplitudes (pA)
        slope               % measure of PSC slope (pA/ms ?)
        baselines           % current during baseline period for pulses
        peakLocs            % peak locations of pulses measured
        excluded            % set true if should be excluded from analyses
        ppr                 % ratio of pulse 2 to pulse 1
        
    end
    
    methods
        
        % Constructor function
        %
        %   Called from 
        %   
        %   
        function obj = PP(data, file_num, intervalstr)
            
            obj.data = data;
            obj.file = file_num;
            obj.interval = intervalstr;
            
            obj.excluded = false;          
            
            
        end % Constructor
        
        
        % measurePP function
        %
        %   Called from 
        %   
        %   
        function measurePP(obj, stimOnset, intervalstr)
            
            %
            % Config params
            %
            
            sr = 20000;
            
            % time of first stim
            % assume stim onset given in seconds
            stimLoc1 = stimOnset * sr;
            
            % Switch on experiment phase
            switch intervalstr
                case 'ms50'
                    stimLoc2 = stimLoc1 + (0.05 * sr);
                case 'ms100'
                    stimLoc2 = stimLoc1 + (0.1 * sr);
            end
            
            % window data to 50 ms before first stim, and 100 ms after last stim
            d = obj.data(stimLoc1-1000:stimLoc2+2000);
            dStim(1) = 1000; % set 'local'/d-relative stim location
            dStim(2) = stimLoc2 - stimLoc1 + dStim(1);
            
            % filter data
            [B,A] = butter(2, 0.1, 'low');
            df = filtfilt(B,A,d);
            
            % convert to pA here - helps making sense of data
            % * 1e12
            df = df * 1e12;
            
            % store limited and filtered data
            obj.dataPP = df;
            
            %
            % Find peaks
            %
            [peaksLoc,peaksMag] = peakfinder(df, (max(df)-min(df))/8, 5, -1);

            % handle number of peaks found
            switch size(peaksLoc,2)
                case 0
                    % zero is the value indicating no peak found
                    PPpeaksLoc = [0,0];
                    PPpeaksMag = [0,0];

                case 1
                    % 1 peak found is a problem - likely 2nd failed
                    PPpeaksLoc = [0,0];
                    PPpeaksMag = [0,0];

                case 2
                    % ideal scenario of 2 peaks
                    PPpeaksLoc = peaksLoc;
                    PPpeaksMag = peaksMag;

                otherwise % more than 2
                    
                    % make the first, the first peak after 1st stim
                    idx = find(peaksLoc > dStim(1),1,'first');
                    PPpeaksLoc(1) = peaksLoc(idx);
                    PPpeaksMag(1) = peaksMag(idx);
                    
                    % make the second, the first after the 2nd stim
                    idx = find(peaksLoc > dStim(2),1,'first');
                    if ~isempty(idx)
                        PPpeaksLoc(2) = peaksLoc(idx);
                        PPpeaksMag(2) = peaksMag(idx);
                    else
                        PPpeaksLoc = [0,0];
                        PPpeaksMag = [0,0];
                    end
                    
            end
            
            obj.peakLocs = PPpeaksLoc;
            

            % if two peaks as expected, process further
            if length(PPpeaksLoc) == 2 && PPpeaksLoc(1) ~= 0
            
                % 
                % For each paired pulse (peak)
                %
                for n=1:length(PPpeaksLoc)

                    % Calculate baselines and thresholds for each pulse
                    % measure baseline period 10 ms prior, and 2 ms after
                    basepA(n) = median(df(dStim(n)-200:dStim(n)+40));
                    obj.baselines(n) = basepA(n);
                    % get variances for thresholding
                    thresh(n) = (max(df(dStim(n)-200:dStim(n)+40))-min(df(dStim(n)-200:dStim(n)+40)))/4;

                    % Calculate the pulse onset; working backwards...
                    % find first location where value is greater than baseline (considering threshold val)
                    idx = find(fliplr(df(1:PPpeaksLoc(n))) > (basepA(n) - (thresh(n)/2)),1,'first');
                    % calcuate location of pulse onset relative to df
                    df_idx = PPpeaksLoc(n) - idx;
                    % calculate 'true' onset from stim
                    % convert to ms
                    obj.onsets(n) = (df_idx - dStim(n)) / 20;

                    % Calculate amplitudes of each peak
                    obj.amplitudes(n) = abs(df(PPpeaksLoc(n)) - basepA(n));


                end


                % compute ratio
                obj.ppr = obj.amplitudes(2) / obj.amplitudes(1);
                
            else
                obj.baselines = [0,0];
                obj.onsets = [0,0];
                obj.amplitudes = [0,0];
                obj.ppr = 0;
            end
            
                    
                    
            
        end % measurePSC
        
    end
    
end

