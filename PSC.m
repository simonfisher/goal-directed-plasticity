%
% PSC class: defines a post-synaptic current, and measures peak
%
% Simon Fisher, 2018
%
classdef PSC < handle
    %PSC Post Synaptic Current measurement
    %   
    
    properties
        data                % PSC waveform data
        dataPSC             % data limited to 100ms around PSC event, and filtered
        file                % number of file that includes this PSC
        phase               % e.g. during baseline, or while on a drug...
        onset               % time in ms from start of stim to start of PSC
        amplitude           % measure of PSC amplitude (pA)
        slope               % measure of PSC slope (pA/ms ?)
        baseline            % current during baseline period to PSC
        peakLoc             % position of peak found
        dStim               % 'local' position of stim relative to dataPSC
        excluded            % set true if should be excluded from analyses
        
    end
    
    methods
        
        % Constructor function
        %
        %   Called from 
        %   
        %   
        function obj = PSC(data, file_num, phase)
            
            obj.data = data;
            obj.file = file_num;
            obj.phase = phase;
            
            obj.excluded = false;
            
            
            
            
        end % Constructor
        
        
        % measurePSC function
        %
        %   Called from 
        %   
        %   
        function measurePSC(obj, stimOnset, phasestr)
            
            % window data to 100 ms around opto stim
            % assume stim onset given in seconds
            stimLoc = stimOnset * 20000;
            d = obj.data(stimLoc-500:stimLoc+3000);
            dStim = 500; % set 'local'/d-relative stim location
            
            % filter data
            [B,A] = butter(2, 0.1, 'low');
            df = filtfilt(B,A,d);
            
            % convert to pA here - helps making sense of data
            % * 1e12
            df = df * 1e12;
            
            % store limited and filtered data
            % baseline normalized to zero
            obj.dataPSC = df - mean(df(1:10));
            
            % measure baseline period 20 ms prior, and 3 ms after
            basepA = mean(df(dStim-400:dStim+60));
            obj.baseline = basepA;
            % get variance for thresholding
            thresh = (max(df(dStim-400:dStim+60))-min(df(dStim-400:dStim+60)))/4;
            
            
            % Switch on experiment phase
            switch phasestr
                
                case 'minus70'
            
                    %
                    % find peak of PSC
                    %
                    
                    % restrict search area
                    df_peak = df(dStim:dStim+600);
                    
                    % call peakfinder()
                    % 'sel' set to default (to consider a peak), and threshold set
                    % to 5 pA (to count as a peak)
                    [peaksLoc,peaksMag] = peakfinder(df_peak, (max(df)-min(df))/3, 5, -1);

        %             % find relevant AP peak - check height
        %             peaksLoc_logical = d(peaksLoc+search_start) > -10;
        %             ap_peak_loc = peaksLoc(peaksLoc_logical);

                    % ensure only 1 peak is found
                    % handle number of peaks found
                    switch size(peaksLoc,2)
                        case 0
                            % zero is the value indicating no peak found
                            PSCpeakLoc = 1;

                        case 1
                            PSCpeakLoc = peaksLoc;

                        otherwise % more than 1
                            % use the biggest (which is most minima)
                            [m,i] = min(peaksMag);
                            PSCpeakLoc = peaksLoc(i);
                    end

                    % calculate PSC onset; working backwards...
                    % find first location where value is greater than
                    % baseline (considering threshold val)
                    idx = find(fliplr(df_peak(1:PSCpeakLoc)) > (basepA - (thresh/2)),1,'first');
                    % calcuate location from stim
                    PSConset = PSCpeakLoc - idx;
                    % convert to ms
                    obj.onset = PSConset / 20;


                    % calculate PSC amplitude
                    obj.amplitude = abs(df_peak(PSCpeakLoc) - basepA);
                    
                    obj.peakLoc = PSCpeakLoc / 20; % this is relative to the stim (ms)


                case 'plus40'
                    
                    % find peak
                    [peaksLoc,peaksMag] = peakfinder(df, (max(df)-min(df))/3, 1, 1);
                    
                    % ensure only 1 peak is found
                    % handle number of peaks found
                    switch size(peaksLoc,2)
                        case 0
                            % zero is the value indicating no peak found
                            PSCpeakLoc = 0;

                        case 1
                            PSCpeakLoc = peaksLoc;

                        otherwise % more than 1
                            % use the biggest (which is most maxima)
                            [m,i] = max(peaksMag);
                            PSCpeakLoc = peaksLoc(i);
                    end

                    % calculate PSC onset
                    % find first location where value is less than baseline
                    idx = find(fliplr(df(dStim:PSCpeakLoc)) < (basepA + (thresh/2)),1,'first');
                    % calcuate location from stim
                    PSConset = length(df(dStim:PSCpeakLoc)) - idx;
                    % convert to ms
                    obj.onset = PSConset / 20;
                    
                    obj.peakLoc = (PSCpeakLoc - dStim) / 20; % this is relative to the stim (ms)
                    
                    % calculate NMDA amplitude
                    % mean of current values 50-60 ms after stim
                    obj.amplitude = mean(df(dStim+1000:dStim+1200)) - basepA;
                    
                    
                    
            end % switch on exp phase
            
            
        end % measurePSC
        
        
        
        
    end
    
end

