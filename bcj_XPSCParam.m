function [Aout] = bcj_XPSCParam(trace,  modeselect, pulseParam, sampleRate, showFig, saveFig)

%Determines parameters (onset, risetime, slope, amplitude, half-width, decay, charge) of the given trace/waveform.
%All parameters also hold the timestamp(s) or amplitude as a extra value to keep track at which timepoint/amplitude the parameter was
%taken, saved at second and later positions in the vector. Can be used to
%plot the positions in the waveform. XPSC refers to detection of both
%inward and outward current, which depending on the internal pipet solution
%ie reversal potential used is EPSC or IPSC, hence the X= Excitatory or    
%Inhibitory.
%
%Input arguments:
%trace = numeric vector of the trace
%
%modeselect = scalar, 0: inward current only present, 1: inward and outward current present, 2: outward current only present, set default 0. 
%       Optional: if input for modeselect is a vector of more than 1 element: Second element is detection threshold, Third element is baseline period in ms (default is 50 ms)
%       For example; modeselect = [0, 5, 100]; 
%       0: Inward current only present is assumed.
%       5: Signals above mean - 5 * standard deviation(SD) are considered an
%       actual signal. 
%       100: The mean and SD are based on the first 100 ms of the trace.
%
%pulseParam = vector consisting of the following parameters: [delay pulse (in ms), width pulse (in ms)]
%       For example; pulseParam = [100, 1];
%       100: stimulation pulse was given at 100 ms after the first
%       datapoint in the trace. If multiple stimulations were given, you
%       can re-run bcj_XPSCParam with different delay pulse lenghts.
%       However be aware that currently the timespan used to detect a XPSC
%
%sampleRate = samplerate at which trace was recorded, set default 10 kHz
%(10000)
%
%showFig = 0: do not plot figure, 1: plot figure;
%
%saveFig = string; print and save option. Input is save name string. 
%       Plot will be saved as a .eps (vectorized image for use in eg. 
%       illustrator).

%First initiate parameters and convert them to the correct units
%(eg. datapoints vs. time)
if isempty(sampleRate)
    sampleRate = 10000;
end
blPeriod = 50; % baseline period in ms to measure background noise and offset
if isempty(modeselect)
    modeselect = [0 3];
  elseif length(modeselect) == 1
    modeselect = modeselect;
    detectThreshold = 3; % in order to detect a XPSC signal has to be beyond 5x SD+mean
  elseif length(modeselect) == 2
    detectThreshold = modeselect(2);
    modeselect= modeselect(1);
  elseif length(modeselect) == 3
    detectThreshold = modeselect(2);
    blPeriod = modeselect(3);
    modeselect=modeselect(1);
end

if isempty(pulseParam)
    pulseParam = [(100/1000)*sampleRate 1]; % default is 100 ms after start(1 ms)
else
    pulseParam(1) = (pulseParam(1)/1000)*sampleRate;
end


XPSCPresentPeriod = (50/1000)*sampleRate ; %time period during which a signal above noise level should be present, default is 50 ms

OutwardPresentPeriod = (300/1000)*sampleRate ; %time period during which the outwardCurrent is quantified
baselinePeriod = (blPeriod/1000)*sampleRate; % a period of 50 ms prior the pulseOnset is used to calculate the baseline level


% Check for signal above noise and normalize trace first
baseline = mean(trace(pulseParam(1)-baselinePeriod: pulseParam(1)));
trace = trace - baseline; % normalize trace

baselineSD = std(trace(pulseParam(1)-baselinePeriod: pulseParam(1)));

%get the maximal/minimal signal within signal period
traceMax = max(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)));
traceMin = min(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)));

if detectThreshold*baselineSD > abs(traceMax) || detectThreshold*baselineSD > abs(traceMin)
    paramStruct.onset = NaN;
    paramStruct.risetime = NaN;
    paramStruct.slope = NaN;
    paramStruct.amplitude = NaN;
    paramStruct.halfwidth = NaN;
    paramStruct.decay = NaN;

end

% Detect presence of inward and outward current
inwardDetectWindow = (15/1000)*sampleRate; %detect inward currents within a 10 ms window after pulse onset
outwardDetectWindow = (20/1000)*sampleRate; % detect outward currents within a 20 ms window after pulse onset

inwardPresent = find(trace(pulseParam(1):pulseParam(1)+inwardDetectWindow) <= (baselineSD*detectThreshold)*-1,1,'first');
outwardPresent = find(trace(pulseParam(1):pulseParam(1)+outwardDetectWindow) >= baselineSD*detectThreshold*2,1,'first'); % to make sure no false outward, multiply detectthreshold by 2

if ~isempty(inwardPresent) && isempty(outwardPresent)
    modeselect= 0;
elseif ~isempty(inwardPresent) && ~isempty(outwardPresent)
    %Turns out that sometimes "outwardCurrents" ie noise, preceed
    %inwardCurrent. To avoid these instances to result in modeselect=1 I
    %added an new conditional for which outwardCurrents cannot preceed
    %inwardCurrent
    if abs(outwardPresent)<abs(inwardPresent)
        %in this case the trace probably has a big noise artefact before
        %the inward current. So look only at the inward current!
        modeselect = 2;
    else
        modeselect= 1;
    end
elseif isempty(inwardPresent) && ~isempty(outwardPresent)
    modeselect= 2;
else
    %warning('trace does not contain a quantifiable inward/outward current')
    'trace does not contain a quantifiable inward/outward current'
    paramStruct.onset = NaN;
    paramStruct.risetime = NaN;
    paramStruct.slope = NaN;
    paramStruct.amplitude = NaN;
    paramStruct.halfwidth = NaN;
    paramStruct.decay = NaN;
    modeselect= 3;
end

% Determine parameters based on the type of signal 

switch modeselect
    case 0 % inward signal only
    amplitude = [traceMin find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin,1,'first')]; % relative to pulse start

    onsetTime = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin*0.1,1,'first')+pulseParam(1)-1; % relative to trace start
    onset = [((onsetTime-pulseParam(1))/sampleRate)*1000 trace(onsetTime)];  % relative to pulse start

    rise20 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin*0.2,1,'first')+pulseParam(1)-1; % relative to trace start
    rise80 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin*0.8,1,'first')+pulseParam(1)-1; % relative to trace start
    risetime = ((rise80 -  rise20)/sampleRate)*1000; %in ms

    slope = (trace(rise80)-trace(rise20))/risetime; % in pA/ms

    rise50 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin*0.5,1,'first')+pulseParam(1)-1; % relative to trace start
    decay50 = find(trace(pulseParam(1)+amplitude(2):(pulseParam(1)+XPSCPresentPeriod)) >= traceMin*0.5,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    halfwidth = ((decay50 - rise50)/sampleRate)*1000;

    decay20 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMin*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay20) %in some cases the decaying phase is longer than the XPSCPresentPeriod, in such case just use the whole length of the trace
        decay20 = find(trace(amplitude(2)+pulseParam(1):end) >= traceMin*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end

    decay80 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMin*0.8,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    decaytime = ((decay20 -  decay80)/sampleRate)*1000;

    %%
    %determine slopes based on derivative of the trace, the peak of derivative trace = slope, within a interval from onset of the signal till 2 or 3 ms after onset

    traceDervInterval2ms = gradient(trace(onsetTime:(onsetTime+0.002*sampleRate)));
    traceDervInterval3ms = gradient(trace(onsetTime:(onsetTime+0.003*sampleRate)));

    slope2ms = min(traceDervInterval2ms); %in pA/ms
    slope3ms = min(traceDervInterval3ms); %in pA/ms

    %calculate the total charge transfer between onset 10% of peak and decay 10% of peak
    decay10 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMin*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay10) %in some cases the decaying phase is longer than the XPSCPresentPeriod, in such case just use the whole length of the trace
        decay10 = find(trace(amplitude(2)+pulseParam(1):end) >= traceMin*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end

    area = trapz(trace(onsetTime:decay10))/(sampleRate/1000);%in pA*ms

    %%
    %Calculate the amplitude based on a sliding averging window of 2 ms
    amplitudeSlide = bcj_slidingPeak(trace, sampleRate, 2, [pulseParam(1) (pulseParam(1)+XPSCPresentPeriod) blPeriod], 1);

    %%

    onsetTimeInward = '';

    case 1 % inward and outward signal
    amplitudeInward = [traceMin find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) <= traceMin,1,'first')+pulseParam(1)-1];% relative to trace start
    amplitude = [traceMax find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax,1,'first')]; % relative to pulse start

    onsetTimeInward = find(trace(pulseParam(1):amplitudeInward(2)) >= traceMin*0.1,1,'last')+pulseParam(1)-1; % relative to trace start
    % in some cases the onsetTimeInward cannot be found, whereas a
    % amplitudeInward is defined. In that case use the current at
    % pulseonset as an offset to add to the traceMin
    if isempty(onsetTimeInward)
      onsetTimeInward = find(trace(pulseParam(1):amplitudeInward(2)) >= traceMin*0.1+trace(pulseParam(1)),1,'last')+pulseParam(1)-1; % relative to trace start  
    end
    onsetInward = [onsetTimeInward-pulseParam(1) trace(onsetTimeInward)];  % relative to pulse start
    onsetTime = find(trace(amplitudeInward(2):(pulseParam(1)+XPSCPresentPeriod)) >= (amplitudeInward(1)+(traceMax-amplitudeInward(1))*0.1),1,'first')+amplitudeInward(2)-1;
    onset = [((onsetTime - onsetTimeInward)/sampleRate)*1000 trace(onsetTime)]; %relative to inward onset

    rise20 = find(trace(amplitudeInward(2):(pulseParam(1)+XPSCPresentPeriod)) >= (amplitudeInward(1)+(traceMax-amplitudeInward(1))*0.2),1,'first')+amplitudeInward(2)-1; % relative to trace start
    rise80 = find(trace(amplitudeInward(2):(pulseParam(1)+XPSCPresentPeriod)) >= (amplitudeInward(1)+(traceMax-amplitudeInward(1))*0.8),1,'first')+amplitudeInward(2)-1;% relative to trace start
    risetime = ((rise80 -  rise20)/sampleRate)*1000;

    slope = (trace(rise80)-trace(rise20))/risetime;

    %rise50 = find(trace(amplitudeInward(2):(pulseParam(1)+XPSCPresentPeriod)) >= (amplitudeInward(1)+(traceMax-amplitudeInward(1))*0.5),1,'first')+amplitudeInward(2)-1; % relative to trace start
    %decay50 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= (amplitudeInward(1)+(traceMax-amplitudeInward(1))*0.5),1,'first')+amplitude(2)+pulseParam(1)-1;% relative to trace start
    rise50 = find(trace(amplitudeInward(2):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.5,1,'first')+amplitudeInward(2)-1; % relative to trace start
    decay50 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.5,1,'first')+amplitude(2)+pulseParam(1)-1;% relative to trace start

    halfwidth = ((decay50 -  rise50)/sampleRate)*1000;

    decay20 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay20)
        decay20 = find(trace(amplitude(2)+pulseParam(1):end) <= traceMax*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end
    decay80 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.8,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    decaytime = ((decay20 -  decay80)/sampleRate)*1000;

    %%
    %determine slopes based on derivative of the trace, the peak of derivative trace = slope, within a interval from onset of the signal till 2 or 3 ms after onset

    traceDervInterval2ms = gradient(trace(onsetTime:(onsetTime+0.002*sampleRate)));
    traceDervInterval3ms = gradient(trace(onsetTime:(onsetTime+0.003*sampleRate)));

    slope2ms = max(traceDervInterval2ms); %in pA/ms
    slope3ms = max(traceDervInterval3ms); %in pA/ms

    %calculate the total charge transfer between onset 10% of peak and decay 10% of peak
    decay10 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay10)
        decay10 = find(trace(amplitude(2)+pulseParam(1):end) >= traceMax*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end

    area = trapz(trace(onsetTime:decay10))/(sampleRate/1000);%in pA*ms
    %%
    %Calculate the amplitude based on a sliding averging window of 2 ms
    amplitudeSlide = bcj_slidingPeak(trace, sampleRate, 2, [pulseParam(1) (pulseParam(1)+XPSCPresentPeriod) blPeriod], 2);
    %%

    case 2 % outward signal only
    amplitude = [traceMax find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax,1,'first')]; % relative to pulse start
    onsetTime = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.1,1,'first')+pulseParam(1)-1; % relative to trace start
    onset = [((onsetTime-pulseParam(1))/sampleRate)*1000 trace(onsetTime)];  % relative to pulse start

    rise20 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.2,1,'first')+pulseParam(1)-1; % relative to trace start
    rise80 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.8,1,'first')+pulseParam(1)-1; % relative to trace start
    risetime = ((rise80 -  rise20)/sampleRate)*1000;

    slope = (trace(rise80)-trace(rise20))/risetime;

    rise50 = find(trace(pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.5,1,'first')+pulseParam(1)-1; % relative to trace start
    decay50 = find(trace(pulseParam(1)+amplitude(2):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.5,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start

    halfwidth = ((decay50 - rise50)/sampleRate)*1000;

    decay20 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay20)
        decay20 = find(trace(amplitude(2)+pulseParam(1):end) <= traceMax*0.2,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end
    decay80 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+OutwardPresentPeriod)) <= traceMax*0.8,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    decaytime = ((decay20 -  decay80)/sampleRate)*1000;

    %%
    %determine slopes based on derivative of the trace, the peak of derivative trace = slope, within a interval from onset of the signal till 2 or 3 ms after onset

    traceDervInterval2ms = gradient(trace(onsetTime:(onsetTime+0.002*sampleRate)));
    traceDervInterval3ms = gradient(trace(onsetTime:(onsetTime+0.003*sampleRate)));

    slope2ms = max(traceDervInterval2ms); %in pA/ms
    slope3ms = max(traceDervInterval3ms); %in pA/ms

    %calculate the total charge transfer between onset 10% of peak and decay 10% of peak
    decay10 = find(trace(amplitude(2)+pulseParam(1):(pulseParam(1)+XPSCPresentPeriod)) >= traceMax*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    if isempty(decay10)
        decay10 = find(trace(amplitude(2)+pulseParam(1):end) >= traceMax*0.1,1,'first')+pulseParam(1)+amplitude(2)-1; % relative to trace start
    else
    end
    area = trapz(trace(onsetTime:decay10))/(sampleRate/1000);%in pA*ms

    %%
    %Calculate the amplitude based on a sliding averging window of 2 ms
    amplitudeSlide = bcj_slidingPeak(trace, sampleRate, 2, [pulseParam(1) (pulseParam(1)+XPSCPresentPeriod) blPeriod], 2);
    %%

    onsetTimeInward = '';

    case 3 % no signal present
    paramStruct.onset = NaN;
    paramStruct.risetime = NaN;
    paramStruct.slope = NaN;
    paramStruct.amplitude = NaN;
    paramStruct.halfwidth = NaN;
    paramStruct.decay = NaN;
    paramStruct.slope2ms = NaN;
    paramStruct.slope3ms = NaN;
    paramStruct.area = NaN;
    paramStruct.amplitudeSlide = NaN;
    Aout = paramStruct;
    return
end

%Save parameters to structure
modeString = {'EPSC-inward','EPSC/IPSC-in/outwardmix','IPSC-outward'};

paramStruct.mode= modeString{modeselect+1};
paramStruct.onset = onset; % in ms, at 10% of the peak, for EPSC is relative to stimulation onset, for IPSC is relative to onset EPSC in case a clear inward EPSC is present, otherwise IPSC onset is also relative to stimulation.
paramStruct.risetime = risetime; % in ms 20% - 80% of peak
paramStruct.slope = slope; % in pA/ms, based on 20% - 80% of peak
paramStruct.amplitude = amplitude; % in pA
paramStruct.halfwidth = halfwidth; % in ms
paramStruct.decay = decaytime; % in ms 20% - 80% of peak
paramStruct.slope2ms = slope2ms; % in pA/ms, based on max(outward)/min(inward) of derivative within 2 ms of onset
paramStruct.slope3ms = slope3ms; % in pA/ms, based on max(outward)/min(inward) of derivative within 3 ms of onset
paramStruct.area = area; %charge transfer between onset 10% of peak till decay 10% of peak
paramStruct.amplitudeSlide = amplitudeSlide; % in pA used with an averaging window of 2ms wide


%generate parameters for plots such that all the determined waveform
%parameters are shown in the waveform as markers.
if isempty(onsetTimeInward)
    onsetTimeInward = 0;
    onsetInwardAmplitude = 0;
else
    onsetInwardAmplitude = trace(onsetTimeInward);
end

if ~isempty(onsetTime)
    onset = onsetTime;
end

onsetPlot = [onset(1)/sampleRate onsetTimeInward/sampleRate; trace(onsetTime) onsetInwardAmplitude];
halfwidthPlot = [rise50/sampleRate decay50/sampleRate; trace(rise50) trace(rise50)];
risePlot = [rise20/sampleRate rise80/sampleRate; trace(rise20) trace(rise80)];
decayPlot = [decay20/sampleRate decay80/sampleRate; trace(decay20) trace(decay80)];

if iscell(saveFig)
   saveFigName = saveFig{1};
elseif saveFig==0
    saveFigName = 'empty';
else
    saveFigName = saveFig;
end
if showFig == 1
    %Plot trace with parameters
        time = (1:length(trace))/sampleRate;
        h = figure;
        hold on
        plot(time, trace, 'k')
        title(saveFigName)
        xlabel('Time (s)')
        ylabel('Current (pA)')

        %plot parameters
        plot(onsetPlot(1,:),onsetPlot(2,:),'b*')
        plot(risePlot(1,:),risePlot(2,:),'go')
        plot(decayPlot(1,:),decayPlot(2,:),'ro')
        plot((amplitude(2)+pulseParam(1))/sampleRate, amplitude(1),'mx')
        plot(halfwidthPlot(1,:),halfwidthPlot(2,:),'cd')

    if saveFig ~= 0
        saveas(h, saveFig,'eps')
    end
end

Aout = paramStruct;
end