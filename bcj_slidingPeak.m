function [Aout] =bcj_slidingPeak(trace, sampleRate, windowSize, peakWindow, XPSCmode)
% Function to determine the peak amplitude using a sliding window of size (windowSize).
% First the trace is calculate based on the sliding window averging, followed by detection of the peak response.
% Difference between peak response and baseline amplitude (first 50 ms[blPeriod] before pulseParam, this can be defined in peakWindow as the third parameter) is used for amplitude
% windowSize : in ms
% peakWindow [start end]: in datapoints
% XPSCmode : scalar, 1 = inward EPSC/IPSC, 2 = outward EPSC/IPSC.

if length(peakWindow) == 3
    blPeriod = peakWindow(3);
else length(peakWindow) == 2
    blPeriod = 50;
end

windowSize = (windowSize/1000)*sampleRate;

slideTrace = movingmean(trace, windowSize,1,1); % uses imported function movingmean (http://www.mathworks.com/matlabcentral/fileexchange/41859-moving-average-function)
slideTraceWindow = slideTrace(peakWindow(1):peakWindow(2));

baselinePeriod = (blPeriod/1000)*sampleRate;

switch XPSCmode
case 1 % inward current
  [amplitude ampPosition] = min(slideTraceWindow);
case 2 % outward current
  [amplitude ampPosition] = max(slideTraceWindow);
end

amplitude = mean(slideTrace(peakWindow(1)-baselinePeriod: peakWindow(1)))-amplitude;

[Aout] = [amplitude ampPosition];

end
