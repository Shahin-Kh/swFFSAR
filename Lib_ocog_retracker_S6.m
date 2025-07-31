function [WL] = Lib_ocog_retracker_S6(waveform, tracker_range, sat_alt, corrections)

percentage=0.5;
skip=0;
ref_bin=50;
FS=395e6;
SOL=299792458;

waveform=waveform/max(waveform);

waveform = waveform .^ 2;
sq_sum = sum(waveform);
waveform = waveform .^ 4;
qa_sum = sum(waveform);

% Compute OCOG threshold
threshold = percentage * sqrt(qa_sum / sq_sum);
ind_first_over = find(waveform > threshold, 1, 'first');

if isempty(ind_first_over) || ind_first_over == 1
    WL = NaN;
    return;
end

% Linear interpolation for sub-bin accuracy
decimal = (waveform(ind_first_over-1) - threshold) / ...
          (waveform(ind_first_over-1) - waveform(ind_first_over));
      
range_bin = skip + ind_first_over - 1 + decimal;


WL=sat_alt-tracker_range-(((range_bin*1)-128)*(SOL/2/FS))-corrections;


