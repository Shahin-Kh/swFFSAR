function ref_wf=Lib_sinc2_ref(P_loc,N)

% Number of bins
% N = 128;

% Generate a vector from -N/2 to N/2
x = linspace(-N/2, N/2, N);

% Compute the normalized sinc squared waveform
sinc_sq_waveform = sinc(x).^2;

% Normalize the waveform
sinc_sq_waveform = sinc_sq_waveform / sum(sinc_sq_waveform);
ref_wf=circshift(sinc_sq_waveform,P_loc-(N/2));
% Plot the waveform
% % % figure;
% % % plot(1:128, ref_wf);
% % % title('Normalized Sinc^2 Waveform');
% % % xlabel('Bins');
% % % ylabel('Amplitude');

end

