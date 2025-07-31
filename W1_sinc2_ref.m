function ref_wf=W1_sinc2_ref(wf,P_loc,N)

normalized_wf=wf/max(wf);
idx1 = find(normalized_wf(1:P_loc) < 0.5, 1, 'last');
idx2 = find(normalized_wf(P_loc:end) < 0.5, 1, 'first')+P_loc;
if(isempty(idx1))
    idx1=1;
end
if(isempty(idx2))
    idx1=1;
end

alpha= 1 / (1*( idx2-idx1));

% Generate a vector from -N/2 to N/2
x = linspace(-N/2, N/2, N-1);
% x = linspace(0, N, N);

% Compute the normalized sinc squared waveform
% alpha = 0.02;  % Try different values like 0.7, 0.3, etc.

% Compute widened sinc^2 function
% sinc_sq_waveform = (sin(pi * alpha * x) ./ (pi * alpha * x)).^2;
 sinc_sq_waveform = sinc(alpha*x).^2;

% Normalize the waveform
sinc_sq_waveform = sinc_sq_waveform / max(sinc_sq_waveform);
ref_wf=circshift(sinc_sq_waveform,round(P_loc-(N/2)));
ref_wf(N)=0;
% Plot the waveform
% % figure;
% % plot(1:127, ref_wf);
% % title('Normalized Sinc^2 Waveform');
% % xlabel('Bins');
% % ylabel('Amplitude');

end

