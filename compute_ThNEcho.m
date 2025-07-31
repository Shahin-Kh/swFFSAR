function ThNEcho=compute_ThNEcho(data,NstartNoise,NendNoise)

    

%        Function -> compute_ThNEcho(data,NstartNoise,NendNoise)
%                    Function providing the Thermal Noise computed from the waveform
% 
%           Input :
%                 data -> waveform data matrix (dimensions are rangeXrecords)
%                 NstartNoise -> value of the range gate from which to start the noise window (counting from 1)
%                 NendNoise -> value of the range gate at which to stop the noise window (counting from 1)
    

    NstartNoise=floor(NstartNoise);
    NendNoise=floor(NendNoise);

    data = sort(data(1:size(data,1) / 2,:),1);
    data(find(data <= 0)) = nan;

    if  (NstartNoise-1<0)
        NstartNoise=0;
    end

    if  (floor(NendNoise>size(data,1)) / 2)
        NendNoise=floor(size(data,1) / 2);
    end

    ThNEcho=median(data(NstartNoise:NendNoise,:),1);
    ThNEcho( find(isnan(ThNEcho)) )= median(ThNEcho);


    end