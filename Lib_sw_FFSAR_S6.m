function [original_wf, regenerated_wf, tracker_range, scalefactor, sat_alt]=Lib_sw_FFSAR_S6(input_file_FFSAR,LAT,LON,Ref_Waveform,segmentation)



%%%%% Select Target Area 
rad=50/110000; %% In order to find closest burst I defined a region with radius of 50 meter around desired point
Target_lat_min=LAT-rad;
Target_lat_max=LAT+rad;
Target_lon_min=LON-rad;
Target_lon_max=LON+rad;


%%%% Configuration
zp=1;                                      %Zeropadding factor in range (integer between 1 and 2)
hamming_range=1;                %Hamming windowing in range (yes:1/no:0)
hamming_az=1;                      %Hamming windowing in azimuth (yes:1/no:0)
posting_rate=640;                  %posting rate in Hz related to the FFSAR along-track multi-looking of the single-looks, the spacing between two multi-looks is equal to `satellite_velocity/posting_rate` meters, (floating number between 20 and 17825)
range_ext_factor=1;          %extension factor in range, processing is done on `128 * range_ext_factor * zp` range gates to be truncated into `128 * zp` central range gates at the end of the process (integer between 1 and 2)


%%%%% Constants
cst.fc=13.575e9;
cst.B=320e6;
cst.Fs=395e6;
cst.c=299792458.0;
cst.n_sample_range=256;
cst.n_pulse_burst=64;
cst.GdB=41.9; % antenna gain in dB


    illumination_time=2.0;        %illumination time in seconds, it defines the time of the synthetic aperture (floating number between 0.08 to 2.3)
    cst.T=32e-6;      %% usable pulse length
    cst.alpha=cst.B/cst.T;
    % cst.BRI=7.2409e-3;
    % cst.PRI=1.1314062e-4;
    cst.BRI=1/(139.26);                %%%1018710*12.5e-9;
    cst.PRI=1/(9.195e3);            %%%4488*12.5e-9;
    cst.theta3dB=0.023387;  %1.34 deg
    cst.GdB=42.2; % antenna gain in dB
    cst.n_sample_range=256;



cst.gamma=sin(cst.theta3dB)^2 / (2*log(2));
%cst.abs_ref_track=44;

cst.tau_offset=0;%77.576/cst.B;
% cst.tau_offset=77.576/cst.B;
cst.tracker_phase_shift=0;%2.567;

cst.sig0_bias_ocean=-2.28;
cst.sig0_bias_ocean_lrm=-2.28;

cst.sig0_bias_ocog_sar=11.51;
cst.sig0_bias_ocog_lrm=0.0;


%%%%% by default values 
n_burst_shift=1;
n_sl_output=180;%floor(n_burst_shift*(cst.BRI/cst.PRI));%906;


%%%%% Coefficients of the hamming window
cst.a1=0.5;
cst.a2=0.5;

%%%%% Variables

cst.t=((-cst.n_sample_range/2):floor(cst.n_sample_range/2)-1)/cst.n_sample_range*cst.T;
cst.tau=(floor(-cst.n_sample_range/2):floor(cst.n_sample_range/2)-1)/cst.B;
cst.eta_burst=(floor(-cst.n_pulse_burst/2):floor(cst.n_pulse_burst/2)-1)*cst.PRI;





%%%%% calculating the number of single-looks (sl) and multi-looks (ml), so that ml is a divisor of sl
n_multilook =1; %floor(n_sl_output*posting_rate*cst.PRI);
n_sl_output =180;% floor(n_sl_output/n_multilook)*n_multilook;

%%%%% convert into pulses or bursts
n_burst_block=floor(illumination_time/cst.BRI);

    burst_time=ncread(input_file_FFSAR,'/data_140/ku/time');
    burst_latitude=ncread(input_file_FFSAR,'/data_140/ku/latitude');
    burst_longitude=ncread(input_file_FFSAR,'/data_140/ku/longitude');
    burst_count_cycle=ncread(input_file_FFSAR,'/data_140/ku/burst_counter');
    surf_type=ncread(input_file_FFSAR,'/data_140/ku/surface_classification_flag');



[burst_size,~] = size(burst_time);

test_bbox=(burst_latitude<=Target_lat_max) & ... 
        (burst_latitude>=Target_lat_min) & ...
        (burst_longitude<=Target_lon_max) & ...
        (burst_longitude>=Target_lon_min);

ind_burst_process=find(test_bbox==1);
burst_bar=ind_burst_process(1:n_burst_shift:end);


previous_burst_time=0;
index_write_nc=0;



% Hamming window in azimuth
    if hamming_az
        hamming_window_az=cst.a1-cst.a2*cos(2*pi*(0:cst.n_pulse_burst-1)/cst.n_pulse_burst);
        % hamming_window_az*=np.sqrt(cst.n_pulse_burst/np.sum(hamming_window_az^2))
        hamming_window_az = hamming_window_az/ sqrt(cst.a1^2+0.5*cst.a2^2);
    end
    % Hamming window in range
    hamming_window_range=cst.a1-cst.a2*cos(2*pi*(0:(cst.n_sample_range-1)*range_ext_factor)/cst.n_sample_range*range_ext_factor);
    % hamming_window_range*=np.sqrt(cst.n_sample_range/np.sum(hamming_window_range^2))
    hamming_window_range = hamming_window_range / sqrt(cst.a1^2+0.5*cst.a2^2);
    
    % Toggles antenna pattern compensation in azimuth
    
    posting_rate = n_multilook/(n_sl_output*cst.PRI);

   % Time and frequency vectors
    t=(-floor(cst.n_sample_range*range_ext_factor/2):floor(cst.n_sample_range*range_ext_factor/2)-1) / (cst.n_sample_range*range_ext_factor)*cst.T;       
    eta=(floor(-n_sl_output/2):(floor(n_sl_output/2)-1))*cst.PRI;
  
     Echo=zeros(cst.n_pulse_burst,cst.n_sample_range,n_burst_block);   
     echo=zeros(cst.n_pulse_burst,cst.n_sample_range,n_burst_block);   




    Out.Multilook_FFSAR=zeros(numel(burst_bar)*n_multilook,cst.n_sample_range*zp);
    Out.alt_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Lat_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Lon_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Doppler_ffsar=zeros(numel(burst_bar)*n_multilook,n_burst_block);
    Out.MSC_ffsar=zeros(numel(burst_bar)*n_multilook,cst.n_sample_range*zp);
    Out.PulsePeak_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.RadialVel_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.RMC_Course_ffsar=zeros(numel(burst_bar)*n_multilook,n_burst_block);
    Out.ScaleFactor_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Seperation_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Time_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Tracker_ffsar=zeros(numel(burst_bar)*n_multilook,1);
    Out.Velocity_ffsar=zeros(numel(burst_bar)*n_multilook,1);



    latBox=[Target_lat_min Target_lat_max Target_lat_max Target_lat_min Target_lat_min];
    lonBox=[Target_lon_min Target_lon_min Target_lon_max Target_lon_max Target_lon_min];

  alt=ncread(input_file_FFSAR,'/data_140/ku/altitude',burst_bar(1),1);
  tracker=ncread(input_file_FFSAR,'/data_140/ku/tracker_range_calibrated_diode',burst_bar(1),1);
   
   

for b=1:numel(burst_bar)





    ind=burst_bar(b);
    deltat=burst_time(ind)-previous_burst_time;
	
    % Checking whether the data buffer needs updating to process burst b is available in the data buffer
   % if b > data_buffer.index_end - floor(n_burst_block/2) or np.isnan(data_buffer.index_end):
        index_start=max(0,ind-floor(n_burst_block/2));
        index_max=min(burst_size,ind+floor(n_burst_block/2));
		
            % Reading the data from index_start to index_max

  time=ncread(input_file_FFSAR,'/data_140/ku/time',index_start,index_max-index_start);
  lat=ncread(input_file_FFSAR,'/data_140/ku/latitude',index_start,index_max-index_start);
  lon=ncread(input_file_FFSAR,'/data_140/ku/longitude',index_start,index_max-index_start);

  xs=ncread(input_file_FFSAR,'/data_140/ku/position_vector',[1 index_start],[1 index_max-index_start]);
  ys=ncread(input_file_FFSAR,'/data_140/ku/position_vector',[2 index_start],[1 index_max-index_start]);
  zs=ncread(input_file_FFSAR,'/data_140/ku/position_vector',[3 index_start],[1 index_max-index_start]);

  vxs=ncread(input_file_FFSAR,'/data_140/ku/velocity_vector',[1 index_start],[1 index_max-index_start]);
  vys=ncread(input_file_FFSAR,'/data_140/ku/velocity_vector',[2 index_start],[1 index_max-index_start]);
  vzs=ncread(input_file_FFSAR,'/data_140/ku/velocity_vector',[3 index_start],[1 index_max-index_start]);

  alt=ncread(input_file_FFSAR,'/data_140/ku/altitude',index_start,index_max-index_start);
  orb_alt_rate=ncread(input_file_FFSAR,'/data_140/ku/altitude_rate',index_start,index_max-index_start);

  h0_applied=ncread(input_file_FFSAR,'/data_140/ku/tm_h0',index_start,index_max-index_start);
  cor2_applied=ncread(input_file_FFSAR,'/data_140/ku/tm_cor2',index_start,index_max-index_start);
  tracker=ncread(input_file_FFSAR,'/data_140/ku/tracker_range_calibrated_diode',index_start,index_max-index_start);
%   cog=ncread(input_file_FFSAR,'cog_cor_l1a_echo_sar_ku',index_start,index_max-index_start);
  agc=ncread(input_file_FFSAR,'/data_140/ku/variable_digital_gain',index_start,index_max-index_start);
%   scale_factor=ncread(input_file_FFSAR,'/data_140/ku/iq_scale_factor',index_start,index_max-index_start);
%   sig0_cal=ncread(input_file_FFSAR,'sig0_cal_ku_l1a_echo_sar_ku',index_start,index_max-index_start);
  i_count=ncread(input_file_FFSAR,'/data_140/ku/i_samples',[1 1 index_start],[256 64 index_max-index_start]);
  q_count=ncread(input_file_FFSAR,'/data_140/ku/q_samples',[1 1 index_start],[256 64 index_max-index_start]);
  burst_count_cycle=ncread(input_file_FFSAR,'/data_140/ku/burst_counter',index_start,index_max-index_start);
  burst_power=ncread(input_file_FFSAR,'/global/ku/burst_phase_array_cor',1,64);
  burst_phase=ncread(input_file_FFSAR,'/global/ku/burst_power_array_cor',1,64);

  pri=ncread(input_file_FFSAR,'/data_140/ku/tm_pri',index_start,index_max-index_start)./4;

 for k=1:index_max-index_start
  % % tracker_count computation
    h=h0_applied(k)+(burst_count_cycle(k)-1)*floor((cor2_applied(k)/4)/16);

    % % h=h0_applied(k)*4*(3.125/64*1e-9)+(burst_count_cycle(k)-1)*floor(((cor2_applied(k)*(3.125/1024*1e-9))/4)/16);

    
    hc=floor(h/(2^8))*(2^8);
    hf=h-hc;
    if hf>127
        hc=hc+256;
        hf=hf-256;
    end
    tracker_count(k)=hc/(2^8);

    % applying calibrations
    echo1=i_count(:,:,k)'+1j*q_count(:,:,k)';





    % AGC correction
    % % echo1=echo1*sqrt(10^(agc(k)/10.0));

        % Tracker experimental correction
        % % echo1=echo1*(exp(1j*tracker_count(k)*cst.tracker_phase_shift-2*1j*pi*cst.fc*2/cst.c*tracker(k)));
        % % echo1=echo1*(exp(1j*tracker_count(k)*cst.tracker_phase_shift-2*1j*pi*cst.fc*2/cst.c*tracker(k)));

        % Range FFT

        % RVP and CAL2 correction
            Echo(:,:,k)=fftshift(fft(echo1,256,2),2);
            Echo(:,:,k)=Echo(:,:,k).*(exp(1j*pi*cst.alpha*(cst.tau-cst.tau_offset).^2));
            % Range IFFT
            echo(:,:,k)=ifft(fftshift(Echo(:,:,k),2),256,2);




   

 end

%%%%%%%%%%%%%%%  FFSAR Function


n_burst=index_max-index_start;
CoB=floor(n_burst/2)+1;

LLA=ecef2lla([xs(CoB),ys(CoB),zs(CoB)],'WGS84');
xyz_focus=lla2ecef([LLA(1),LLA(2),LLA(3)-tracker(CoB)]);

elev_focus=alt(CoB)-tracker(CoB);
        
radial_velocity=-(vxs(CoB)*(xyz_focus(:,1)-xs(CoB)) + vys(CoB)*(xyz_focus(:,2)-ys(CoB)) + vzs(CoB)*(xyz_focus(:,3)-zs(CoB)))/tracker(CoB);

    velocity=sqrt(vxs(CoB)^2+vys(CoB)^2+vzs(CoB)^2);
    spacing=1.2*velocity/posting_rate;
        
    rs=sqrt((xs-xyz_focus(1)).^2 + (ys-xyz_focus(2)).^2 + (zs-xyz_focus(3)).^2);
        
    fds=-2/cst.c*cst.fc*(vxs.*(xyz_focus(1)-xs) + vys.*(xyz_focus(2)-ys) + vzs.*(xyz_focus(3)-zs))./rs;
   
        doppler_freqs=(fds - 2*cst.fc/cst.c*orb_alt_rate');

    rmc_coarse=((rs-rs(CoB)-tracker+tracker(CoB) - fds*cst.c/2/cst.alpha)/cst.c*2*cst.B);
      

       % =============================================================================
        % Interpolation of the position and the velocity        
        % =============================================================================
        
        % Time tag of the pulses
        eta_pulse=time-time(CoB)+cst.eta_burst;
        % eta_pulse=np.arange(-n_burst//2,n_burst//2)[:,None]*cst.BRI + cst.eta_burst[None,:]
        
        % Time tag of the bursts
        eta_burst=eta_pulse(:,(cst.n_pulse_burst/2+1));
        
        % Interpolation functions
        fx=(polyfit(eta_burst,xs,4));
        fy=(polyfit(eta_burst,ys,4));
        fz=(polyfit(eta_burst,zs,4));
        fvx=(polyfit(eta_burst,vxs,4));
        fvy=(polyfit(eta_burst,vys,4));
        fvz=(polyfit(eta_burst,vzs,4));
        
        list_ind_burst=1:n_burst;
        n_burst_process=numel(list_ind_burst);

x_pulse=polyval(fx, reshape(eta_pulse',1,[])); 
y_pulse=polyval(fy, reshape(eta_pulse',1,[])); 
z_pulse=polyval(fz, reshape(eta_pulse',1,[])); 

vx_pulse=polyval(fvx, reshape(eta_pulse',1,[])); 
vy_pulse=polyval(fvy, reshape(eta_pulse',1,[])); 
vz_pulse=polyval(fvz, reshape(eta_pulse',1,[])); 

        tracker_pulse=repelem(tracker,cst.n_pulse_burst,1);


% =============================================================================
% Extending the window and range IFFT        
% =============================================================================

% data_ext=zeros(cst.n_pulse_burst,floor(cst.n_sample_range*range_ext_factor/2)+floor(cst.n_sample_range/2),n_burst_process);
% data=zeros(cst.n_pulse_burst,cst.n_sample_range*range_ext_factor,n_burst);
%     
if hamming_az
    data_ext=(Echo.*hamming_window_az');
else
    data_ext=(Echo);
end   

            data=ifft(fftshift(data_ext,2),256,2);

if hamming_range
    data=data.*hamming_window_range;
end



% =============================================================================
% Focusing points        
% =============================================================================
xx=polyval(fx, reshape(eta',1,[]));    
yy=polyval(fy, reshape(eta',1,[]));    
zz=polyval(fz, reshape(eta',1,[]));    

LLA_focous=ecef2lla([xx',yy',zz'],'WGS84'); %%Lat Lon Alt
xyz_n_focus=lla2ecef([LLA_focous(:,1),LLA_focous(:,2),ones(n_sl_output,1)*elev_focus]);



% =============================================================================
% Loop over focusing points        
% =============================================================================

diff_eta=repelem(eta(floor(n_sl_output/n_multilook/2)+1:floor(n_sl_output/n_multilook):end),floor(n_sl_output/n_multilook));
H=zeros(n_burst*cst.n_pulse_burst,cst.n_sample_range*range_ext_factor);
slc=zeros(n_sl_output,cst.n_sample_range*zp);
% multilook=zeros(n_multilook,cst.n_sample_range*zp);
% msc=zeros(n_multilook,cst.n_sample_range*zp);     
        time_multilook=zeros(n_multilook);
        lat_multilook=zeros(n_multilook);
        lon_multilook=zeros(n_multilook);
        tracker_multilook=zeros(n_multilook);
        alt_multilook=zeros(n_multilook);



                % PLRM : the PRLM is part of the code here for the estimation of the thermal noise level 1b output variable used in level 2 retrackers
         h_plrm=zeros(4,cst.n_pulse_burst,cst.n_sample_range);
%         self.echo_foc_plrm=np.empty((4,cst.n_pulse_burst,cst.n_sample_range),dtype=np.complex)
%         self.stack_complex_plrm=np.empty((4,cst.n_pulse_burst,cst.n_sample_range*self.zp),dtype=np.complex)
%         self.stack_plrm=np.empty((4,cst.n_pulse_burst,cst.n_sample_range*self.zp),dtype=np.float)
%         self.multilook_plrm=np.empty(cst.n_sample_range*self.zp,dtype=np.float)
  rmc_plrm=zeros(4,cst.n_pulse_burst);



            data1 = permute(data,[1 3 2]);
            data = reshape(data1,[],256,1);







for j=1:n_sl_output
 

    r_p=sqrt((x_pulse-xyz_n_focus(j,1)).^2 + (y_pulse-xyz_n_focus(j,2)).^2 + (z_pulse-xyz_n_focus(j,3)).^2);
    fd_p=-2/cst.c*cst.fc*(vx_pulse .* (xyz_n_focus(j,1) - x_pulse)...
                          +vy_pulse .* (xyz_n_focus(j,2) - y_pulse)...
                          +vz_pulse .* (xyz_n_focus(j,3) - z_pulse)) ./ r_p;

    corr_doppler=-fd_p*cst.c/cst.alpha/2;
    delay_range=0;
    corr_range=r_p' - LLA_focous(j,3) + elev_focus + corr_doppler' - tracker_pulse + tracker(floor(n_burst/2)+1) - delay_range - (LLA_focous(floor(n_sl_output/2)+1,3)-LLA_focous(j,3)) - diff_eta(j)*radial_velocity;

    H=exp(-2*1j*pi*cst.alpha*t*2/cst.c.*corr_range + 1j*2*pi*2/cst.c*cst.fc*r_p'); 
    tmp=(data.*H);
    % slc = fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);
    % % tmp2=fftshift(fft((tmp),cst.n_sample_range*range_ext_factor*zp),2);
    % % % slc(j,:) = fftshift(fft(nansum(data.*H),cst.n_sample_range*range_ext_factor*zp),2);

% % % %     for mg=1:size(tmp,1)
% % % %         tmp1=fftshift(fft((tmp(mg,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % %         power_matrix(mg,:)=tmp1.*conj(tmp1);
% % % %     end
% % % % power_matrix_total(:,:,j)=power_matrix;
    % 
    % for mg=1:size(tmp,1)
    %     tot_pow=zeros(1,128);
    %     for mh=1:size(tmp,1)
    %         real_cross_term= abs(real(tmp(mg,:)) .* real(tmp(mh,:)));
    %         imag_cross_term= abs(imag(tmp(mg,:)) .* imag(tmp(mh,:)));
    %         tot_pow = tot_pow + (real_cross_term - imag_cross_term);
    %     end
    %     power_matrix_direct(mg,:,j)=tot_pow;
    % 
    % end

% total_power = [];
%     total_power = sum(power_matrix);

% n_signals = size(tmp, 1);
% for i = 1:n_signals
%       % Cross-term contributions
%     for j = i+1:n_signals
%         % Compute real and imaginary cross-terms
%         real_cross_term = 2 * real(power_matrix(i,:)) .* real(power_matrix(j,:));
%         imag_cross_term = 2 * imag(power_matrix(i,:)) .* imag(power_matrix(j,:));
% 
%         total_power = total_power + sum(real_cross_term + imag_cross_term);
%     end
% end







n = 64*180;%%11520;
m = segmentation;
P = zeros(m, size(tmp, 2)); % Initialize P with the appropriate size

for i = 1:m
    start_idx = (i-1)*n/m + 1;
    if i < m
        end_idx = i*n/m;
    else
        end_idx = n;
    end
    P(i,:) = fftshift(fft(nansum(tmp(start_idx:end_idx, :)), cst.n_sample_range*range_ext_factor*zp), 1);
    P(i,:) = P(i,:) .* conj(P(i,:));
    temp1(i,:,j)=P(i,:);
end


iy = 1;


for ix = 1:m
    start_idx = (ix-1)*n/m + 1;
    if ix < m
        end_idx = ix*n/m;
    else
        end_idx = n;
    end
    
    % del = finddelay((P(ix,:) / max(P(ix,:))), (Ref_Waveform / max(Ref_Waveform)));
    % a = circshift((P(ix,:) / max(P(ix,:))), del);
    % corr = corrcoef((Ref_Waveform / max(Ref_Waveform)), a);
    % correlation = corr(1,2);
    % 
    % if (correlation > 0.66 && abs(del) < 2)
    %     misfit=sqrt((1/1)*sum((P(ix,:) / max(P(ix,:)))-Ref_Waveform).^2);
    %     Total_power(((iy-1)*n/m) + 1: (iy*n/m), :) = tmp(start_idx:end_idx, :) .* 1000 * (1/misfit^2);
    %     iy = iy + 1;
    % end

% % for i=-1:1
% %   a = circshift((P(ix,:) / max(P(ix,:))), i);
% %   tmp=corrcoef((Ref_Waveform / max(Ref_Waveform)), a);
% %   cor(i+2) = tmp(1,2);
% %   misfit(i+2)=sqrt((1/1)*sum((a / max(a))-Ref_Waveform).^2);
% % 
% % end
% a = circshift((P(ix,:) / max(P(ix,:))), 0);
  % misfit=sqrt((1/1)*sum((a / max(a))-Ref_Waveform).^2);
    [correlation, lags] = xcorr(P(ix,:) / max(P(ix,:)), Ref_Waveform / max(Ref_Waveform), 'coeff');
    [max_corr, max_idx] = max(correlation); % Maximum correlation value
    lag(ix) = abs(lags(max_idx));
    RR(ix)=max_corr*100;
    misfit(ix)=sqrt(sum((P(ix,:) / max(P(ix,:)))-Ref_Waveform).^2);
    % tmp1=corrcoef((Ref_Waveform / max(Ref_Waveform)), P(ix,:) / max(P(ix,:)));
    R1(ix)=correlation(lags==0)*100;
    if(lag(ix)<3)
        R1(ix)=RR(ix);
        lag(ix)=0;
    end

    w1=1;%(R1(ix)/100)^1.6;  
w1=1.01^R1(ix);    w1=(w1-1)/((1.01^100)-1); 
w2=(1/(min(misfit(ix))));   w2= (w2-(1/128))/(10^6-(1/128));
w3=(1/(1+(1.4^abs(lag(ix)))));  w3= (w3-(1/(1+(1.4^64))))/(0.5-(1/(1+(1.4^64))));

wgt(ix)= w1  * 1 * w3;


% del = finddelay((P(ix,:) / max(P(ix,:))), (Ref_Waveform / max(Ref_Waveform)));
    %         wgt(ix)= 2^R1(ix)  * (1/(min(misfit(ix)))) * (1/(1+2^abs(lag(ix))));
    %         % wgt(ix)=R1 * RR * (1/(min(misfit))) * (1/(1+2^abs(lag)));
    % 
    % if(lag(ix)<4)
    %     wgt(ix)=1000000000;
    % end


            Total_power(((iy-1)*n/m) + 1: (iy*n/m), :) = tmp(start_idx:end_idx, :) .* wgt(ix);

        % % if ( abs(del) < 2)
        % misfit=sqrt((1/1)*sum((P(ix,:) / max(P(ix,:)))-Ref_Waveform).^2);
        % % % if (del~=0)
        % % %     Total_power(((iy-1)*n/m) + 1: (iy*n/m), :) = tmp(start_idx:end_idx, :) .* 1 * (1/(1*abs(del)*min(misfit)));
        % % % else
        % % %     Total_power(((iy-1)*n/m) + 1: (iy*n/m), :) = tmp(start_idx:end_idx, :) .* 1 * (1/(min(misfit)));
        % % % end
        iy = iy + 1;
        % % end


end
Total_power=Total_power/(sum(wgt));

% % for ix = 1:m
% %     start_idx = (ix-1)*n/m + 1;
% %     if ix < m
% %         end_idx = ix*n/m;
% %     else
% %         end_idx = n;
% %     end
% % 
% %     del = finddelay((P(ix,:) / max(P(ix,:))), (Ref_Waveform / max(Ref_Waveform)));
% %     a = circshift((P(ix,:) / max(P(ix,:))), del);
% %     corr = corrcoef((Ref_Waveform / max(Ref_Waveform)), a);
% %     correlation = corr(1,2);
% % 
% %     if (correlation > 0.66 && abs(del) < 2)
% %         Total_power(((iy-1)*n/m) + 1: (iy*n/m), :) = tmp(start_idx:end_idx, :);
% %         iy = iy + 1;
% %     end
% % end


% % % P(1,:)=fftshift(fft(nansum(tmp(1:11520/5-1,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % P(1,:)=P(1,:).*conj(P(1,:));
% % % 
% % % P(2,:)=fftshift(fft(nansum(tmp(11520/5:2*11520/5-1,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % P(2,:)=P(2,:).*conj(P(2,:));
% % % 
% % % P(3,:)=fftshift(fft(nansum(tmp(2*11520/5:3*11520/5-1,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % P(3,:)=P(3,:).*conj(P(3,:));
% % % 
% % % P(4,:)=fftshift(fft(nansum(tmp(3*11520/5:4*11520/5-1,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % P(4,:)=P(4,:).*conj(P(4,:));
% % % 
% % % P(5,:)=fftshift(fft(nansum(tmp(4*11520/5:11520,:)),cst.n_sample_range*range_ext_factor*zp),2);
% % % P(5,:)=P(5,:).*conj(P(5,:));

% PT=fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);

% % accepted_parts=[];
% % % iy=1;
% % % Total_power=[];
% % % for ix=1:5
% % %     del=finddelay((P(ix,:)/max(P(ix,:))),(Ref_Waveform/max(Ref_Waveform)));
% % %     a=circshift((P(ix,:)/max(P(ix,:))),del);
% % %     corr=corrcoef((Ref_Waveform/max(Ref_Waveform)),a);
% % %     correlation=corr(1,2);
% % %         if(correlation>0.66 && abs(del)<2)
% % %             Total_power((((iy-1)*11520)/5)+1:((iy*11520)/5),:)=tmp((((ix-1)*11520)/5)+1:((ix*11520)/5),:);
% % %             iy=iy+1;
% % %         end
% % % end
% % % 
if(~isempty(Total_power))
    slc(j,:)=fftshift(fft(nansum(Total_power),cst.n_sample_range*range_ext_factor*zp),1);
else
    % slc(j,:)=fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);
    slc(j,:)=0;
end


    slc1(j,:)=fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),1);


end




        slc=slc*10^(-agc(floor(n_burst/2)+1)/20) / sqrt(cst.n_pulse_burst) / (cst.n_sample_range);
        aa=squeeze(reshape(real(slc.*conj(slc)),n_multilook,[],cst.n_sample_range*zp));
        aa( ~any(aa,2), : ) = [];  %rows
       % Scaling the power
        PT=mean(aa)'*10^(-agc(floor(n_burst/2)+1)/20) / sqrt(cst.n_pulse_burst) / (cst.n_sample_range);

        slc1=slc1*10^(-agc(floor(n_burst/2)+1)/20) / sqrt(cst.n_pulse_burst) / (cst.n_sample_range);
        bb=squeeze(reshape(real(slc1.*conj(slc1)),n_multilook,[],cst.n_sample_range*zp));
        bb( ~any(bb,2), : ) = [];  %rows
       % Scaling the power
        PT1=mean(bb)'*10^(-agc(floor(n_burst/2)+1)/20) / sqrt(cst.n_pulse_burst) / (cst.n_sample_range);



    % aa=squeeze(reshape((real(slc.*conj(slc))),n_multilook,[],cst.n_sample_range*zp));

% % for ix=1:180
% %     del(ix)=finddelay((aa(ix,:)/max(aa(ix,:))),(Ref_Waveform/max(Ref_Waveform)));
% %     a=circshift((aa(ix,:)/max(aa(ix,:))),del);
% %     corr=corrcoef((Ref_Waveform/max(Ref_Waveform)),a);
% %     correlation(ix)=corr(1,2);
% %      if(correlation(ix)>0.66 && abs(del(ix))<3)
% %         weights(ix)=1;
% %      else
% %         weights(ix)=0;
% %      end
% % end




      


for ix=1:n_sl_output


    [correlation, lags] = xcorr(aa(ix,:) / max(aa(ix,:)), Ref_Waveform / max(Ref_Waveform), 'coeff');
    [max_corr, max_idx] = max(correlation); % Maximum correlation value
    lag(ix) = abs(lags(max_idx));
    RR(ix)=max_corr*100;
    misfit(ix)=sqrt(sum((aa(ix,:) / max(aa(ix,:)))-Ref_Waveform).^2);
    % tmp1=corrcoef((Ref_Waveform / max(Ref_Waveform)), aa(ix,:) / max(aa(ix,:)));
    R1(ix)=correlation(lags==0)*100;



% % %     if(lag(ix)<3)
% % %         R1(ix)=RR(ix);
% % %         lag(ix)=0;
% % %     end
% % %     w1=(R1(ix)/100)^1.4;  
% % % % w1=1.01^R1(ix);    w1=(w1-1)/((1.01^100)-1); 
% % % w2=(1/(min(misfit(ix))));   w2= (w2-(1/128))/(10^6-(1/128));
% % % w3=(1/(1+(1.4^abs(lag(ix)))));  w3= (w3-(1/(1+(1.4^64))))/(0.5-(1/(1+(1.4^64))));
% % % wgt(ix)= w1  * 1 * w3;
% % % 




if(lag(ix)<3)
        R1(ix)=RR(ix);
        lag(ix)=0;
end

    w1=(R1(ix)/100)^1.4;  
% w1=1.01^R1(ix);    w1=(w1-1)/((1.01^100)-1); 
w2=(1/(min(misfit(ix))));   w2= (w2-(1/128))/(10^6-(1/128));
w3=(1/(1+(1.4^abs(lag(ix)))));  w3= (w3-(1/(1+(1.4^64))))/(0.5-(1/(1+(1.4^64))));



wgt(ix)= w1  * 1 * w3;


end
% tmp = (1 ./ (abs(del)'.*misfit2')) .* aa;
tmp = wgt' .* aa;

% tmp=(1/(misfit2')^2).*aa;
% tmp( ~any(tmp,2), : ) = [];  %rows
PT2=sum(tmp,1)/(sum(wgt));







for ix=1:n_sl_output


    [correlation, lags] = xcorr(bb(ix,:) / max(bb(ix,:)), Ref_Waveform / max(Ref_Waveform), 'coeff');
    [max_corr, max_idx] = max(correlation); % Maximum correlation value
    lag(ix) = abs(lags(max_idx));
    RR(ix)=max_corr*100;
    misfit(ix)=sqrt(sum((bb(ix,:) / max(bb(ix,:)))-Ref_Waveform).^2);
    % tmp1=corrcoef((Ref_Waveform / max(Ref_Waveform)), bb(ix,:) / max(bb(ix,:)));
    R1(ix)=correlation(lags==0)*100;
    RR(ix)=correlation(max_idx)*100;
    if(lag<2)
        R1(ix)=RR(ix);
    end

            wgt(ix)= R1(ix)  * (1/(min(misfit(ix)))) * (1/(1+2^abs(lag(ix))));
            % wgt(ix)=R1 * RR * (1/(min(misfit))) * (1/(1+2^abs(lag)));


end
% tmp = (1 ./ (abs(del)'.*misfit2')) .* aa;
tmp = wgt' .* bb;

% tmp=(1/(misfit2')^2).*aa;
% tmp( ~any(tmp,2), : ) = [];  %rows
PT3=sum(tmp,1)/(sum(wgt));







        % Time tag of multilooks
        time_multilook=time(floor(n_burst/2)+1)+eta(floor(n_sl_output/n_multilook/2)+1:floor(n_sl_output/n_multilook):end);

        % Interplation of the latitude and longitude of the multilooks
        lat_multilook=interp1(eta_burst,lat,time_multilook-time(floor(n_burst/2)+1),'spline');
        lon_multilook=interp1(eta_burst,lon,time_multilook-time(floor(n_burst/2)+1),'spline');

        % Adding the COG to the tracker
        % alt_multilook[:]=data_block.alt[n_burst//2]
        alt_multilook=alt(floor(n_burst/2)+1) + radial_velocity*eta(floor(n_sl_output/n_multilook/2):floor(n_sl_output/n_multilook):end);




     % Scale factor
        % scalefactor=agc(floor(n_burst/2)+1) + sig0_cal(floor(n_burst/2)+1) + 30*log10(alt(floor(n_burst/2)+1)) + 10*log10(1+alt(floor(n_burst/2)+1)/6378137.0) + 10*log10(4*16*pi^2*cst.fc^2/cst.c^3/256/cst.Fs) -2*cst.GdB;
        
        % Pulse peakiness
        pulse_peakiness=max(PT,[],2)./mean(PT,2) * 85/128;


original_wf=PT1';
sat_alt=alt(CoB);
tracker_range=tracker(CoB);%+cog(CoB);
regenerated_wf=PT2;
scalefactor=0;


end









