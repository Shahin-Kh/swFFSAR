function [original_wf, regenerated_wf, tracker_range, scalefactor, sat_alt]=Lib_sw_FFSAR(input_file_FFSAR,LAT,LON,Ref_Waveform,segmentation)

%%%%% by default values 
n_sl_output=180;
n_multilook=1;
%%%

%%%% Configuration
illumination_time=2.3;            %illumination time in seconds, it defines the time of the synthetic aperture (floating number between 0.08 to 2.3)
n_burst_block=floor(illumination_time*78.5307);     %% 180 bursts for sentinel 3 over each point on ground  % floor(illumination_time*cst.BRF)

%%%% Hamming window
hamming_az=1;                     %Hamming windowing in azimuth (yes:1/no:0)
zp=1;                             %Zeropadding factor in range (integer between 1 and 2)

hamming_range=1;                  %Hamming windowing in range (yes:1/no:0)
range_ext_factor=1;               %extension factor in range, processing is done on `128 * range_ext_factor * zp` range gates to be truncated into `128 * zp` central range gates at the end of the process (integer between 1 and 2)

cst.a1=0.5;         % Coefficients of the hamming window
cst.a2=0.5;         % Coefficients of the hamming window
%%%

%%%%% Constants S3

cst.c=299792458.0;
cst.n_sample_range=128;
cst.n_pulse_burst=64;
cst.fc=13.575e9;
cst.B=320e6;
cst.Fs=320e6;
cst.T=48.95e-6;  %%Pulse Duration > for CryoSat:44.8e-6
cst.alpha=cst.B/cst.T;
cst.BRI=1018710*12.5e-9;    %   78.74 hz
cst.PRI=4488*12.5e-9;   % 17.825 Mhz
cst.GdB=41.9; % antenna gain in dB

cst.tau_offset=77.576/cst.B;
cst.tracker_phase_shift=2.567;  %   0.82 * pi      

%%%%% Variables
cst.t = (-cst.n_sample_range/2:1:cst.n_sample_range/2-1) * (cst.T/cst.n_sample_range);
cst.tau = (-cst.n_sample_range/2:1:cst.n_sample_range/2-1) * (1/cst.B);
cst.eta_burst = (-cst.n_pulse_burst/2:1:cst.n_pulse_burst/2-1) * cst.PRI;

% Hamming window in azimuth
    if hamming_az
        hamming_window_az=cst.a1-cst.a2*cos(2*pi*(0:cst.n_pulse_burst-1)/cst.n_pulse_burst);
        % hamming_window_az*=np.sqrt(cst.n_pulse_burst/np.sum(hamming_window_az^2))
        hamming_window_az = hamming_window_az/ sqrt(cst.a1^2+0.5*cst.a2^2);
    end
% Hamming window in range
    if hamming_range
        hamming_window_range=cst.a1-cst.a2*cos(2*pi*(0:(cst.n_sample_range-1)*range_ext_factor)/cst.n_sample_range*range_ext_factor);
        % hamming_window_range*=np.sqrt(cst.n_sample_range/np.sum(hamming_window_range^2))
        hamming_window_range = hamming_window_range / sqrt(cst.a1^2+0.5*cst.a2^2);
    end
    % Toggles antenna pattern compensation in azimuth
    
    posting_rate = n_multilook/(n_sl_output*cst.PRI);

   % Time and frequency vectors
    t=(-floor(cst.n_sample_range*range_ext_factor/2):floor(cst.n_sample_range*range_ext_factor/2)-1) / (cst.n_sample_range*range_ext_factor)*cst.T;       
    eta=(floor(-n_sl_output/2):(floor(n_sl_output/2)-1))*cst.PRI;

    posting_rate = n_multilook/(n_sl_output*cst.PRI);

rad=50/110000; %% In order to find closest burst I defined a region with radius of 50 meter around desired point
Target_lat_min=LAT-rad;
Target_lat_max=LAT+rad;
Target_lon_min=LON-rad;
Target_lon_max=LON+rad;

%%% Extracting the coordination of all of the bursts
burst_latitude=ncread(input_file_FFSAR,'lat_l1a_echo_sar_ku');
burst_longitude=ncread(input_file_FFSAR,'lon_l1a_echo_sar_ku');

%%% Finding the nearest burst 
test_bbox=(burst_latitude<=Target_lat_max) & ... 
        (burst_latitude>=Target_lat_min) & ...
        (burst_longitude<=Target_lon_max) & ...
        (burst_longitude>=Target_lon_min);

ind_burst_process=find(test_bbox==1);
burst_bar=ind_burst_process(1:end);

ind=burst_bar(1);
burst_size=numel(burst_latitude);
% Checking whether the data buffer needs updating to process burst b is available in the data buffer
% if b > data_buffer.index_end - floor(n_burst_block/2) or np.isnan(data_buffer.index_end):
    index_start=max(0,ind-floor(n_burst_block/2));
    index_max=min(burst_size,ind+floor(n_burst_block/2));


%%%%% CAL2
%  CAL2=[0.34307632   0.36192444   0.3939476    0.43641937   0.48737532          0.54574859   0.60562903   0.66930281   0.72598263   0.78652351          0.83173867   0.87549745   0.91067044   0.93888334   0.95649523          0.96975496   0.97447969   0.972329     0.9655552    0.96352306          0.95490339   0.94523378   0.93588593   0.92487849   0.91724103          0.91414201   0.91424362   0.91209294   0.92015377   0.92382855          0.93051769   0.94084774   0.95011092   0.95884913   0.96679142          0.9708049    0.97954311   0.98064385   0.98257438   0.98155831          0.97903507   0.97190564   0.96975496   0.96218523   0.95202452          0.9478925    0.93866319   0.92989111   0.92521719   0.92144079          0.91697008   0.92105129   0.92393016   0.92794364   0.93358284          0.94050905   0.94873923   0.95805321   0.96772282   0.97290478          0.97803594   0.98355659   0.9852331    0.99166822   0.98591048          0.98909417   0.98145671   0.97791739   0.96763814   0.96186347          0.95163503   0.94660548   0.93801968   0.93180471   0.927334            0.92667355   0.92299876   0.93163537   0.930179     0.94223637          0.94570795   0.9569186    0.95952651   0.97066942   0.97617314          0.98448798   0.99385277   0.993277     0.9956817    0.99561396          0.99276896   0.98785795   0.97986486   0.97138067   0.96582615          0.95492032   0.94892551   0.94252426   0.93407394   0.9336167           0.93300706   0.93674959   0.94062759   0.94819732   0.95439535          0.96536892   0.97691826   0.98875548   0.99695179   1.0                  0.9999492    0.99505512   0.98447105   0.96518264   0.93387072          0.90130565   0.85976529   0.80376285   0.75106264   0.68832027          0.62268209   0.55982117   0.50131242   0.44847674   0.4024826           0.36854583   0.34620921   0.33818225];
%  gprw_mean = circshift(CAL2,1);
% This Vector of callibration is fixed over entire orbit. The reason why I
% used try and catch is because after January 2018 Eumetsat made a change
% in format of this vector in L1A data products to avoid repitation 
try
    gprw_mean=squeeze(ncread(input_file_FFSAR,'gprw_meas_ku_l1a_echo_sar_ku',[1 1 1],[128 1 1]))';
catch
    gprw_mean=squeeze(ncread(input_file_FFSAR,'gprw_meas_ku_l1a_echo_sar_ku',1,128))';
end


  agc=ncread(input_file_FFSAR,'agc_ku_l1a_echo_sar_ku',index_start,index_max-index_start);
  cog=ncread(input_file_FFSAR,'cog_cor_l1a_echo_sar_ku',index_start,index_max-index_start);


  time=ncread(input_file_FFSAR,'time_l1a_echo_sar_ku',index_start,index_max-index_start);
  lat=ncread(input_file_FFSAR,'lat_l1a_echo_sar_ku',index_start,index_max-index_start);
  lon=ncread(input_file_FFSAR,'lon_l1a_echo_sar_ku',index_start,index_max-index_start);
  xs=ncread(input_file_FFSAR,'x_pos_l1a_echo_sar_ku',index_start,index_max-index_start);
  ys=ncread(input_file_FFSAR,'y_pos_l1a_echo_sar_ku',index_start,index_max-index_start);
  zs=ncread(input_file_FFSAR,'z_pos_l1a_echo_sar_ku',index_start,index_max-index_start);
  vxs=ncread(input_file_FFSAR,'x_vel_l1a_echo_sar_ku',index_start,index_max-index_start);
  vys=ncread(input_file_FFSAR,'y_vel_l1a_echo_sar_ku',index_start,index_max-index_start);
  vzs=ncread(input_file_FFSAR,'z_vel_l1a_echo_sar_ku',index_start,index_max-index_start);
  alt=ncread(input_file_FFSAR,'alt_l1a_echo_sar_ku',index_start,index_max-index_start);
  orb_alt_rate=ncread(input_file_FFSAR,'orb_alt_rate_l1a_echo_sar_ku',index_start,index_max-index_start);

  h0_applied=ncread(input_file_FFSAR,'h0_applied_l1a_echo_sar_ku',index_start,index_max-index_start);
  cor2_applied=ncread(input_file_FFSAR,'cor2_applied_l1a_echo_sar_ku',index_start,index_max-index_start);
  tracker=ncread(input_file_FFSAR,'range_ku_l1a_echo_sar_ku',index_start,index_max-index_start);
  sig0_cal=ncread(input_file_FFSAR,'sig0_cal_ku_l1a_echo_sar_ku',index_start,index_max-index_start);
  
  i_count=ncread(input_file_FFSAR,'i_meas_ku_l1a_echo_sar_ku',[1 1 index_start],[128 64 index_max-index_start]);
  q_count=ncread(input_file_FFSAR,'q_meas_ku_l1a_echo_sar_ku',[1 1 index_start],[128 64 index_max-index_start]);

  burst_count_cycle=ncread(input_file_FFSAR,'burst_count_cycle_l1a_echo_sar_ku',index_start,index_max-index_start);

try
  burst_power=ncread(input_file_FFSAR,'burst_power_cor_ku_l1a_echo_sar_ku',1,64);
  burst_phase=ncread(input_file_FFSAR,'burst_phase_cor_ku_l1a_echo_sar_ku',1,64);
catch 
  burst_power=ncread(input_file_FFSAR,'burst_power_cor_ku_l1a_echo_sar_ku',[1 index_start],[64 1]);
  burst_phase=ncread(input_file_FFSAR,'burst_phase_cor_ku_l1a_echo_sar_ku',[1 index_start],[64 1]);
end
  %%%%%     Calibrate one burst and compute the tracker count and frequency domain pulses. 

% % test=ncread(input_file_FFSAR,'cor2_applied_l1a_echo_sar_ku',index_start,index_max-index_start);
% % figure;plot(test)
% % 
% % test=ncread(input_file_FFSAR,'scale_factor_c_l1a_echo_sar_ku',index_start,index_max-index_start);
% % figure;plot(test)
% % 


Echo=zeros(cst.n_pulse_burst,cst.n_sample_range,n_burst_block);   
echo=zeros(cst.n_pulse_burst,cst.n_sample_range,n_burst_block);   

for k=1:180     % number of single looks = 180
  % tracker_count computation
    h=h0_applied(k)+(burst_count_cycle(k)-1)*floor((cor2_applied(k)/4)/16);
    hc=floor(h/(2^8))*(2^8);
    hf=h-hc;
    if hf>127
        hc=hc+256;
        hf=hf-256;
    end
    tracker_count(k)=hc/(2^8);

    % applying calibrations
    echo1=i_count(:,:,k)'+1j*q_count(:,:,k)';      


    % pi shift between bursts       exp(1j*pi)=-1       < for Sentinel 3
    if (mod(burst_count_cycle(k) , 2) == 1)
        echo1=-1*echo1;
    end



    % AGC correction
    echo1=echo1*sqrt(10^(agc(k)/10.0));

    % Tracker experimental correction
    echo1=echo1*(exp(1j*tracker_count(k)*cst.tracker_phase_shift-2*1j*pi*cst.fc*2/cst.c*tracker(k)));
    
    % Range FFT
    Echo(:,:,k)=fftshift(fft(echo1,128,2),2);
% % Residual Video Phase (RVP) correction, the Residual Range Phase (RRP) correction:
% % The RVP correction needs not to be applied to S6 pulse data, as the RVP is a residual phase of the deramp-on-receive radars CS2 and S3 only.
    % RVP and CAL2 correction
    Echo(:,:,k)=Echo(:,:,k).*(exp(1j*pi*cst.alpha*(cst.tau-cst.tau_offset).^2) ./ sqrt(gprw_mean));
    
    % Range IFFT
    echo(:,:,k)=ifft(fftshift(Echo(:,:,k),2),128,2);
  end

 %%%%%%%%%%%%%%%  FFSAR Function

    n_burst=180;                % number of bursts  180
    CoB=floor(n_burst/2)+1;     % Central Burst index   91
    
    LLA=ecef2lla([xs(CoB),ys(CoB),zs(CoB)],'WGS84');                    % Location of Satellite at central Burst
    xyz_focus=lla2ecef([LLA(1),LLA(2),LLA(3)-tracker(CoB)]);            % ECEF coordination of focusing point
    
    elev_focus=alt(CoB)-tracker(CoB);                                   % approximate altitute of focusing point
            
    radial_velocity=-(vxs(CoB)*(xyz_focus(:,1)-xs(CoB)) + vys(CoB)*(xyz_focus(:,2)-ys(CoB)) + vzs(CoB)*(xyz_focus(:,3)-zs(CoB)))/tracker(CoB);      %% radial_velocity of satellite

    velocity=sqrt(vxs(CoB)^2+vys(CoB)^2+vzs(CoB)^2);
    spacing=1.2*velocity/posting_rate;                  %% measurements displacements   
        
    rs=sqrt((xs-xyz_focus(1)).^2 + (ys-xyz_focus(2)).^2 + (zs-xyz_focus(3)).^2);    %% vector from focusing point to center of each burst
        
    fds=-2/cst.c*cst.fc*(vxs.*(xyz_focus(1)-xs) + vys.*(xyz_focus(2)-ys) + vzs.*(xyz_focus(3)-zs))./rs;
        
    doppler_freqs=(fds - 2*cst.fc/cst.c*orb_alt_rate);      %% Doppler value for each burst

     
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

    x_pulse=polyval(fx, reshape(eta_pulse',1,[])); 
    y_pulse=polyval(fy, reshape(eta_pulse',1,[])); 
    z_pulse=polyval(fz, reshape(eta_pulse',1,[])); 
    
    vx_pulse=polyval(fvx, reshape(eta_pulse',1,[])); 
    vy_pulse=polyval(fvy, reshape(eta_pulse',1,[])); 
    vz_pulse=polyval(fvz, reshape(eta_pulse',1,[])); 
 %% Here We've just calculated the position and velocity of satellite at the spot that pulse is generated

        tracker_pulse=repelem(tracker,cst.n_pulse_burst,1);



    if hamming_az
        data_ext=(Echo.*hamming_window_az');
    else
        data_ext=Echo;
    end   
data=ifft(fftshift(data_ext,2),128,2);

    if hamming_range
        data=data.*hamming_window_range;
    end


% =============================================================================
% Focusing points        
% =============================================================================
xx=polyval(fx, reshape(eta',1,[]));    
yy=polyval(fy, reshape(eta',1,[]));    
zz=polyval(fz, reshape(eta',1,[]));    
% % % % % % % % % % % % % % % % % % % % % % dis=sqrt((xx(180)-xx(1)).^2+(yy(180)-yy(1)).^2+(zz(180)-zz(1)).^2)
LLA_focous=ecef2lla([xx',yy',zz'],'WGS84'); %%Lat Lon Alt
xyz_n_focus=lla2ecef([LLA_focous(:,1),LLA_focous(:,2),ones(n_sl_output,1)*elev_focus]);


diff_eta=repelem(eta(floor(n_sl_output/n_multilook/2)+1:floor(n_sl_output/n_multilook):end),floor(n_sl_output/n_multilook));
H=zeros(n_burst*cst.n_pulse_burst,cst.n_sample_range*range_ext_factor);
slc=zeros(n_sl_output,cst.n_sample_range*zp);
slc1=zeros(n_sl_output,cst.n_sample_range*zp);

data1 = permute(data,[1 3 2]);
data = reshape(data1,[],128,1);


RIP=0;

j=(n_sl_output/2);
DDM_Ave=zeros(180,128);
%  for j=1:n_sl_output

r_p=sqrt((x_pulse-xyz_n_focus(j,1)).^2 + (y_pulse-xyz_n_focus(j,2)).^2 + (z_pulse-xyz_n_focus(j,3)).^2);
            fd_p=-2/cst.c*cst.fc*(vx_pulse .* (xyz_n_focus(j,1) - x_pulse)...
                                  +vy_pulse .* (xyz_n_focus(j,2) - y_pulse)...
                                  +vz_pulse .* (xyz_n_focus(j,3) - z_pulse)) ./ r_p;
    
            corr_doppler=-fd_p*cst.c/cst.alpha/2;
            delay_range=0;%-0.5;
            corr_range=r_p' - LLA_focous(j,3) + elev_focus + corr_doppler' - tracker_pulse + tracker(floor(n_burst/2)+1) - delay_range - (LLA_focous(floor(n_sl_output/2)+1,3)-LLA_focous(j,3)) - diff_eta(j)*radial_velocity;

            t_comp_h=t;
            corr_range_comp_h=corr_range;
            r_p_comp_h=r_p;
            fc=cst.fc;
            alpha=cst.alpha;
            c=cst.c;

    H=exp(-2*1j*pi*alpha*t_comp_h*2/c.*corr_range_comp_h + 1j*2*pi*2/c*fc*r_p_comp_h');

    aa1= permute(data_ext,[1 3 2]);
    aa2 = reshape(aa1,[],128,1);
    WW=(aa2.*H);
a1=reshape(WW,64,180,128);
DDM=squeeze(sum(a1,1));
% imagesc(real(a2))
% imagesc((real(a2).^2+imag(a2).^2),[1 10e8])

 DDM_Ave=DDM_Ave+sqrt(real(DDM).^2+imag(DDM).^2);
%  end
% % imagesc(DDM_Ave,[0 quantile(DDM_Ave(:),.95)])
% % mesh(DDM_Ave)
% % view(18,41)
% % % % % % % % % % % % % % % % % % % % % % dis=sqrt((xyz_n_focus(180,1)-xyz_n_focus(1,1)).^2+(xyz_n_focus(180,2)-xyz_n_focus(1,2)).^2+(xyz_n_focus(180,3)-xyz_n_focus(1,3)).^2)



% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % [m, n] = size(DDM_Ave);
% % % % % % 
% % % % % % % Step 2: Perform k-means Clustering
% % % % % % % Number of clusters
% % % % % % numClusters = 3;
% % % % % % 
% % % % % % % Perform k-means clustering
% % % % % % [idx, C] = kmeans(DDM_Ave, numClusters, 'Replicates', 5);
% % % % % % 
% % % % % % % Reshape the clustered data back into the original 3D image dimensions
% % % % % % classifiedImage = reshape(idx, m, n);
% % % % % % 
% % % % % % % Step 3: Display the Classified Image
% % % % % % % Display the classified image
% % % % % % figure;
% % % % % % imshow(classifiedImage, []);
% % % % % % title('Classified Image');
% % % % % % 
% % % % % % % Step 4: Overlay the Boundaries
% % % % % % % Find boundaries of the classified regions
% % % % % % boundaries = edge(DDM_Ave, 'Canny');
% % % % % % 
% % % % % % % Overlay boundaries on the original image
% % % % % % % Here, I'm assuming the 3D image is a color image (RGB)
% % % % % % originalImage = mat2gray(DDM_Ave); % Convert to grayscale for display
% % % % % % 
% % % % % % % Display the original image with boundaries
% % % % % % figure;
% % % % % % imshow(originalImage);
% % % % % % hold on;
% % % % % % % Overlay the boundaries in red color
% % % % % % visboundaries(boundaries, 'Color', 'r');
% % % % % % title('Classified Image with Boundaries');
% % % % % % hold off;
% % % % % % 
% % % % % % 
% % % % % % 








power_matrix_total=[];
for j=1:180




    r_p=sqrt((x_pulse-xyz_n_focus(j,1)).^2 + (y_pulse-xyz_n_focus(j,2)).^2 + (z_pulse-xyz_n_focus(j,3)).^2);
    fd_p=-2/cst.c*cst.fc*(vx_pulse .* (xyz_n_focus(j,1) - x_pulse)...
                          +vy_pulse .* (xyz_n_focus(j,2) - y_pulse)...
                          +vz_pulse .* (xyz_n_focus(j,3) - z_pulse)) ./ r_p;

    corr_doppler=-fd_p*cst.c/cst.alpha/2;
    delay_range=0;%-0.5;
    corr_range=r_p' - LLA_focous(j,3) + elev_focus + corr_doppler' - tracker_pulse + tracker(floor(n_burst/2)+1) - delay_range - (LLA_focous(floor(n_sl_output/2)+1,3)-LLA_focous(j,3)) - diff_eta(j)*radial_velocity;

    H=exp(-2*1j*pi*cst.alpha*t*2/c.*corr_range + 1j*2*pi*2/cst.c*cst.fc*r_p'); 
    tmp=(data.*H);
    % slc = fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);
    % % tmp2=fftshift(fft((tmp),cst.n_sample_range*range_ext_factor*zp),2);
    % % % slc(j,:) = fftshift(fft(nansum(data.*H),cst.n_sample_range*range_ext_factor*zp),2);

    for mg=1:size(tmp,1)
        tmp1=fftshift(fft((tmp(mg,:)),cst.n_sample_range*range_ext_factor*zp),2);
        power_matrix(mg,:)=tmp1.*conj(tmp1);
    end
power_matrix_total(:,:,j)=power_matrix;
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







n = 11520;
m = segmentation;
P = zeros(m, size(tmp, 2)); % Initialize P with the appropriate size

for i = 1:m
    start_idx = (i-1)*n/m + 1;
    if i < m
        end_idx = i*n/m;
    else
        end_idx = n;
    end
    P(i,:) = fftshift(fft(nansum(tmp(start_idx:end_idx, :)), cst.n_sample_range*range_ext_factor*zp), 2);
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
    w1=(R1(ix)/100)^1.4;  
% w1=1.01^R1(ix);    w1=(w1-1)/((1.01^100)-1); 
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
    slc(j,:)=fftshift(fft(nansum(Total_power),cst.n_sample_range*range_ext_factor*zp),2);
else
    % slc(j,:)=fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);
    slc(j,:)=0;
end


    slc1(j,:)=fftshift(fft(nansum(tmp),cst.n_sample_range*range_ext_factor*zp),2);


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




      


for ix=1:180


    [correlation, lags] = xcorr(aa(ix,:) / max(aa(ix,:)), Ref_Waveform / max(Ref_Waveform), 'coeff');
    [max_corr, max_idx] = max(correlation); % Maximum correlation value
    lag(ix) = abs(lags(max_idx));
    RR(ix)=max_corr*100;
    misfit(ix)=sqrt(sum((aa(ix,:) / max(aa(ix,:)))-Ref_Waveform).^2);
    % tmp1=corrcoef((Ref_Waveform / max(Ref_Waveform)), aa(ix,:) / max(aa(ix,:)));
    R1(ix)=correlation(lags==0)*100;



    if(lag(ix)<3)
        R1(ix)=RR(ix);
        lag(ix)=0;
    end
    w1=(R1(ix)/100)^1.4;  
% w1=1.01^R1(ix);    w1=(w1-1)/((1.01^100)-1); 
w2=(1/(min(misfit(ix))));   w2= (w2-(1/128))/(10^6-(1/128));
w3=(1/(1+(1.4^abs(lag(ix)))));  w3= (w3-(1/(1+(1.4^64))))/(0.5-(1/(1+(1.4^64))));
wgt(ix)= w1  * 1 * w3;





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







for ix=1:180


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
        scalefactor=agc(floor(n_burst/2)+1) + sig0_cal(floor(n_burst/2)+1) + 30*log10(alt(floor(n_burst/2)+1)) + 10*log10(1+alt(floor(n_burst/2)+1)/6378137.0) + 10*log10(4*16*pi^2*cst.fc^2/cst.c^3/256/cst.Fs) -2*cst.GdB;
        
        % Pulse peakiness
        pulse_peakiness=max(PT,[],2)./mean(PT,2) * 85/128;


original_wf=PT1';
sat_alt=alt(CoB);
tracker_range=tracker(CoB);
regenerated_wf=PT2;



end


