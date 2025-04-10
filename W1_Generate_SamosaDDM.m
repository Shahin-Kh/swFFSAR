function genWf=W1_Generate_SamosaDDM(epoch_ns,SWH,GEO)


% %         """
% %          Private Method -> __Generate_SamosaDDM
% % 
% %              Input :
% %                     - self : class self
% %                     - epoch_ns : input epoch given in nanoseconds
% %                     - SWH : input SWH in meter
% %                     - tau : time array  (giving the time of each range gate of the waveform, tau=0 is given at the reference gate)
% %                     - l : input Doppler Beam index
% %                     - GEO : structure with fields ... (see method Retrack_Samosa)
% % 
% %              Output:
% % 
% %                     - DDM (Delay-Doppler Map)
% % 
% %         """








    CST.c0=299792458.0;                   %%% speed of light in m/sec
    CST.R_e=6378137.0;                    %%% Reference Ellipsoid Earh Radius in m
    CST.f_e=1/298.257223563;              %%% Reference Ellipsoid Earth Flatness
    CST.gamma_3_4=1.2254167024651779;     %%% Gamma Function Value at 3/4

    RDB.Np_burst=64;                        % number of pulses per burst
    RDB.Npulse=128;                         % number of the range gates per pulse (without zero-padding)
    RDB.PRF_SAR=17825;                      % Pulse Repetition Frequency in SAR mode , given in Hz
    RDB.BRI=0.0127;                         % Burst Repetition Interval, given in sec
    RDB.f_0=13.575e9;                       % Carrier Frequency in Hz
    RDB.Bs=320e6;                           % Sampled Bandwidth in Hz
    RDB.theta_3x=(pi/180)*(1.0766);         % (rad) Antenna 3 dB beamwidth (along-track)
    RDB.theta_3y=(pi/180)*(1.2016);         % (rad) Antenna 3 dB beamwidth (cross-track)

    LUT.F0='LUT_F0.txt';                                                                           %%% filename of the F0 LUT
    LUT.F1='LUT_F1.txt';                                                                           %%% filename of the F1 LUT
    LUT.alphap_noweight='alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_NOHAMMING).txt';          %%% filename of the alphap LUT ( case no weighting)
    LUT.alphap_weight='alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_HAMMING).txt';              %%% filename of the alphap LUT ( case weighting)
    LUT.alphapower_noweight='alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt';      %%% filename of the alpha power LUT ( case no weighting)
    LUT.alphapower_weight='alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt';        %%% filename of the alpha power LUT ( case weighting)

    
    wf_zp=1;
    Nstart = RDB.Npulse * wf_zp;
    Nend  = RDB.Npulse * wf_zp;
    
    dt = 1. / (RDB.Bs * wf_zp);                                    %%%%% time sampling step for the array tau, it includes the zero-padding factor
    tau=-(Nstart/2)*dt:dt:((Nend-1)/2)*dt;



            fileID = fopen('LUT_F0.txt','r');
            F0_LUT= fscanf(fileID,'%f %f',[2 Inf])';

            fileID = fopen('LUT_F1.txt','r');
            F1_LUT = fscanf(fileID,'%f %f',[2 Inf])';

            fileID = fopen('alphap_table_SEN3_09_Nov_2017.txt','r');
            alphap_LUT_NoWght = fscanf(fileID,'%f %f',[2 Inf])';

            fileID = fopen('alphap_table_SEN3_09_Nov_2017.txt','r');
            alphap_LUT_Wght = fscanf(fileID,'%f %f',[2 Inf])';

            fileID = fopen('alphap_table_SEN3_09_Nov_2017.txt','r');
            alphapower_LUT_NoWght = fscanf(fileID,'%f %f',[2 Inf])';
            fileID = fopen('alphap_table_SEN3_09_Nov_2017.txt','r');
            alphap_LUT_NoWght = fscanf(fileID,'%f %f',[2 Inf])';



        CST=CST;
        RDB=RDB;
%         OPT=OPT;
%         CONF = None;

        RDB.PRI_SAR = 1./RDB.PRF_SAR;
        RDB.lambda_0 = CST.c0 / RDB.f_0;
        RDB.dfa = RDB.PRF_SAR / RDB.Np_burst;
        CST.ecc_e = sqrt((2. - CST.f_e)* CST.f_e); % Earth Eccentricty
        CST.b_e = CST.R_e* sqrt(1. - CST.ecc_e^2);

        max_model=1;
        sucess=true;



        % % GEO.LAT=42.842;%lat_20_hr_ku(k);                              %%% latitude in degree for the waveform under iteration
        % % GEO.LON=270.99;%lon_20_hr_ku(k);                              %%% longitude in degree between -180, 180 for the waveform under iteration
        % % GEO.Height=8.0872e+05;%alt_20_hr_ku(k);                           %%% Orbit Height in meter for the waveform under iteration
        % % GEO.Vs=7.5411e+03;%sat_vel_vec_20_hr_ku(k);                       %%% Satellite Velocity in m/s for the waveform under iteration
        % % GEO.Hrate=12.26;%orb_alt_rate_20_hr_ku(k);                   %%% Orbit Height rate in m/s for the waveform under iteration
        % % GEO.Pitch=-8.7266e-06%;off_nadir_pitch_angle_str_20_hr_ku(k);      %%% Altimeter Reference Frame Pitch in radiant
        % % GEO.Roll=1.7453e-06;%off_nadir_roll_angle_str_20_hr_ku(k);        %%% Altimeter Reference Frame Roll in radiant
        GEO.nu=0;                                                         %%% Inverse of the mean square slope
        GEO.track_sign=0;                                                 %%% if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=False in CONF


look_angle_start_20_hr_ku=-0.0094;
look_angle_stop_20_hr_ku=0.0094;
stack_number_after_weighting_20_hr_ku=180;


        LookAngles=90-linspace( (look_angle_start_20_hr_ku)*(180/pi),(look_angle_stop_20_hr_ku)*(180/pi),...
                             floor(stack_number_after_weighting_20_hr_ku) );

        DopFreqs = (2*GEO.Vs / RDB.lambda_0) * cosd( (LookAngles) );
CONF.beamsamp_factor=1;
        BeamIndex=round(CONF.beamsamp_factor*DopFreqs / RDB.dfa)/CONF.beamsamp_factor;  
        span=find(diff(BeamIndex,1)==0);
        BeamIndex(span)=[];
l=BeamIndex;
CONF.flag_slope=0;










        epoch_sec = epoch_ns * 1e-9 ;   %%% epoch (convert back in second)

        earth_radius = sqrt(CST.R_e^2.0 * (cos((GEO.LAT)))^2 + CST.b_e^ 2.0 *(sin((GEO.LAT)))^2);

        tau=tau-epoch_sec;             %%% tau and epoch are both given in seconds
        Dk = (tau*RDB.Bs);
        kappa = (1. + GEO.Height/earth_radius);

        alpha_x = 8. * log(2.0)  / (GEO.Height^2 * RDB.theta_3x^2);
        alpha_y = 8. * log(2.0)  / (GEO.Height^2 * RDB.theta_3y^2);



             [~ ,ind] = (min(abs(alphap_LUT_NoWght(:, 1) - SWH)));
             alpha_p = alphap_LUT_NoWght(ind, 2);
             Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR);

            [~,ind]=min(abs(alphapower_LUT_NoWght(:, 1) - SWH));
            alpha_power=alphapower_LUT_NoWght(ind, 2);

        Ly = sqrt(CST.c0 * GEO.Height / (kappa * RDB.Bs));
        Lz = CST.c0 / (2. * RDB.Bs);
        Lg = kappa/ (2.*GEO.Height*alpha_y);

        sigma_s = (SWH/ (4. * Lz));
        sigma_z = (SWH/ 4.);

        yk = 0 * Dk;
        yk(find(Dk > 0)) = Ly*sqrt(Dk(find(Dk > 0)));

        xl = Lx * l;

        orbit_slope = GEO.track_sign*((CST.R_e^2 - CST.b_e^2)/ (2. * earth_radius^2))* sin((2. * GEO.LAT)) - ...
                      (-GEO.Hrate/GEO.Vs);

        ls = CONF.flag_slope*orbit_slope* GEO.Height/ (kappa* Lx);

        gl = 1./ sqrt(alpha_p^2 + 4. * (alpha_p^2) * (Lx/ Ly)^4 * (l - ls).^2 + sign(SWH)* (SWH/ (4. * Lz))^2);

        csi = gl .* Dk';

        z = 1. / 4. * csi .^ 2;

        xp =  + GEO.Height * GEO.Pitch;
        yp =  - GEO.Height * GEO.Roll;

        Gamma_0 = exp(-alpha_y*yp^2 - alpha_x* (xl - xp).^2 - xl.^2*GEO.nu/ GEO.Height^2 -...
                      (alpha_y + GEO.nu/GEO.Height^2)*yk'.^2).* cosh(2.*alpha_y*yp*yk');
[~,tt]=size(z);
        T_kappa = zeros(size(z));
        T_kappa(find(Dk > 0), :)  = repmat(( (1 + GEO.nu/ ((GEO.Height^2)*alpha_y)) - yp./(Ly* sqrt(Dk(find(Dk > 0)))).* tanh(2*alpha_y*yp* Ly*sqrt(Dk(find(Dk > 0)))))',1,tt);
        T_kappa(find(Dk <= 0),:) = (1 + GEO.nu/ ((GEO.Height^2)*alpha_y)) - 2*alpha_y * yp^2;

        csi_max_F0 = max(F0_LUT(:,1));
        csi_min_F0 = min(F0_LUT(:, 1));
        clip_F0=(csi>=csi_min_F0 & csi<=csi_max_F0);
tmp=csi';
ttmp=size(z);
        f0 = zeros(size(z));
        Index = floor((numel(F0_LUT(:,1)) - 1) * ((tmp(clip_F0') - csi_min_F0) / (csi_max_F0 - csi_min_F0)))+1;
       
        try
        
        tmp2=((tmp(clip_F0') - F0_LUT(Index,1)).*((F0_LUT(Index + 1,2) - F0_LUT(Index,2)) ./(F0_LUT(Index + 1,1) - F0_LUT(Index,1))) + F0_LUT(Index,2));
        
        catch

            2+2

        end


% 
% f0(clip_F0)=tmp2;
k=1;
        for(i1=1:ttmp(1))
            for(i2=1:ttmp(2))
                if(clip_F0(i1,i2)==1)
                    f0(i1,i2)=tmp2(k);
                    k=k+1;
                end
            end
        end
        
        
       
        f0(find(csi>csi_max_F0))=1/2*sqrt(pi)/(z(find(csi>csi_max_F0))).^(1/4)*(1+3/(32*z(find(csi>csi_max_F0)))+105./(2048*(z((csi>csi_max_F0))).^2) + 10395/(196608*(z(find(csi>csi_max_F0))).^3));
        f0(find(csi == 0)) = (1./2.)*(pi*2^(3./4.))/(2.*CST.gamma_3_4);
        f0(find(csi < csi_min_F0))=0;

        csi_max_F1 = max(F1_LUT(:,1));
        csi_min_F1 = min(F1_LUT(:, 1));
        clip_F1=(csi>=csi_min_F1 & csi<=csi_max_F1);



        f1 = zeros(size(z));
        Index = floor((numel(F1_LUT(:,1)) - 1) * ((tmp(clip_F1') - csi_min_F1) / (csi_max_F1 - csi_min_F1)))+1;
        tmp2=(tmp(clip_F1') - F1_LUT(Index,1)).*((F1_LUT(Index + 1,2) - F1_LUT(Index,2)) ./(F1_LUT(Index + 1,1) - F1_LUT(Index,1))) + F1_LUT(Index,2);
        
   k=1;
        for(i1=1:ttmp(1))
            for(i2=1:ttmp(2))
                if(clip_F1(i1,i2)==1)
                    f1(i1,i2)=tmp2(k);
                    k=k+1;
                end
            end
        end     

        f1(find(csi>csi_max_F1))=1/2 * 1/4 *sqrt(pi)/(z(find(csi>csi_max_F1))).^(3/4);
        f1(find(csi == 0)) = -(1./2.)*(2. ^ (3. / 4. ))*(CST.gamma_3_4/2);
        f1(find(csi < csi_min_F1))=0;



        f = (f0 + sigma_z ./ Lg .* T_kappa .* gl .* sigma_s .* f1);

        const=sqrt(2.*pi*alpha_power^4);

        DDM = const*sqrt(gl).*Gamma_0.*f;
        DDM=DDM;







      earth_radius = sqrt(CST.R_e ^ 2.0 * (cosd((GEO.LAT))) ^ 2 + CST.b_e ^ 2.0 * (sind((GEO.LAT))) ^ 2);
        kappa = (1. + GEO.Height / earth_radius);

            Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR);
            MaskRanges_demin = GEO.Height * ( sqrt(1 + (kappa * ( (Lx * BeamIndex)  / GEO.Height).^2)) - 1 );

        dr=CST.c0/(2*RDB.Bs*wf_zp);
        R  = repmat( MaskRanges_demin, [128, 1] )';
        Dr = repmat( dr * (128:-1:1),[numel(BeamIndex),1] );


    DDM(find(R >= Dr)) = 0;
        DDM_Data=DDM;

        Pr = sum(DDM, 2) ./ numel(BeamIndex);

        max_model = max(Pr);
        g_max_model=max_model;
Pu=1;
GEO.ThN_norm=0.0012;
        Pr = Pu * (Pr/ max(Pr)) + GEO.ThN_norm;
            genWf=Pr;








        end