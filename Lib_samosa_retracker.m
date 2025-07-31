       function [WL] = Lib_samosa_retracker(wf,tracker_range, sat_alt, correction,GEO)

% %         """
% %          Public Method -> Retrack_Samosa -> method to retrack a SAR (Unfocused) Altimetry waveform by SAMOSA or SAMOSA+ retracker
% % 
% %              Input :
% %                     - self : class self
% %                     - tau : time array  (giving the time of each range gate of the waveform, tau=0 is given at the reference gate)
% %                     - wf : input waveform
% %                     - LookAngles : input Look Angle Array of each Doppler Beams in degrees, set it to None if you dont have this input ( generated looks will be
% %                                    counted between N_Look_min and N_Look_max, both given in CONF)
% %                     - MaskRanges : input Mask Range Array in meter, set it to None if you dont have this input (in this case they will be
% %                                     autonomously computed by the library)
% %                     - GEO  : structure with fields: LAT (latitude in deg), LON (longitude in deg), Height (Orbit Height in m),
% %                              Vs (Satellite Velocity in m/sec), Hrate (Orbit Height Rate in m/sec), Pitch (Altimeter Pitch in radiant),
% %                              Roll (Altimeter Roll in radiant), nu (inverse of mean square slope), ThN (Thermal Noise)
% %                              track_sign (if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=0 in CONF )
% %                     - CONF : structure with fields: flag_slope (flag to include in the model the slope of orbit and surface), wf_weighted (set it to True if waveform is weighted)
% %                              beamsamp_factor (1 means only one beam per resolution cell is generated in the DDM), N_Look_min (number of the first Look to generate in the DDM),
% %                              N_Look_max (number of the last Look to generate in the DDM), guess_epoch (first-guess epoch in second), guess_swh (first-guess swh in m),
% %                              guess_pu (first-guess Pu), guess_nu (first-guess nu), lb_epoch (lower bound on epoch in sec), lb_swh (lower bound in swh in m),
% %                              lb_pu (lower bound on Pu), lb_nu (lower bound on nu), ub_epoch (upper bound on epoch in sec), ub_swh (upper bound in swh in m),
% %                              ub_pu (upper bound on Pu), ub_nu (upper bound on nu), rtk_type (it can be 'samosa' to retrack the waveform with SAMOSA retracker
% %                              or it can be 'samosa+' to retrack the waveform with SAMOSA+ retracker)
% % 
% %              Output :
% % 
% %                     - epoch in seconds
% %                     - SWH in meter
% %                     - Amplitude Pu
% %                     - misfit
% %                     - ocean-like flag (1 means openocean, 0 means non-openocean)
% % 
% %         """
    ref_bin=44;
    FS=320e6;
    SOL=299792458;
    tau_44=-6.5625e-08;



    CST.c0=299792458.0;                   %%% speed of light in m/sec
    CST.R_e=6378137.0;                    %%% Reference Ellipsoid Earh Radius in m
    CST.f_e=1/298.257223563;              %%% Reference Ellipsoid Earth Flatness
    CST.gamma_3_4=1.2254167024651779;     %%% Gamma Function Value at 3/4

        MaskRanges=0;

        CONF.flag_slope = false;                    %%% flag True commands to include in the model the slope of orbit and surface (this effect usually is included in LookAngles Array)
        CONF.beamsamp_factor = 1;                   %%% 1 means only one beam per resolution cell is generated in the DDM, the other ones are decimated
        CONF.wf_weighted = false;                   %%% flag True if the waveform under iteration is weighted
        CONF.N_Look_min = -90;                      %%% number of the first Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.N_Look_max = 90;                       %%% number of the last Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.guess_swh = 2;                         %%% first-guess SWH in meter
        CONF.guess_pu = 1;                          %%% first-guess Pu
        CONF.guess_nu = 2;                          %%% first-guess nu (only used in second step of SAMOSA+)
        CONF.lb_epoch = 0;                          %%% lower bound on epoch in sec. If set to None, lower bound will be set to the first time in input array tau
        CONF.lb_swh = -0.5;                         %%% lower bound on SWH in m
        CONF.lb_pu = 0.2;                           %%% lower bound on Pu
        CONF.lb_nu = 0;                             %%% lower bound on nu (only used in second step of SAMOSA+)
        CONF.ub_epoch = 0;                          %%% upper bound on epoch in sec. If set to None, upper bound will be set to the last time in input array tau
        CONF.ub_swh = 30;                           %%% upper bound on SWH in m
        CONF.ub_pu = 1.5;                           %%% upper bound on Pu
        CONF.ub_nu = 1e9;                           %%% upper bound on nu (only used in second step of SAMOSA+)
        CONF.rtk_type = 'samosa+';                  %%% choose between 'samosa' or 'samosa+'
        CONF.wght_factor= 1.4705;                   %%% widening factor of PTR main lobe after Weighting Window Application


    wf_zp=1;
    RDB.Np_burst=64;                        % number of pulses per burst
    RDB.Npulse=128;                         % number of the range gates per pulse (without zero-padding)
    RDB.PRF_SAR=17825;                      % Pulse Repetition Frequency in SAR mode , given in Hz
    RDB.BRI=0.0127;                         % Burst Repetition Interval, given in sec
    RDB.f_0=13.575e9;                       % Carrier Frequency in Hz
    RDB.Bs=320e6;                           % Sampled Bandwidth in Hz
    RDB.theta_3x=(pi/180)*(1.0766);         % (rad) Antenna 3 dB beamwidth (along-track)
    RDB.theta_3y=(pi/180)*(1.2016);         % (rad) Antenna 3 dB beamwidth (cross-track)





        RDB.PRI_SAR = 1./RDB.PRF_SAR;
        RDB.lambda_0 = CST.c0 / RDB.f_0;
        RDB.dfa = RDB.PRF_SAR / RDB.Np_burst;
        CST.ecc_e = sqrt((2. - CST.f_e)* CST.f_e); % Earth Eccentricty
        CST.b_e = CST.R_e* sqrt(1. - CST.ecc_e^2);

    LUT.F0='LUT_F0.txt';                                                                           %%% filename of the F0 LUT
    LUT.F1='LUT_F1.txt';                                                                           %%% filename of the F1 LUT
    LUT.alphap_noweight='alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_NOHAMMING).txt';          %%% filename of the alphap LUT ( case no weighting)
    LUT.alphap_weight='alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_HAMMING).txt';              %%% filename of the alphap LUT ( case weighting)
    LUT.alphapower_noweight='alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt';      %%% filename of the alpha power LUT ( case no weighting)
    LUT.alphapower_weight='alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt';        %%% filename of the alpha power LUT ( case weighting)

global F0_LUT F1_LUT alphap_LUT_NoWght alphap_LUT_Wght alphapower_LUT_NoWght alphap_LUT_NoWght
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

    Nstart = RDB.Npulse * wf_zp;
    Nend  = RDB.Npulse * wf_zp;
    dt = 1. / (RDB.Bs * wf_zp);                                    %%%%% time sampling step for the array tau, it includes the zero-padding factor
    tau=-(Nstart/2)*dt:dt:((Nend-1)/2)*dt;

look_angle_start_20_hr_ku=-0.0094;
look_angle_stop_20_hr_ku=0.0094;
stack_number_after_weighting_20_hr_ku=180;
LookAngles=90-linspace( (look_angle_start_20_hr_ku)*(180/pi),(look_angle_stop_20_hr_ku)*(180/pi),...
                     floor(stack_number_after_weighting_20_hr_ku) );

    [pks,locs,widths,proms]=findpeaks(wf);
    threshold=0.2*max(pks);
    pks1=pks(pks>threshold);
    loc=locs(ismember(pks,pks1));

    % [~,tt]=min(abs(loc-PLoc));

    CONF.guess_epoch=tau(loc(1));


% % % % 
        if (CONF.lb_epoch == 0)
            CONF.lb_epoch=tau(1)*1e9;
        end
        if (CONF.ub_epoch == 0)
            CONF.ub_epoch=tau(end)*1e9;
        end
% % % % 
% % % %         if CONF.flag_slope:
% % % %             CONF.flag_slope=1
% % % %         else:
% % % %             CONF.flag_slope=0

    NstartNoise = 2;    %%% noise range gate counting from 1, no oversampling
    NendNoise   = 40;    %%% noise range gate counting from 1, no oversampling
    GEO.ThN=compute_ThNEcho(wf',NstartNoise*wf_zp,NendNoise*wf_zp);        %%%% computing Thermal Noise from the waveform matric


        CONF.lb = [CONF.lb_epoch, CONF.lb_swh, CONF.lb_pu];
        CONF.ub = [CONF.ub_epoch, CONF.ub_swh, CONF.ub_pu];
        guess_triplet=[CONF.guess_epoch*1e9,CONF.guess_swh,CONF.guess_pu];
        CONF.step = 1;
        wf_norm = wf ./ max(wf);
                E = -1 * sum(wf_norm .^ 2 .* log2(wf_norm .^ 2), 2);
        PP = 1 ./ sum(wf_norm, 2);
        wf_zp = numel(wf_norm) / RDB.Npulse;
        GEO.ThN_norm = GEO.ThN / max(wf);


% %  [out,resnorm,residual,exitflag,output]  = lsqnonlin(@(guess_triplet)Compute_Residuals(guess_triplet, tau, wf_norm, LookAngles, MaskRanges, GEO),guess_triplet,CONF.lb,CONF.ub);
% % g_max_model=max(residual);

% %             out=scipy.optimize.least_squares(__Compute_Residuals,guess_triplet, bounds=(CONF.lb,CONF.ub),loss=OPT.loss,
% %                                              method=OPT.method, ftol=OPT.ftol, xtol=OPT.xtol, gtol=OPT.gtol,
% %                                              max_nfev=OPT.max_nfev,args=(tau,wf_norm, LookAngles, MaskRanges, GEO))

% %             swh = out(2);
% %             misfit = sqrt(1. / (numel(tau)) * sum(residual .^ 2)) * 100;

            if (CONF.rtk_type == "samosa+")

                if (1 || E * PP < 0.68 || E * PP > 0.78 || (100 * PP) * wf_zp > 8 || (E / misfit) / wf_zp < 4)

                    CONF.lb = [CONF.lb_epoch, CONF.lb_nu, CONF.lb_pu];
                    CONF.ub = [CONF.ub_epoch, CONF.ub_nu, CONF.ub_pu];
                    CONF.step = 2;

                    guess_triplet = [CONF.guess_epoch * 1e9, CONF.guess_nu, CONF.guess_pu];

 [out,resnorm,residual,exitflag,output]  = lsqnonlin(@(guess_triplet)Compute_Residuals(guess_triplet, tau, wf_norm, LookAngles, MaskRanges, GEO),guess_triplet,CONF.lb,CONF.ub);
g_max_model=max(residual);

            swh = out(2);
            misfit = sqrt(1. / (numel(tau)) * sum(residual .^ 2)) * 100;

                end
            end
                                                       

        Pu=(out(3)*max(wf)/g_max_model);
%         out.model=residual+wf_norm;
        oceanlike_flag=~(E * PP < 0.68 || E * PP > 0.78 || (100 * PP) * wf_zp > 8 || (E / misfit) / wf_zp < 4);


        epoch_sec=out(1)*1e-9;
        SWH=swh;
        Pu=Pu;
        misfit=misfit;
        oceanlike_flag=oceanlike_flag;



WL=sat_alt-tracker_range - (epoch_sec-tau_44)*SOL/2 -correction;



        function residuals=Compute_Residuals(guess_triplet, tau, wf_norm, LookAngles, MaskRanges, GEO)

% % %         """
% % %          Private Method -> __Compute_Residuals
% % % 
% % %              Input :
% % % 
% % %                     - self : class self
% % %                     - guess_triplet : triplet of guess epoch (in ns), SWH (or nu for second step of SAMOSA+), and Pu
% % %                     - tau : time array  (giving the time of each range gate of the waveform, tau=0 is given at the reference gate)
% % %                     - wf_norm : input normalized waveform
% % %                     - LookAngles : input Look Angle Array of each Doppler Beam
% % %                     - MaskRanges : input Mask Range Array
% % %                     - GEO : structure with fields ... (see method Retrack_Samosa)
% % % 
% % %              Output :
% % % 
% % %                     - residuals : residuals between model waveform and data waveform
% % %         """

        wf_zp=numel(tau)/RDB.Npulse;
        dr=CST.c0/(2*RDB.Bs*wf_zp);

        earth_radius = sqrt(CST.R_e ^ 2.0 * (cosd((GEO.LAT))) ^ 2 + CST.b_e ^ 2.0 * (sind((GEO.LAT))) ^ 2);
        kappa = (1. + GEO.Height / earth_radius);

        if LookAngles==0

            dtheta = GEO.Vs* RDB.BRI/ ( GEO.Height * kappa );
            Theta1 = pi / 2 + dtheta * CONF.N_Look_min;
            Theta2 = pi / 2 + dtheta * CONF.N_Look_max;
            LookAngles = (180/pi)*(Theta1:dtheta:Theta2);
        end

        DopFreqs = (2*GEO.Vs / RDB.lambda_0) * cosd( (LookAngles) );
        BeamIndex=round(CONF.beamsamp_factor*DopFreqs / RDB.dfa)/CONF.beamsamp_factor;
        span=find(diff(BeamIndex,1)==0);
        BeamIndex(span)=[];
        if (CONF.rtk_type=="samosa" || CONF.step == 1 )

            epoch_ns = guess_triplet(1);
            SWH = guess_triplet(2);
            Pu = guess_triplet(3);

% %             DDM=__Generate_SamosaDDM(epoch_ns,SWH,tau, BeamIndex, GEO)
        DDM=Generate_SamosaDDM(epoch_ns,SWH,tau, BeamIndex, GEO)';
%         DDM=Generate_SamosaDDM(self,tau(60),SWH,tau, BeamIndex, GEO);


        elseif (CONF.rtk_type=="samosa+" && CONF.step == 2)

            epoch_ns = guess_triplet(1);
            SWH = guess_triplet(2);
            Pu = guess_triplet(3);

%             DDM = __Generate_SamosaDDM(epoch_ns, 0, tau, BeamIndex, GEO)
            DDM=Generate_SamosaDDM(epoch_ns, 0, tau, BeamIndex, GEO)';
%         DDM=Generate_SamosaDDM(self,tau(60),SWH,tau, BeamIndex, GEO);

        else

%             print('  SAMOSA Retracker Generation given in input ' + CONF.rtk_type + ' not recognized')
%             return np.nan * np.ones(np.shape(wf_norm))
        end
        if  (MaskRanges == 0) 

            Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR);
            MaskRanges_demin = GEO.Height * ( sqrt(1 + (kappa * ( (Lx * BeamIndex)  / GEO.Height).^2)) - 1 );

        else

% %             MaskRanges = MaskRanges(span)=[];
% %             MaskRanges_demin = MaskRanges - min(MaskRanges)
        end
        R  = repmat( MaskRanges_demin, [numel(wf_norm), 1] )';
        Dr = repmat( dr * (numel(wf_norm):-1:1),[numel(BeamIndex),1] );

        global DDM_Data;

        DDM(find(R >= Dr)) = 0;
        DDM_Data=DDM;

        Pr = sum(DDM', 2) ./ numel(BeamIndex);

        max_model = max(Pr);
        g_max_model=max_model;

        Pr = Pu * (Pr/ max(Pr)) + GEO.ThN_norm;
            genWf=Pr;
        residuals=Pr-squeeze(wf_norm)';


        end



        function DDM=Generate_SamosaDDM(epoch_ns,SWH,tau,l,GEO)

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

        epoch_sec = epoch_ns * 1e-9 ;   %%% epoch (convert back in second)

        earth_radius = sqrt(CST.R_e^2.0 * (cos((GEO.LAT)))^2 + CST.b_e^ 2.0 *(sin((GEO.LAT)))^2);

        tau=tau-epoch_sec;             %%% tau and epoch are both given in seconds
        Dk = (tau*RDB.Bs);
        kappa = (1. + GEO.Height/earth_radius);

        alpha_x = 8. * log(2.0)  / (GEO.Height^2 * RDB.theta_3x^2);
        alpha_y = 8. * log(2.0)  / (GEO.Height^2 * RDB.theta_3y^2);

        if (CONF.wf_weighted==1) 

            if (CONF.step == 1) 

%              s  ind = np.argmin(abs(alphap_LUT_Wght[:, 0] - SWH))
%              s   alpha_p = alphap_LUT_Wght[:, 1][ind]
%              s   Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR)
% 
%              d   ind=np.argmin(abs(alphapower_LUT_Wght[:, 0] - SWH))
%              d   alpha_power=alphapower_LUT_Wght[:, 1][ind]
            elseif (CONF.step == 2) 

                %ind = np.argmin(abs(alphap_LUT_Wght[:, 0] - SWH))
                %alpha_p = alphap_LUT_Wght[:, 1][ind]
                %Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR)

%                 ind = np.argmin(abs(alphap_LUT_NoWght[:, 0] - SWH))
%                 alpha_p = alphap_LUT_NoWght[:, 1][ind]
%                 Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR) * CONF.wght_factor
% 
%                 ind=np.argmin(abs(alphapower_LUT_Wght[:, 0] - SWH))
%                 alpha_power=alphapower_LUT_Wght[:, 1][ind]
            end
        elseif (CONF.wf_weighted==0)

             [~ ,ind] = (min(abs(alphap_LUT_NoWght(:, 1) - SWH)));
             alpha_p = alphap_LUT_NoWght(ind, 2);
             Lx = CST.c0 * GEO.Height / (2. * GEO.Vs * RDB.f_0 * RDB.Np_burst * RDB.PRI_SAR);

            [~,ind]=min(abs(alphapower_LUT_NoWght(:, 1) - SWH));
            alpha_power=alphapower_LUT_NoWght(ind, 2);
        else

%             print('  Waveform Weighting Flag given in input ' + CONF.wf_weighted + ' not recognized')
%             return np.nan * np.ones( (len(Dk),len(l)) )
        end
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

        end










        end

