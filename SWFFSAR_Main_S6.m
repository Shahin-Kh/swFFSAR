clear;
clc;
close all;
warning off;

ref_bin=50;
Fs=395e6;
SOL=299792458;


%% ======================= Load Data Products ==========================
% CaseStudies="Boca Res NR Truckee CA - 10344490";
% CaseStudies="LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500";
% CaseStudies="Rathbun Lake near Rathbun, IA - 06903880";
% CaseStudies="Itaparica reservoir";
% CaseStudies="Lake Tahoe near Tahoe City CA - 10337000";
% CaseStudies="Lake Koshkonong Near Newville, WI - 05427235";
% CaseStudies="Columbia River at Stevenson, WA - 14128600";         %Good
% CaseStudies="Ohio River at Cannelton Dam at Cannelton, IN - 03303280";
% CaseStudies="Susquehanna River at Harrisburg, PA - 01570500";   %Good
% CaseStudies="Cheyenne River near Eagle Butte SD - 06439500";
% CaseStudies="COOSA River at Leesburg - 02399500";
% CaseStudies="Mohawk River at Vischer Ferry Dam NY 01356000";
% CaseStudies="Mendota";
% CaseStudies="Monona";
% CaseStudies="Waubesa";
CaseStudies="Kegonsa";
% 
% USGS_Downloader("425715089164700",'E:\FFSARm Paper\Git_Code\')

if(CaseStudies=="Monona")
    Monona='F:\Sentinel-6-Level1\L2_AltBundle\Sentinel6Monona\Sentinel6Monona_p2501t178\Sentinel6Monona_p2501t178o110sr2000wo80or99.mat';
    load(Monona);
    l2_dataset=Sentinel6Monona_p2501t178o110sr2000wo80or99;
    % % [L1A, L1B, L2]=Load_S6MF_L1A_Data_products_Mendota();
    load('05429000.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products("Mendota");
    outlier_cycles=[];
elseif(CaseStudies=="Waubesa")
    Waubesa='F:\Sentinel-6-Level1\L2_AltBundle\Sentinel6LakeWaubesa\Sentinel6LakeWaubesa_p2501t178\Sentinel6LakeWaubesa_p2501t178o110sr2000wo80or99.mat';
    load(Waubesa);
    l2_dataset=Sentinel6LakeWaubesa_p2501t178o110sr2000wo80or99;
    % % [L1A, L1B, L2]=Load_S6MF_L1A_Data_products_Mendota();
    load('05429485.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products("Mendota");
    outlier_cycles=[];
elseif(CaseStudies=="Mendota")
    Mendota='F:\Sentinel-6-Level1\L2_AltBundle\Sentinel6Mendota\Sentinel6Mendota_p2501t178\Sentinel6Mendota_p2501t178o110sr3000wo80or99.mat';
    load(Mendota);
    l2_dataset=Sentinel6Mendota_p2501t178o110sr3000wo80or99;
    % % [L1A, L1B, L2]=Load_S6MF_L1A_Data_products_Mendota();
    load('05428000.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products("Mendota");
    outlier_cycles=[];
elseif(CaseStudies=="Kegonsa")
    Kegonsa='F:\Sentinel-3-Level1\Altbundle data\Sentinel6LakeKegonsa\Sentinel6LakeKegonsa_p2501t178\Sentinel6LakeKegonsa_p2501t178o110sr3000wo80or99.mat';
    load(Kegonsa);
    l2_dataset=Sentinel6LakeKegonsa_p2501t178o110sr3000wo80or99;
    % % [L1A, L1B, L2]=Load_S6MF_L1A_Data_products_Mendota();
    load('425715089164700.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products("Mendota");
    outlier_cycles=[];
end





% Load L1A Data Products

%% =============== Main Loop: Process Each Available Cycle =============
for cycle_idx = 24:25
    current_cycle = l1a_cycles(cycle_idx);

    if ~ismember(current_cycle, l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Cycle  ) || ismember(current_cycle, outlier_cycles)
        continue;
    end
current_cycle
    % Find pointers to data for current cycle
    ptrs = find(current_cycle == l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Cycle);

    % Extract L2 Dataset Parameters
    % % latitudes     = l2_dataset.ObjVS.Gen.Sat.Lat.Hi.C.Signal(ptrs);
    % % longitudes    = l2_dataset.ObjVS.Gen.Sat.Lon.Hi.C.Signal(ptrs);
    % % xgm           = l2_dataset.ObjVS.Gen.GeoH.XGM2019(ptrs);
    % % dry_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal(ptrs);
    % % wet_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal(ptrs);
    % % iono_alt      = l2_dataset.ObjVS.Gen.Cor.Atm.Iono.AltPLRM.Hi.Ku.Signal(ptrs);

    latitudes     = l2_dataset.ObjVS.Raw.Sat.Lat.Hi.Ku.Signal(ptrs);
    longitudes    = l2_dataset.ObjVS.Raw.Sat.Lon.Hi.Ku.Signal(ptrs);
    xgm           = l2_dataset.ObjVS.Gen.GeoH.XGM2019(ptrs);
    dry_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal(ptrs);
    wet_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal(ptrs);
    % iono_alt      = l2_dataset.ObjVS.Gen.Cor.Atm.Iono.AltPLRM.Hi.Ku.Signal(ptrs);
    iono_alt      = l2_dataset.ObjVS.Gen.Cor.Atm.Iono.Alt.Hi.Signal(ptrs);

    inv_bar       = zeros(numel(ptrs),1);%l2_dataset.ObjVS.Gen.Cor.Atm.InvBar.ECMWF.Hi.Signal(ptrs);

    ocean_tide    = zeros(numel(ptrs),1);%l2_dataset.ObjVS.Gen.Cor.Target.OceanTide.Sol1.Hi.Signal(ptrs);
    earth_tide    = l2_dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal(ptrs);
    pole_tide     = l2_dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal(ptrs);
    
    waveforms     = l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal(ptrs, :);
    peakiness     = 60;%l2_dataset.ObjVS.Raw.Mes.Wvf.Pknss.SAR.Hi.Ku.Signal(ptrs);
    altitude      = l2_dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal(ptrs);
    range         = l2_dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal(ptrs);
    range_ocog    = l2_dataset.ObjVS.Raw.Mes.Rng.OCOG.SAR.Hi.Ku.Signal(ptrs);

     correction = dry_tropo + wet_tropo + ...
                     iono_alt + ...
                     ocean_tide + earth_tide + pole_tide + ...
                     xgm;


for ii=1:numel(ptrs)
    if(isnan(range_ocog(ii)))
        [~,range_ocog(ii)] =Lib_ocog_retracker_Aux(waveforms(ii,:)/max(waveforms(ii,:)), range(ii), 0, 0);
    end
end

if isnan(nanmean(range_ocog))
    continue;
end

    raw_elevation_OCOG = altitude - range_ocog - correction;
    raw_elevation = altitude - range;

% % % shifts=round((raw_elevation-mean(raw_elevation))/(SOL/(Fs*2)));
% % % for i=1:numel(shifts)
% % %     shifted_wf(i,:)=circshift(waveforms(i,:),-shifts(i));
% % % end

    
    % Time processing
    time_stamps = datenum(Lib_Tag_Time_AltBundle(l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time(ptrs)));
    time_days = floor(time_stamps);

    % Read L1A file for current cycle
    current_l1a_file = string(l1a_files(cycle_idx));
    l1a_lat=ncread(current_l1a_file,'/data_140/ku/latitude');
    l1a_lon=ncread(current_l1a_file,'/data_140/ku/longitude');



    % Define bounding box around the L2 points
    margin = 0.0032;
    lat_min = min(latitudes) - margin;
    lat_max = max(latitudes) + margin;
    lon_min = min(longitudes) - margin;
    lon_max = max(longitudes) + margin;

    % Find matching L1A points in the bounding box
    in_bbox = (l1a_lat >= lat_min & l1a_lat <= lat_max) & ...
              (l1a_lon >= lon_min & l1a_lon <= lon_max);
    burst_indices = find(in_bbox);

    % Time of burst (first index)
    time_burst=ncread(current_l1a_file,'/data_140/ku/time',burst_indices(1), 1);
    burst_datetime = Tag_Time_S6(time_burst);

    bounds = [lat_min, lat_max, lon_min, lon_max];


   
        peak_locs=[];
        for pt_idx = 1:numel(ptrs)
                [pks, locs, widths, proms] = findpeaks(waveforms(pt_idx,:), 'SortStr', 'descend');
                peak_locs(pt_idx)=locs(1);
        end



        % % % estimated_WL_elevation = raw_elevation - correction - ...
        % % %                       ((peak_locs' - ref_bin)) * (1/Fs) * (SOL / 2);

        estimated_WL_elevation = raw_elevation_OCOG  ;


   for kk=1:numel(ptrs)
     [m k(kk)]=min(CartesianDistance( [l1a_lat l1a_lon] , [latitudes(kk) longitudes(kk)] ));
   
        alt(kk)=ncread(current_l1a_file,'data_140/ku/altitude',k(kk),1);
        orb_alt_rate(kk)=ncread(current_l1a_file,'data_140/ku/altitude_rate',k(kk),1);
        V(kk)=norm(ncread(current_l1a_file,'data_140/ku/velocity_vector',[1 k(kk)],[3 1]));

        PITCH(kk)=ncread(current_l1a_file,'data_140/ku/off_nadir_pitch_angle_pf',k(kk),1);
        ROLL(kk)=ncread(current_l1a_file,'data_140/ku/off_nadir_roll_angle_pf',k(kk),1);
   end



    % =================== Process Each Point in Cycle ===================
    for pt_idx = 1:numel(ptrs)
        fprintf("Processing waveform %d/%d\n", pt_idx, numel(ptrs));
     [m k(pt_idx)]=min(CartesianDistance( [l1a_lat l1a_lon] , [latitudes(pt_idx) longitudes(pt_idx)] ));


        % Peak detection
        [pks, locs, widths, proms] = findpeaks(waveforms(pt_idx,:), 'SortStr', 'descend');
        peak_thresh = 0.01 * max(pks);
        valid_peaks = pks(pks > peak_thresh);
        locs = locs(ismember(pks, valid_peaks));

        % Apply corrections to elevation
        correction = dry_tropo(pt_idx) + wet_tropo(pt_idx) + ...
                     iono_alt(pt_idx) + ...
                     ocean_tide(pt_idx) + earth_tide(pt_idx) + pole_tide(pt_idx) + ...
                     xgm(pt_idx);
                 
        corrected_elevation = raw_elevation(pt_idx) - correction - ...
                              (locs - ref_bin) * (1/Fs) * (SOL / 2);
        % Match with in-situ data
        
        [~ ,time_match]=min(abs(time_days(pt_idx) - insitu_ref(:,4)));
        
        
    Result.ref(pt_idx) = insitu_ref(time_match, 5);

        % Choose peak location
        % % % if numel(locs) > 1
        % % %     % % % [~, closest_idx] = min(abs(corrected_elevation - mean(estimated_WL_elevation)));
        % % %     % peak_loc = locs(closest_idx);
        % % %     peak_loc =(raw_elevation(pt_idx)-median(estimated_WL_elevation))/((1/Fs) * (SOL / 4))+ref_bin;
        % % %     [~, closest_idx] = min(abs(peak_loc - locs));
        % % %     peak_loc = locs(closest_idx);
        % % % 
        % % % else
        % % %     % % peak_loc = locs;
        % % %     peak_loc =(raw_elevation(pt_idx)-median(estimated_WL_elevation))/((1/Fs) * (SOL / 4))+ref_bin;
        % % %     [~, closest_idx] = min(abs(peak_loc - locs));
        % % %     peak_loc = locs(closest_idx);
        % % % end


% % % data = estimated_WL_elevation(~isnan(estimated_WL_elevation));
% % % 
% % % % Compute the first and third quartiles
% % % Q1 = prctile(data, 30);
% % % Q3 = prctile(data, 70);
% % % IQR = Q3 - Q1;
% % % 
% % % % Define outlier thresholds
% % % lowerBound = Q1 - 1.5 * IQR;
% % % upperBound = Q3 + 1.5 * IQR;
% % % 
% % % % Filter out outliers
% % % filtered_data = data(data >= lowerBound & data <= upperBound);
% % % 
% % % % Compute the mean without outliers
% % % mean_without_outliers = mean(filtered_data);

peakiness = calculate_peakiness_S6(waveforms(pt_idx,:));
            peak_loc =(raw_elevation(pt_idx)-correction-nanmedian(estimated_WL_elevation))/((1/Fs) * (SOL / 4))+256;

            % peak_loc =-(raw_elevation(pt_idx)-correction-nanmedian(estimated_WL_elevation))/((1/Fs) * (SOL / 4))+15;
              [~, closest_idx] = min(abs(peak_loc - locs));
              peak_loc = locs(closest_idx);

         epoch_ns=-28.1250;
        SWH=2;
        GEO.LAT=latitudes((pt_idx));%lat_20_hr_ku(k);                              %%% latitude in degree for the waveform under iteration
        GEO.LON=longitudes((pt_idx));%lon_20_hr_ku(k);                              %%% longitude in degree between -180, 180 for the waveform under iteration
        GEO.Height=alt(pt_idx);%alt_20_hr_ku(k);                           %%% Orbit Height in meter for the waveform under iteration
        GEO.Vs=V(pt_idx);%sat_vel_vec_20_hr_ku(k);                       %%% Satellite Velocity in m/s for the waveform under iteration
        GEO.Hrate=orb_alt_rate(pt_idx);%orb_alt_rate_20_hr_ku(k);                   %%% Orbit Height rate in m/s for the waveform under iteration
        GEO.Pitch=PITCH(pt_idx);%;off_nadir_pitch_angle_str_20_hr_ku(k);      %%% Altimeter Reference Frame Pitch in radiant
        GEO.Roll=ROLL(pt_idx);%off_nadir_roll_angle_str_20_hr_ku(k);        %%% Altimeter Reference Frame Roll in radiant

      
        GEO.nu=0;                                                         %%% Inverse of the mean square slope
        GEO.track_sign=0;                                                 %%% if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=False in CONF







        % Generate reference waveform
        if peakiness > 20
            RefWave = W1_sinc2_ref(waveforms(pt_idx,:),peak_loc/2,256);
            % % % figure;hold on;plot(RefWave/,max(RefWave));plot(waveforms(pt_idx,:)/max(waveforms(pt_idx,:)))
            % RefWave = downsample(RefWave,2);
            % RefWave = Lib_sinc2_ref(peak_loc, 128);
        else

            % % % % alt(pt_idx)=ncread(current_l1a_file,'alt_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % orb_alt_rate(pt_idx)=ncread(current_l1a_file,'orb_alt_rate_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % vxs(pt_idx)=ncread(current_l1a_file,'x_vel_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % vys(pt_idx)=ncread(current_l1a_file,'y_vel_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % vzs(pt_idx)=ncread(current_l1a_file,'z_vel_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % V(pt_idx)=norm([vxs(pt_idx),vys(pt_idx),vzs(pt_idx)]);
            % % % % PITCH(pt_idx)=ncread(current_l1a_file,'pitch_sral_mispointing_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % ROLL(pt_idx)=ncread(current_l1a_file,'roll_sral_mispointing_l1a_echo_sar_ku',k(pt_idx),1);
            % % % % epoch_ns = -28.1250;
            % % % % SWH = 2;
            % % % % GEO.LAT = l1a_lat(pt_idx);
            % % % % GEO.LON = l1a_lon(pt_idx);
            % % % % GEO.Height = altitude(pt_idx);
            % % % % GEO.Vs = V(pt_idx);              % Define or load this if missing
            % % % % GEO.Hrate = orb_alt_rate(pt_idx); % Define or load this if missing
            % % % % GEO.Pitch = PITCH(pt_idx);        % Define or load this if missing
            % % % % GEO.Roll = ROLL(pt_idx);          % Define or load this if missing
            % % % % GEO.nu = 0;
            % % % % GEO.track_sign = 0;
            % % % % 
            % % % % 
            % % % % 
            % % % % RefWave = W1_Generate_SamosaDDM((-66 + (peak_loc/2)) * 3.125, SWH, GEO)';
    RefWave=W1_Generate_SamosaDDM_S6(((peak_loc/2)-128)*2.5316455,SWH,GEO)';
% % figure;hold on;plot(RefWave/max(RefWave));plot(waveforms(pt_idx,:)/max(waveforms(pt_idx,:)))
                if(isnan(RefWave(1,1)))
                    RefWave = W1_sinc2_ref(waveforms(pt_idx,:),peak_loc/2,256);
                    RefWave = downsample(RefWave,2);
                    % RefWave = Lib_sinc2_ref(peak_loc, 128);
                end

        end

        % Run waveform processing
         % [Wf,Wf_Enh,tracker,scale_factor,SatAlt]=W2_FFSAR_S6(current_l1a_file,latitudes(pt_idx), longitudes(pt_idx),RefWave);

        [original_wf, regenerated_wf, tracker_range, scalefactor, sat_alt] = Lib_sw_FFSAR_S6(current_l1a_file, latitudes(pt_idx), longitudes(pt_idx), RefWave, 7);
        % figure;hold on;plot(original_wf/max(original_wf));plot(regenerated_wf/max(regenerated_wf));plot(RefWave/max(RefWave))
        % Plot results
                % OrgWave = downsample(waveforms(pt_idx,:),2);
   
        Result.WF_FFSAR(pt_idx,:)=original_wf;
        Result.WF_swFFSAR(pt_idx,:)=regenerated_wf;

% % figure;hold on;plot(OrgWave/max(OrgWave));plot(original_wf/max(original_wf));plot(regenerated_wf/max(regenerated_wf));plot(RefWave/max(RefWave))
        % ploting_waveform(original_wf/max(original_wf), regenerated_wf/max(regenerated_wf), RefWave);
        pause(0.01);

        Result.ORG_OCOG_WL(pt_idx)=Lib_ocog_retracker_S6(original_wf/max(original_wf), tracker_range, sat_alt, correction);
        Result.ORG_SAMOSA_WL(pt_idx)=Lib_samosa_retracker_S6((original_wf/max(original_wf)),tracker_range, sat_alt, correction,GEO);
        Result.REG_OCOG_WL(pt_idx)=Lib_ocog_retracker_S6(regenerated_wf/max(regenerated_wf), tracker_range, sat_alt, correction);
        Result.REG_SAMOSA_WL(pt_idx)=Lib_samosa_retracker_S6((regenerated_wf/max(regenerated_wf)),tracker_range, sat_alt, correction,GEO);



        fclose('all');
    end
bias=median(Result.ref-Result.REG_SAMOSA_WL);
% figure;hold on;plot(Result.ORG_SAMOSA_WL);plot(Result.REG_SAMOSA_WL);plot(Result.ref-bias);
pause(0.01);
    save(fullfile(CaseStudies,strcat(sprintf('%03d', l1a_cycles(cycle_idx)),'_',string(burst_datetime))),"Result");
    clear Result
    clc
    close all
    fclose('all');

end
% CaseStudies="Columbia River at Stevenson, WA - 14128600"
Lib_plot_results(CaseStudies,"ALL","OCOG") 
Lib_plot_results(CaseStudies,"Cycle","SAMOSA") 

% CaseStudies="Mendota";
% CaseStudies="Monona";
% CaseStudies="Waubesa";
% CaseStudies="Kegonsa";

CaseStudies="Boca Res NR Truckee CA - 10344490";
CaseStudies="LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500";
CaseStudies="Rathbun Lake near Rathbun, IA - 06903880";
CaseStudies="Itaparica reservoir";
CaseStudies="Lake Tahoe near Tahoe City CA - 10337000";
CaseStudies="Lake Koshkonong Near Newville, WI - 05427235";
CaseStudies="Columbia River at Stevenson, WA - 14128600";         %Good
CaseStudies="Ohio River at Cannelton Dam at Cannelton, IN - 03303280";
CaseStudies="Susquehanna River at Harrisburg, PA - 01570500";   %Good
CaseStudies="Cheyenne River near Eagle Butte SD - 06439500";
CaseStudies="COOSA River at Leesburg - 02399500";
CaseStudies="Mohawk River at Vischer Ferry Dam NY 01356000";