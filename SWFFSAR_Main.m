clear;
clc;
close all;
warning off;

%% ======================= Load Data Products ==========================
CaseStudy = "Lake Rathbun";

% Load L2 Data Products
load('Rathbun_p3101t635o001wo85or99.mat');
l2_dataset = Rathbun_p3101t635o001wo85or99;

% Load In-situ Measurements
load('Rathbun.mat');
insitu_ref = insitu_data;

ref_bin=44;
Fs=320e6;
SOL=299792458;

% Load L1A Data Products
[l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products_Rathbun();

%% =============== Main Loop: Process Each Available Cycle =============
for cycle_idx = 36:numel(l1a_cycles)

    current_cycle = l1a_cycles(cycle_idx);

    if ~ismember(current_cycle, l2_dataset.ObjVS.Gen.Sat.Lat.Hi.C.Cycle)
        continue;
    end

    % Find pointers to data for current cycle
    ptrs = find(current_cycle == l2_dataset.ObjVS.Gen.Sat.Lat.Hi.C.Cycle);

    % Extract L2 Dataset Parameters
    latitudes     = l2_dataset.ObjVS.Gen.Sat.Lat.Hi.C.Signal(ptrs);
    longitudes    = l2_dataset.ObjVS.Gen.Sat.Lon.Hi.C.Signal(ptrs);
    xgm           = l2_dataset.ObjVS.Gen.GeoH.XGM2019(ptrs);
    dry_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal(ptrs);
    wet_tropo     = l2_dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal(ptrs);
    iono_alt      = l2_dataset.ObjVS.Gen.Cor.Atm.Iono.AltPLRM.Hi.Ku.Signal(ptrs);
    inv_bar       = l2_dataset.ObjVS.Gen.Cor.Atm.InvBar.ECMWF.Hi.Signal(ptrs);

    ocean_tide    = l2_dataset.ObjVS.Gen.Cor.Target.OceanTide.Sol1.Hi.Signal(ptrs);
    earth_tide    = l2_dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal(ptrs);
    pole_tide     = l2_dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal(ptrs);
    
    waveforms     = l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal(ptrs, :);
    peakiness     = l2_dataset.ObjVS.Raw.Mes.Wvf.Pknss.SAR.Hi.Ku.Signal(ptrs);
    altitude      = l2_dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal(ptrs);
    range         = l2_dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal(ptrs);
    raw_elevation = altitude - range;
    
    % Time processing
    time_stamps = datenum(Lib_Tag_Time_AltBundle(l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time(ptrs)));
    time_days = floor(time_stamps);

    % Read L1A file for current cycle
    current_l1a_file = string(l1a_files(cycle_idx));
    l1a_lat = ncread(current_l1a_file, 'lat_l1a_echo_sar_ku');
    l1a_lon = ncread(current_l1a_file, 'lon_l1a_echo_sar_ku');

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
    time_burst = ncread(current_l1a_file, 'UTC_day_l1a_echo_sar_ku', burst_indices(1), 1);
    burst_datetime = Tag_Time(time_burst);

    bounds = [lat_min, lat_max, lon_min, lon_max];


        correction = dry_tropo + wet_tropo + ...
                     iono_alt + ...
                     ocean_tide + earth_tide + pole_tide + ...
                     xgm;
        peak_locs=[];
        for pt_idx = 1:numel(ptrs)
                [pks, locs, widths, proms] = findpeaks(waveforms(pt_idx,:), 'SortStr', 'descend');
                peak_locs(pt_idx)=locs(1);
        end



        estimated_WL_elevation = raw_elevation - correction - ...
                              (peak_locs' - ref_bin) * (1/Fs) * (SOL / 2);


    % =================== Process Each Point in Cycle ===================
    for pt_idx = 1:numel(ptrs)
        fprintf("Processing waveform %d/%d\n", pt_idx, numel(ptrs));
     [m k(pt_idx)]=min(CartesianDistance( [l1a_lat l1a_lon] , [latitudes(pt_idx) longitudes(pt_idx)] ));


        % Peak detection
        [pks, locs, widths, proms] = findpeaks(waveforms(pt_idx,:), 'SortStr', 'descend');
        peak_thresh = 0.2 * max(pks);
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
        time_match = find(time_days(pt_idx) == insitu_ref(:,4));
        Result.ref(pt_idx) = insitu_ref(time_match, 5);

        % Choose peak location
        if numel(locs) > 1
            [~, closest_idx] = min(abs(corrected_elevation - median(estimated_WL_elevation)));
            peak_loc = locs(closest_idx);
        else
            peak_loc = locs;
        end

        % Generate reference waveform
        if peakiness(pt_idx) > 20
            RefWave = Lib_sinc2_ref(peak_loc, 128);
        else

            alt(pt_idx)=ncread(current_l1a_file,'alt_l1a_echo_sar_ku',k(pt_idx),1);
            orb_alt_rate(pt_idx)=ncread(current_l1a_file,'orb_alt_rate_l1a_echo_sar_ku',k(pt_idx),1);
            vxs(pt_idx)=ncread(current_l1a_file,'x_vel_l1a_echo_sar_ku',k(pt_idx),1);
            vys(pt_idx)=ncread(current_l1a_file,'y_vel_l1a_echo_sar_ku',k(pt_idx),1);
            vzs(pt_idx)=ncread(current_l1a_file,'z_vel_l1a_echo_sar_ku',k(pt_idx),1);
            V(pt_idx)=norm([vxs(pt_idx),vys(pt_idx),vzs(pt_idx)]);
            PITCH(pt_idx)=ncread(current_l1a_file,'pitch_sral_mispointing_l1a_echo_sar_ku',k(pt_idx),1);
            ROLL(pt_idx)=ncread(current_l1a_file,'roll_sral_mispointing_l1a_echo_sar_ku',k(pt_idx),1);
            epoch_ns = -28.1250;
            SWH = 2;
            GEO.LAT = l1a_lat(pt_idx);
            GEO.LON = l1a_lon(pt_idx);
            GEO.Height = altitude(pt_idx);
            GEO.Vs = V(pt_idx);              % Define or load this if missing
            GEO.Hrate = orb_alt_rate(pt_idx); % Define or load this if missing
            GEO.Pitch = PITCH(pt_idx);        % Define or load this if missing
            GEO.Roll = ROLL(pt_idx);          % Define or load this if missing
            GEO.nu = 0;
            GEO.track_sign = 0;



            RefWave = W1_Generate_SamosaDDM((-66 + peak_loc) * 3.125, SWH, GEO)';
        end

        % Run waveform processing
        [original_wf, regenerated_wf, tracker_range, scalefactor, sat_alt] = ...
            Lib_sw_FFSAR(current_l1a_file, latitudes(pt_idx), longitudes(pt_idx), RefWave, 5);

        % Plot results
        % ploting_waveform(original_wf/max(original_wf), regenerated_wf/max(regenerated_wf), RefWave);

        Result.ORG_OCOG_WL(pt_idx)=Lib_ocog_retracker(original_wf, tracker_range, sat_alt, correction);
        % Result.ORG_SAMOSA_WL(pt_idx)=SAMOSA(original_wf);
        Result.REG_OCOG_WL(pt_idx)=Lib_ocog_retracker(regenerated_wf, tracker_range, sat_alt, correction);
        % Result.REG_SAMOSA_WL(pt_idx)=SAMOSA(regenerated_wf);

        fclose('all');
    end

    save(fullfile(CaseStudy,strcat(sprintf('%03d', l1a_cycles(cycle_idx)),'_',string(burst_datetime))),"Result");
    clear Result
    clc
    close all
    fclose('all');

end

Lib_plot_results

