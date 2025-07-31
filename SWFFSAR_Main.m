clear;
clc;
% close all;
warning off;

ref_bin=87;
Fs=320e6;
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
CaseStudies="Susquehanna River at Harrisburg, PA - 01570500";   %Good
% CaseStudies="Cheyenne River near Eagle Butte SD - 06439500";
% CaseStudies="COOSA River at Leesburg - 02399500";
% CaseStudies="Mohawk River at Vischer Ferry Dam NY 01356000";
% CaseStudies="Chattahoochee River - 0234296910";
% CaseStudies="Lake Abiquiu";
% 
if(CaseStudies=="COOSA River at Leesburg - 02399500")
    COOSA='F:\Sentinel-3-Level1\Altbundle data\S3B COOSA River\x2685S3BCOOSARIVERATLEESBURGAL_p3201t549o110sr3000wo40or99.mat';
    load(COOSA);
    l2_dataset=x2685S3BCOOSARIVERATLEESBURGAL_p3201t549o110sr3000wo40or99;
    load('02399500.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];
elseif(CaseStudies=="Rathbun Lake near Rathbun, IA - 06903880")
    load('S3ARathbun_p3101t635o110sr2000wo70or99.mat');
    l2_dataset = S3ARathbun_p3101t635o110sr2000wo70or99;
    load('06903880.mat')
    insitu_ref=USGS_Data{1,2};
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];
elseif(CaseStudies=="Lake Tahoe near Tahoe City CA - 10337000")
    Tahoe='F:\Sentinel-3-Level1\Altbundle data\S3BTahoe\S3BTahoe_p3201t409o110sr20000wo70or99.mat';
    load(Tahoe);
    l2_dataset=S3BTahoe_p3201t409o110sr20000wo70or99;
    outlier_cycles=[];

    load('10337000.mat')
    insitu_ref=USGS_Data{1,2};
    USGS_POS=[39.18071, -120.1193368];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);

elseif(CaseStudies=="Susquehanna River at Harrisburg, PA - 01570500")
    SusquehannaRiver='F:\Sentinel-3-Level1\Altbundle data\x1003S3A01570500SusquehannaRiverAt\x1003S3A01570500SusquehannaRiverAt_p3101t263o110sr3000wo70or99.mat';
    load(SusquehannaRiver);
    l2_dataset=x1003S3A01570500SusquehannaRiverAt_p3101t263o110sr3000wo70or99;
    % % outlier_cycles=1:40;
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];

    load('01570500.mat')
    insitu_ref=USGS_Data{1,2};
    USGS_POS=[40.25481164, -76.8860846];
elseif(CaseStudies=="Columbia River at Stevenson, WA - 14128600")
    ColombiaRiver='F:\Sentinel-3-Level1\Altbundle data\x9828S3A14128600COLUMBIARIVERATST\x9828S3A14128600COLUMBIARIVERATST_p3101t140o110sr3000wo30or0.mat';
    load(ColombiaRiver);
    l2_dataset=x9828S3A14128600COLUMBIARIVERATST_p3101t140o110sr3000wo30or0;

    load('14128600.mat')
    insitu_ref=USGS_Data{1,2};
    USGS_POS=[45.6992835, -121.8684125];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];

    % % outlier_cycles=1:40;
elseif(CaseStudies=="Boca Res NR Truckee CA - 10344490")
    Boca='F:\Sentinel-3-Level1\Altbundle data\S3BBoca\S3BBoca_p3201t409o110sr3000wo70or99.mat';
    load(Boca);
    l2_dataset=S3BBoca_p3201t409o110sr3000wo70or99;

    load('10344490.mat')
    insitu_ref=USGS_Data{1,2};
    USGS_POS=[39.3887, -120.0962];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];

    % % outlier_cycles=1:40;
elseif(CaseStudies=="LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500")
    MOHAVE='F:\Sentinel-3-Level1\Altbundle data\S3AMohave\S3AMohave_p3101t95o110sr3000wo80or99.mat';
    load(MOHAVE);
    l2_dataset=S3AMohave_p3101t95o110sr3000wo80or99;
    load('09422500.mat')
    insitu_ref=USGS_Data{2,2};
    USGS_POS=[ 	35.1963472,	-114.5702028];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];

elseif(CaseStudies=="Lake Koshkonong Near Newville, WI - 05427235")
    Koshkonong='F:\Sentinel-3-Level1\Altbundle data\S3AKoshkonong\S3AKoshkonong_p3101t435o110sr3000wo70or99.mat';
    load(Koshkonong);
    l2_dataset=S3AKoshkonong_p3101t435o110sr3000wo70or99;

    load('05427235.mat')
    insitu_ref=USGS_Data{2,2};
    USGS_POS=[	42.85750645, 	-88.9409429];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    outlier_cycles=[];

elseif(CaseStudies=="Itaparica reservoir")
    load('S3AItaparica_p3101t375o110sr2500wo70or99.mat');
    load('BelemDoSaoFrancisco.mat')
    insitu_ref=BelemDoSaoFrancisco;
    l2_dataset=S3AItaparica_p3101t375o110sr2500wo70or99;
    USGS_POS=[-9.139, -38.312];
    outlier_cycles=[4 9 16 18 28 39 41];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
elseif(CaseStudies=="Lake Abiquiu")
    Abiquiu='F:\Sentinel-3-Level1\Altbundle data\S3A Abiquiu\S3AAbiquiu_p3101t351o110sr3000wo80or99.mat';
    load(Abiquiu);
    l2_dataset=S3AAbiquiu_p3101t351o110sr3000wo80or99;
    USGS_POS=[36.244 , -106.436];
    outlier_cycles=[];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    load('Data_USGSst_08286900_2015_2024.mat')
    insitu_ref=Data_USGSst_08286900_2015_2024;
    insitu_ref(:,4)=insitu_ref(:,1);
    insitu_ref(:,5)=insitu_ref(:,2);
elseif(CaseStudies=="Chattahoochee River - 0234296910")
    Chattahoochee='F:\Sentinel-3-Level1\Altbundle data\Chattahoochee_p3201t549\Chattahoochee_p3201t549o110sr3000wo90or99.mat';
    load(Chattahoochee);
    l2_dataset=Chattahoochee_p3201t549o110sr3000wo90or99;
    USGS_POS=[31.908 , -85.141];
    outlier_cycles=[];
    [l1a_files, l1a_cycles] = Lib_Load_S3A_L1A_Data_products(CaseStudies);
    load('Data_USGSst_08286900_2015_2024.mat')
    insitu_ref=Data_USGSst_08286900_2015_2024;
    insitu_ref(:,4)=insitu_ref(:,1);
    insitu_ref(:,5)=insitu_ref(:,2);
end





% Load L1A Data Products
meas_number=1;
%% =============== Main Loop: Process Each Available Cycle =============
for cycle_idx = 1:numel(l1a_cycles)
    current_cycle = l1a_cycles(cycle_idx);

    if ~ismember(current_cycle, l2_dataset.ObjVS.Gen.Sat.Lat.Hi.C.Cycle) || ismember(current_cycle, outlier_cycles)
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
    inv_bar       = zeros(numel(ptrs),1);%l2_dataset.ObjVS.Gen.Cor.Atm.InvBar.ECMWF.Hi.Signal(ptrs);

    ocean_tide    = zeros(numel(ptrs),1);%l2_dataset.ObjVS.Gen.Cor.Target.OceanTide.Sol1.Hi.Signal(ptrs);
    earth_tide    = l2_dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal(ptrs);
    pole_tide     = l2_dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal(ptrs);
    
    waveforms     = l2_dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal(ptrs, :);
    peakiness     = l2_dataset.ObjVS.Raw.Mes.Wvf.Pknss.SAR.Hi.Ku.Signal(ptrs);
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


   
        peak_locs=[];
        for pt_idx = 1:numel(ptrs)
                [pks, locs, widths, proms] = findpeaks(waveforms(pt_idx,:), 'SortStr', 'descend');
                peak_locs(pt_idx)=locs(1);
        end



        % % % estimated_WL_elevation = raw_elevation - correction - ...
        % % %                       ((peak_locs' - ref_bin)) * (1/Fs) * (SOL / 2);

        estimated_WL_elevation = raw_elevation_OCOG  ;
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



            peak_loc =(raw_elevation(pt_idx)-correction-nanmedian(estimated_WL_elevation))/((1/Fs) * (SOL / 4))+ref_bin;
              [~, closest_idx] = min(abs(peak_loc - locs));
              peak_loc = locs(closest_idx);
            
              
              
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


        % Generate reference waveform
        if peakiness(pt_idx) > 20
            RefWave = W1_sinc2_ref(waveforms(pt_idx,:),peak_loc,256);
            RefWave = downsample(RefWave,2);
            % RefWave = Lib_sinc2_ref(peak_loc, 128);
        else

            


            RefWave = W1_Generate_SamosaDDM((-66 + (peak_loc/2)) * 3.125, SWH, GEO)';

                if(isnan(RefWave(1,1)))
                    RefWave = W1_sinc2_ref(waveforms(pt_idx,:),peak_loc,256);
                    RefWave = downsample(RefWave,2);
                    % RefWave = Lib_sinc2_ref(peak_loc, 128);
                end

        end

        % Run waveform processing
        [original_wf, regenerated_wf, tracker_range, scalefactor, sat_alt] = Lib_sw_FFSAR_S3(current_l1a_file, latitudes(pt_idx), longitudes(pt_idx), RefWave, 5);
        % figure;hold on;plot(original_wf/max(original_wf));plot(regenerated_wf/max(regenerated_wf));plot(RefWave/max(RefWave))
        % Plot results
        Result.WF_FFSAR(pt_idx,:)=original_wf;
        Result.WF_swFFSAR(pt_idx,:)=regenerated_wf;


        % ploting_waveform(original_wf/max(original_wf), regenerated_wf/max(regenerated_wf), RefWave);
        pause(0.01);
try
        Result.ORG_OCOG_WL(pt_idx)=Lib_ocog_retracker(original_wf/max(original_wf), tracker_range, sat_alt, correction);
            Result.ORG_SAMOSA_WL(pt_idx)=Lib_samosa_retracker((original_wf/max(original_wf)),tracker_range, sat_alt, correction,GEO);
        Result.REG_OCOG_WL(pt_idx)=Lib_ocog_retracker(regenerated_wf/max(regenerated_wf), tracker_range, sat_alt, correction);
            Result.REG_SAMOSA_WL(pt_idx)=Lib_samosa_retracker(regenerated_wf/max(regenerated_wf),tracker_range, sat_alt, correction,GEO);
catch
    continue
end

      meas_number=meas_number+1;
      




        fclose('all');
    end

% figure;hold on;plot(Result.ORG_OCOG_WL);plot(Result.REG_OCOG_WL);plot(Result.ref);
!!pause(0.01);
    save(fullfile(CaseStudies,strcat(sprintf('%03d', l1a_cycles(cycle_idx)),'_',string(burst_datetime))),"Result");
    clear Result
    clc
    close all
    fclose('all');

end
% CaseStudies="Columbia River at Stevenson, WA - 14128600"
Lib_plot_results(CaseStudies,"ALL","OCOG") 
Lib_plot_results(CaseStudies,"Cycle","SAMOSA") 


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


Rathbun_Radargrams
Lib_plot_Radargrams("Rathbun_Radargrams") 
