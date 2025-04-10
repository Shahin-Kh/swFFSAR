function DataSet=W1_DataSet_Explore2(dataset)

Sentinel3_Data=dataset;

DataSet.times=dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time;

%% Corrections

Dry_Tropo=dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal;
Wet_Tropo=dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal;
Ionosphere=dataset.ObjVS.Gen.Cor.Atm.Iono.AltPLRM.Hi.Ku.Signal;
Inverse_Barometer=0;%dataset.ObjVS.Gen.Cor.Atm.InvBar.ECMWF.Hi.Signal;

SS_Bias=0;%dataset.ObjVS.Gen.Cor.Target.SSBias.Hi.Ku.Signal;
Ocean_Tide=dataset.ObjVS.Gen.Cor.Target.OceanTide.Sol1.Hi.Signal;
Earth_Tide=dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal;
Polar_Tide=dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal;

Geoid_Correction=Sentinel3_Data.ObjVS.Gen.GeoH.XGM2019;


   DataSet.SumOfCorrections = ...
   Dry_Tropo + Wet_Tropo + Ionosphere + Inverse_Barometer + ...
   SS_Bias + Ocean_Tide + Earth_Tide + Polar_Tide + Geoid_Correction;


DataSet.wf=normalize(Sentinel3_Data.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal,'scale',2);
DataSet.Krts=kurtosis(DataSet.wf,1,2);
DataSet.sknss=skewness(DataSet.wf,1,2);



DataSet.Tracker_Range=Sentinel3_Data.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal;

DataSet.OCOG_ReTracker_Range=Sentinel3_Data.ObjVS.Raw.Mes.Rng.OCOG.SAR.Hi.Ku.Signal;
DataSet.Satellite_Altitude=Sentinel3_Data.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal;
DataSet.SSH_OCOG = DataSet.Satellite_Altitude - DataSet.OCOG_ReTracker_Range - DataSet.SumOfCorrections;

DataSet.Latitude_cycle=Sentinel3_Data.ObjVS.Raw.Sat.Lat.Hi.Ku.Signal;
DataSet.Longitude_cycle=Sentinel3_Data.ObjVS.Raw.Sat.Lon.Hi.Ku.Signal-360;
DataSet.BackScatterCoff=dataset.ObjVS.Raw.Mes.Wvf.Sig0.OCOG.SAR.Hi.Ku.Signal;
% DataSet.Peakiness=dataset.ObjVS.Raw.Mes.Wvf.Pknss.Hi.Ku.Signal;


for cnt=1:numel(DataSet.wf(:,1))
    [pks,locs,widths,proms]=findpeaks((DataSet.wf(cnt,:)));
    threshold=0.6*max(pks);
    pks1=pks(pks>threshold);
    loc=locs(ismember(pks,pks1));
    width=widths(ismember(pks,pks1));
    DataSet.firstPeak(cnt,:)=pks1(1);
    DataSet.firstPeakLoc(cnt,:)=loc(1);

    [MaxPeak locMax]=max(pks1);
    DataSet.MaxPeak(cnt,:)=MaxPeak;
    DataSet.MaxPeakLoc(cnt,:)=loc(locMax);
    DataSet.MaxPeakWidth(cnt,:)=width(locMax);
    
    
    DataSet.NoPeaks(cnt,:)=numel(pks1);
    DataSet.peakplacement(cnt,:)=DataSet.MaxPeakLoc(cnt,:)-DataSet.firstPeakLoc(cnt,:);
    
    [NoisLev NoisLevNormalized]=NoiseLevel(DataSet.wf(cnt,:));
    DataSet.NoisLev(cnt,:)=NoisLev;
    DataSet.NoisLevNormalized(cnt,:)=NoisLevNormalized;

    DataSet.COG(cnt,:)=sum((1:128).*(DataSet.wf(cnt,:).^2))/(sum(DataSet.wf(cnt,:).^2));
    DataSet.Peakiness(cnt,:)=max(DataSet.wf(cnt,:))/sum(DataSet.wf(cnt,:));

    
end





%% AGC
% Signal_AGC=Sentinel3_Data.ObjVS.Raw.Mes.Wvf.AGC.SAR.Lo.Ku.Signal;
% DataSet.Signal_AGC=Signal_AGC;
DataSet.cycls=dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Cycle;
cycls=unique(dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Cycle);
AGC_mean=nanmean(dataset.ObjVS.Raw.Mes.Wvf.AGC.SAR.Lo.Ku.Signal  );
for i=1:numel(cycls)
    ptrs=find(cycls(i)==dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Cycle);
        if(ismember(cycls(i),dataset.ObjVS.Raw.Mes.Wvf.AGC.SAR.Lo.Ku.Cycle))
            k=find(cycls(i)==dataset.ObjVS.Raw.Mes.Wvf.AGC.SAR.Lo.Ku.Cycle);
            tmp=dataset.ObjVS.Raw.Mes.Wvf.AGC.SAR.Lo.Ku.Signal(k);
            DataSet.AGC(ptrs)=tmp(1);
            if(isnan(tmp(1)))
                DataSet.AGC(ptrs)=AGC_mean;%mean(DataSet.AGC);
            end

        else
            DataSet.AGC(ptrs)=AGC_mean;%mean(DataSet.AGC);
        end
end


DataSet.AGC(ptrs)=DataSet.AGC(ptrs)';

%% Itaparica Reservoir
% % % % % % load('BelemDoSaoFrancisco.mat')
% % % % % % ref=BelemDoSaoFrancisco;
% % % % % % 
% % % % % % for cnt=1:numel(DataSet.wf(:,1))
% % % % % % 
% % % % % %         time=floor(datenum(Tag_Time_AltBundle(dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time(cnt))));
% % % % % %         ptr=find(time==ref(:,4));
% % % % % %     try
% % % % % %         DataSet.ref(cnt,1)=ref(ptr,5);
% % % % % %     catch
% % % % % %         DataSet.ref(cnt,1)=0;
% % % % % %     end
% % % % % % 
% % % % % % 
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % DataSet.Y=zeros(numel(DataSet.wf(:,1)),1);
% % % % % % DataSet.Error=abs(DataSet.ref-DataSet.SSH_OCOG);
% % % % % % 
% % % % % % for cnt=1:numel(DataSet.wf(:,1))
% % % % % %     if(DataSet.Error(cnt)>2)% || (DataSet.MaxPeakLoc(cnt,:)>10)%  || DataSet.peakplacement(cnt,:)>5) 
% % % % % %         DataSet.Y(cnt,1)=1;
% % % % % %         if(DataSet.ref(cnt)==0)
% % % % % %             if(DataSet.SSH_OCOG(cnt)<300)
% % % % % %                 DataSet.Y(cnt,1)=1;
% % % % % %             else
% % % % % %                  DataSet.Y(cnt,1)=0;
% % % % % %             end
% % % % % %         end
% % % % % %     end
% % % % % % end


for cnt=1:numel(cycls)
ptrs=find(cycls(cnt)==dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Cycle);

            % RetrackError=(DataSet.Tracker_Range(ptrs)-((DataSet.MaxPeakLoc(ptrs)-44)/320e6)*299792458/2);
RetrackError=DataSet.Satellite_Altitude(ptrs)-(DataSet.Tracker_Range(ptrs)-((DataSet.firstPeakLoc(ptrs)-44)/320e6)*299792458/2)-DataSet.SumOfCorrections(ptrs);
DataSet.DeviationFromMean(ptrs)=abs(RetrackError-mean(RetrackError));

% DataSet.DeviationFromMean(ptrs)=abs(DataSet.MaxPeakLoc(ptrs)-mean(DataSet.MaxPeakLoc(ptrs)));
end
DataSet.DeviationFromMean=DataSet.DeviationFromMean';


end

