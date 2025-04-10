clear all
clc

OCOG = 1;
SAMOSA = 0;


folder = "Lake Rathbun";

% Get a list of all .mat files in the folder
filePattern = fullfile(folder, '*.mat');
matFiles = dir(filePattern);

% Initialize temporary variables for the current case study
In_Situ = [];
SSH_OCOG_Wf_ORG = [];
SSH_OCOG_Wf_ENH = [];
SSH_INSITU=[];

% Loop through each .mat file
for k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    load(fullFileName);

        SSH_OCOG_Wf_ORG=[SSH_OCOG_Wf_ORG Result.ORG_OCOG_WL];
        SSH_OCOG_Wf_ENH=[SSH_OCOG_Wf_ENH Result.REG_OCOG_WL];
        SSH_INSITU=[SSH_INSITU Result.ref];
end

           
bias=nanmedian(SSH_OCOG_Wf_ENH-SSH_INSITU);
Metric_ORG = Metrics(SSH_OCOG_Wf_ORG',SSH_INSITU'+bias);
Metric_ENH = Metrics(SSH_OCOG_Wf_ENH',SSH_INSITU'+bias);



figure;
hold on;
plot(SSH_OCOG_Wf_ORG)
plot(SSH_OCOG_Wf_ENH)
plot(SSH_INSITU+bias)

