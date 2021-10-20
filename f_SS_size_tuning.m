function f_SS_size_tuning(myKsDir)
% Surronud suppression and specific size tunning



% myKsDir                 = 'D:\matwork\Data\M145_2019_07_12';

session_name='size';

outputDir     = fullfile(myKsDir,session_name);
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

load([myKsDir,'\',session_name,'_spike_data.mat'])
inputfile_stm = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)

pdfs=spike_data.pdfs;
xp=spike_data.xp;
% cid=spike_data.evoked_cids;
evoked_cids = cell2mat(spike_data.evoked_cids);
% stim_size =[5:5:30 35:10:55];
% stim_size =[4 6 8 10:5:30 35:10:55 60];
% stim_sizes= [2.5 5:5:45];  %[2.5 5:5:45 55 65];
stim_sizes = stimulus_sizes{2};
interpolated_size = min(stim_sizes):2:max(stim_sizes);
stim_angle=stimulus_angles{2};
Ncond=length(stim_angle)*length(stim_sizes);         %number of condition
[~,ind] = min(abs(interpolated_size-30));


time_window=[0 0.4];                        % averaging time  window
heatmapDim1=length(stim_sizes);
heatmapDim2=length(stim_angle);






%%%
win_ind = find(xp>=time_window(1) & xp<time_window(2));

Nc=length(pdfs);       %number of clusters
sr           = zeros(1,Nc);
sr30         = zeros(1,Nc);
pref_size    = zeros(1,Nc);
size_tuning  = [];
category      = ones(1,Nc);

for neuron=1:Nc
    heatmap=zeros(heatmapDim1,heatmapDim2);
    histo=pdfs{neuron};
    for cond=1:Ncond%%% loop over all condition for selected cluster
        heatmap(cond)=mean(histo(cond,win_ind(1):win_ind(end)));
    end
    heatmap(:,1)=smoothdata(heatmap(:,1),'gaussian',5);%%% In 'smoothdata' the SD is fixed to 1/5th of window width
    heatmap(:,2)=smoothdata(heatmap(:,2),'gaussian',5);%%% In 'smoothdata' the SD is fixed to 1/5th of window width
    heatmap1 = interp1(stim_sizes,heatmap(:,1),interpolated_size,'spline');
    heatmap2 = interp1(stim_sizes,heatmap(:,2),interpolated_size,'spline');
    [m1,I1] = max(heatmap1);
    [m2,I2] = max(heatmap2);
    sr1=(m1-heatmap1(end))/m1;
    sr2=(m2-heatmap2(end))/m2;
    if sr1>=sr2
        sr(neuron)=sr1;
        pref_size(neuron) = interpolated_size(I1);
    else
        sr(neuron)=sr2;
        pref_size(neuron) = interpolated_size(I2);
    end
    m1=max(heatmap1(1:ind)); 
    m2=max(heatmap2(1:ind));
    sr1=(m1-heatmap1(end))/m1;
    sr2=(m2-heatmap2(end))/m2;
    if sr1>=sr2
        sr30(neuron) = sr1;
        tunning=heatmap1;
    else
        sr30(neuron) = sr2;
        tunning=heatmap2;
    end
    
    size_tuning  = [size_tuning;tunning];
    if sr30(neuron)>= -.05 && sr30(neuron)<= 0.05
        category(neuron)=2;
    elseif sr30(neuron)<-0.05
        category(neuron)=3;
    end

end
SR=vertcat(evoked_cids,sr);
SR30=vertcat(evoked_cids,sr30);
if length(heatmap1) == length(interpolated_size)
    stim_sizes = interpolated_size;
end
S_S_data.stim_sizes  = stim_sizes;
S_S_data.SR          = SR;
S_S_data.SR30        = SR30;
S_S_data.pref_size   = pref_size;
S_S_data.size_tuning = size_tuning;
S_S_data.category     = category;
save([myKsDir,'\',session_name,'_S_S_data.mat'],'S_S_data')
