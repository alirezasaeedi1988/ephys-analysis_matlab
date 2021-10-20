function f_spikeFeature(myKsDir)


% myKsDir       = 'D:\matwork\Data\M145_2019_07_09';  

%% Loading kilosort/phy data
sp = loadKSdir(myKsDir);

% good_cid = sp.cgs==2;
good_cid = sp.cids;%(good_cid); for all cluster

    
    
%% Computing some useful properties of the spikes and templates
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
clusterDepths=zeros(length(good_cid),1);
for i=1:length(good_cid)%%%loop over clusters
    
    cid  = sp.spikeCluster==good_cid(i);
    spikes_in_cluster =  spikeDepths(cid);
    clusterDepths(i)=mean(spikes_in_cluster);
end
tempPerClu = findTempForEachClu( sp.spikeCluster, sp.spikeTemplates);
tempPerClu = rmmissing(tempPerClu);
% ma = max(abs(tempPeakWF),[],2);
% tempPeakWF=tempPeakWF./ma;
[mi,mi_ind]= min(tempPeakWF,[],2);
% ma = zeros(length(mi),1);

tpl = zeros(length(mi),1);
for i=1:length(mi)
    [~,tpl(i)]= max(tempPeakWF(i,mi_ind(i):end));
end

tpl=tpl(tempPerClu+1);

ma = max(tempPeakWF,[],2);
tpr=ma./abs(mi);
tpr=tpr(tempPerClu+1);
Ntpr=(ma-mi)./(ma+mi);
Ntpr=Ntpr(tempPerClu+1);

% [spikeWidths, tempWidths, clusterWidths] = computeSpikeWidths(tempsUnW, sp.spikeTemplates, sp.spikeCluster);

% clusterDepths = clusterAverage(sp.clu, spikeDepths);

%% Loading raw waveforms

% gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
% apD = dir(fullfile(myKsDir, '*.dat')); % AP band file from spikeGLX specifically
% gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
% gwfparams.dataType = 'int16';             % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];               % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                     % Number of waveforms per unit to pull out
% gwfparams.spikeTimes = ceil(sp.spikeTime*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = sp.spikeCluster;
% 
% wf = getWaveForms(gwfparams);
sp.spikeAmps   = spikeAmps;
sp.spikeDepths = spikeDepths;
% sp.spikeWidths = spikeWidths;

sp.templateYpos = templateYpos;
sp.tempAmps     = tempAmps;
sp.tempsUnW     = tempsUnW;
sp.tempDur      = tempDur;
sp.tempPeakWF   = tempPeakWF;
% sp.tempWidths   = tempWidths;

sp.clusDepths    = clusterDepths;
sp.clusTPlatency = tpl;
sp.clusTPratio   = tpr;
sp.NclusTPratio   = Ntpr;

% sp.clusWidths    = clusterWidths;
sp.tempPerClu   = tempPerClu;
save([myKsDir,'\spike_Feature.mat'],'sp')

