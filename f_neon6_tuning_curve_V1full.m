function f_neon6_tuning_curve_V1full(myKsDir)

session_name='Fneon';

outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end

load([myKsDir,'\',session_name,'_spike_data']) %#ok<LOAD>
pdfs=spike_data.pdfs;
xp=spike_data.xp;
% cid=spike_data.cluster_id;
cid=spike_data.evoked_cids;

stim_len=spike_data.stim_len;




tuningYlabel='Average firing rate';
tuningXlabel='Grating angles';

stim_angle=[0 45 90 135 180 225 270 315];
stim_type=[1 2];
neon_type=[1 2 3];
Ncond=length(stim_type)*length(stim_angle)*length(neon_type);     %number of condition

heatmapDim1=length(stim_angle);
heatmapDim2=length(stim_type)*length(neon_type);
       

time_window=[0 stim_len/2];                 % averaging time  window 




win_ind = find(xp>=time_window(1) & xp<time_window(2));

Nc=length(pdfs);       %number of clusters

for neuron=1:Nc %%%loop over clusters
    heatmap=zeros(heatmapDim1,heatmapDim2);
    histo=pdfs{neuron};
    for cond=1:Ncond%%% loop over all condition for selected cluster
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% average firing rate
        heatmap(cond)=mean(histo(cond,win_ind(1):win_ind(end)));
    end
    figure('units','normalized','outerposition',[0 0 0.5 0.6]);set(gcf, 'Visible', 'off');
    plot(heatmap(:,1),'-.bd','LineWidth',2);
    hold on
    plot(heatmap(:,2),'-r*','LineWidth',2);
    plot(heatmap(:,3),'-m*','LineWidth',2);
    plot(heatmap(:,4),'-kd','LineWidth',2);
    plot(heatmap(:,5),'-c*','LineWidth',2);
    plot(heatmap(:,6),'-yd','LineWidth',2);
%     legend('neon-color','physical-1','physical-2','ctrl','Location','northwest')
    legend('neon-color','luminance','Ring block','square block','dynamic circle','Lum+circle','Location','northwest')

    title(['direction tuning curve cluster:',num2str(cid{neuron})]);
    xlabel(tuningXlabel);ylabel(tuningYlabel);
    
%     set(gca,'XTick',stim_angle)
        set(gca,'xticklabel',stim_angle)
        ax = gca;
        ax.FontSize = 18;
    saveas(gcf,[outputDir,'\',session_name,'_tuning curve_win-',num2str((time_window(2)-time_window(1))*1000),'ms_cluster-',num2str(cid{neuron}),'.png'])
    close gcf

end
    
