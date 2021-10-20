function f_neonopto_plot_pref_raster(myKsDir,cond_preferred_angle,psth_locs,spike_data,neuron,...
    neuron_latency5p,neuron_latency3,bin_size,area_name,session_name,f1f0,cell_chan,p_values)

neuron_latency5p=neuron_latency5p./1000;
neuron_latency3=neuron_latency3./1000;


% session_name='neon';

outputDir           = fullfile(myKsDir,session_name);
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

% load([myKsDir,'\',session_name,'_spike_data.mat'])
smoothed_pdf=spike_data.smoothed_pdf{neuron};

xp=spike_data.xp;
spike_in_neuron=spike_data.spikes{neuron};
stim_len=spike_data.stim_len;
Ntrial_per_cond=spike_data.trial_per_cond{neuron}(psth_locs);
id=spike_data.evoked_cids{neuron};
% Ntrial_per_cond=trial_per_cond{neuron};
% spike_in_neuron= spikes{neuron};
% smoothed_pdf=smoothed_pdfs{neuron};
%% plot options

Nhp=5; Nvp=2;
raster_loc=[1 2 8 7 5 6 3 4];




if ~isempty(spike_in_neuron)
    figure('units','normalized','outerposition',[0 0 0.8 1],'Visible','off');
    histo=smoothed_pdf(psth_locs,:);
    
    for cond=1:length(psth_locs)%%% loop over all condition for selected cluster and preferred angle
        spike_in_cond=spike_in_neuron(spike_in_neuron(:,3)==psth_locs(cond),1:2);
        subplot(Nhp,Nvp,raster_loc(cond))
        scatter(spike_in_cond(:,1),spike_in_cond(:,2),2,'black','fill') %try this instead of loop
        xlim([xp(1) xp(end)]);
        ylim([0 Ntrial_per_cond(cond)])
        ylabel('trial')
        ax = gca;
        ax.FontSize = 12;
        rectangle('Position',[0 0 stim_len Ntrial_per_cond(cond)],'EdgeColor','g','LineWidth',2)
        if ~isnan(neuron_latency5p(cond))
            xline(neuron_latency5p(cond),'-.r',{'L5%',num2str(neuron_latency5p(cond))},'LineWidth',1,'LabelVerticalAlignment','top','LabelOrientation','horizontal');
        end
        if p_values(cond)<.01 && neuron_latency3(cond)<.3
            xline(neuron_latency3(cond),'-.b',{'L3',num2str(neuron_latency3(cond))},'LineWidth',1,'LabelVerticalAlignment','bottom','LabelOrientation','horizontal');
        end
        switch cond
            case 1
                title(['NCS ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 2
                title(['LDG ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 3
                title(['NCS+ring ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 4
                title(['NCS+square ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 5
                title(['NCS+circle ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 6
                title(['LDG+circle ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 7
                title(['NCS+light ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])
            case 8
                title(['LDG+light ','f1f0=',num2str(f1f0(cond),'%4.2f'),' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')])                
        end
        
        %%%psth plot
        sp=histo(cond,:);
        s1=length(sp);
        m  = s1 - mod(s1, bin_size);
        sp=mean(reshape(sp(1:m),bin_size,[]));
        xpp=mean(reshape(xp(1:m),bin_size,[]));
        switch cond
            case 1
                subplot(Nhp,Nvp,9)
                hold on
                plot(xpp,sp,'-b','LineWidth',2)
            case 2
                subplot(Nhp,Nvp,9)
                plot(xpp,sp,'-r','LineWidth',2)
            case 7
                subplot(Nhp,Nvp,9)
                plot(xpp,sp,'-.c','LineWidth',2);
            case 8
                subplot(Nhp,Nvp,9)
                plot(xpp,sp,'-.m','LineWidth',2);
                ylabel('Firing rate')
                legend('NCS','LDG','NCS+light','LDG+light')
                rectangle('Position',[0 0 stim_len max(histo(:))],'EdgeColor','g','LineWidth',2)
                xlim([xp(1) (xp(end)+xp(end)/4)]);
                ylim([min(histo(:)) max(histo(:))])
                xlabel('time (second)')
                ax = gca;
                ax.FontSize = 12;
            case 3
                subplot(Nhp,Nvp,10)
                hold on
                plot(xpp,sp,'-m','LineWidth',2)
            case 4
                subplot(Nhp,Nvp,10)
                plot(xpp,sp,'-k','LineWidth',2);
            case 5
                subplot(Nhp,Nvp,10)
                hold on
                plot(xpp,sp,'-c','LineWidth',2);
            case 6
                subplot(Nhp,Nvp,10)
                plot(xpp,sp,'-y','LineWidth',2);
                ylabel('Firing rate')
                legend('ring block','square block','circle','LDG+ circle')
                rectangle('Position',[0 0 stim_len max(histo(:))],'EdgeColor','g','LineWidth',2)
                xlim([xp(1) (xp(end)+xp(end)/4)]);
                ylim([min(histo(:)) max(histo(:))])
                xlabel('time (second)')
                ax = gca;
                ax.FontSize = 12;
        end
    end
    sgtitle(['ch:',num2str(cell_chan),' ',area_name,' neuron:',num2str(id)],'FontSize',28)
    saveas(gcf,[outputDir,'\ch',num2str(cell_chan),'_preferred raster_cluster-',num2str(id),'.png'])
        %% save in different folders
    if p_values(1)<0.001 && p_values(2)<0.001
    saveas(gcf,[myKsDir(1:end-15),'pv0.001\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005
        saveas(gcf,[myKsDir(1:end-15),'pv0.005\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    if p_values(1)<0.01 && p_values(2)<0.01
        saveas(gcf,[myKsDir(1:end-15),'pv0.01\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    
    if p_values(1)<0.001 && p_values(2)<0.001 && f1f0(2)<0
    saveas(gcf,[myKsDir(1:end-15),'pv0.001_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005 && f1f0(2)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.005_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    if p_values(1)<0.01 && p_values(2)<0.01 && f1f0(2)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.01_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    
    if p_values(1)<0.001 && p_values(2)<0.001 && f1f0(2)<0 && f1f0(1)<0
    saveas(gcf,[myKsDir(1:end-15),'pv0.001_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005 && f1f0(2)<0 && f1f0(1)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.005_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end 
    if p_values(1)<0.01 && p_values(2)<0.01 && f1f0(2)<0 && f1f0(1)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.01_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.png'])
    end   
    close gcf
end


