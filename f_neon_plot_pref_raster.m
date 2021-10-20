function f_neon_plot_pref_raster(myKsDir,preferred_angle,psth_locs,spike_data,neuron, neuron_latency4,neuron_latency2,f1f0,bin_size)

neuron_latency4=neuron_latency4./1000;
neuron_latency2=neuron_latency2./1000;


session_name='neon';

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

        Nhp=3; Nvp=2;
        raster_loc=[ 1 2 4 3];
        



if ~isempty(spike_in_neuron)
    figure('units','normalized','outerposition',[0 0 0.8 1],'Visible','off');
    histo=smoothed_pdf(psth_locs,:);

    for cond=1:4%%% loop over all condition for selected cluster and preferred angle
        spike_in_cond=spike_in_neuron(spike_in_neuron(:,3)==psth_locs(cond),1:2);
        subplot(Nhp,Nvp,raster_loc(cond))
        for k=1:Ntrial_per_cond(cond) %% loop over all trails in a condition for raster plot
            Tri=spike_in_cond(spike_in_cond(:,2)==k,1);
            if ~isempty(Tri)
                y=ones(1,length(Tri)).*(k);
                scatter(Tri,y,4,'black','fill')
                hold on
            end
        end
%         set(gca,'xticklabel',[]);
%         set(gca,'YTick',[]);
        xlim([xp(1) xp(end)]);
        ylim([0 Ntrial_per_cond(cond)])
        ylabel('trial')
        ax = gca;
        ax.FontSize = 22;
        rectangle('Position',[0 0 stim_len k],'EdgeColor','g','LineWidth',2)
        if ~isnan(neuron_latency4(cond))
            xline(neuron_latency4(cond),'-.r',{'L4',num2str(neuron_latency4(cond))},'LineWidth',1,'LabelVerticalAlignment','top','LabelOrientation','horizontal');
        end
        if ~isnan(neuron_latency2(cond))
            xline(neuron_latency2(cond),'-.b',{'L2',num2str(neuron_latency2(cond))},'LineWidth',1,'LabelVerticalAlignment','bottom','LabelOrientation','horizontal');
        end
        switch cond
            case 1
                title(['NIG, , f1f0=',num2str(f1f0(cond))])
%                 title(['neon, f1f0=',num2str(f1f0(cond))])
            case 2
                title(['LDG-1, f1f0=',num2str(f1f0(cond))])
%                 title('Luminance defined grating-1')
            case 3
%                 title('small circle')
                title(['LDG-2, f1f0=',num2str(f1f0(cond))])
            case 4
                title('control')
        end

        %%%psth plot
        subplot(Nhp,Nvp,5:6)
        hold on
        sp=histo(cond,:);
        s1=length(sp);
        m  = s1 - mod(s1, bin_size);
        sp=mean(reshape(sp(1:m),bin_size,[]));
        xpp=mean(reshape(xp(1:m),bin_size,[]));
        switch cond
            case 1
                plot(xpp,sp,'-.b','LineWidth',2)
            case 2
                plot(xpp,sp,'-r','LineWidth',2)
            case 3
%                 plot(xpp,sp,'-m','LineWidth',2)
            case 4
                plot(xpp,sp,'-k','LineWidth',2);
        end
        legend('NIG','LDG','ctrl')
        legend('boxoff')
        rectangle('Position',[0 0 stim_len max(histo(:))],'EdgeColor','g','LineWidth',2)
%         xticks(xp(1):.05:xp(end))
%         set(gca,'xticklabel',[]);
%         set(gca,'YTick',[]);
        xlim([xp(1) xp(end)]);
        ylim([min(histo(:)) max(histo(:))])
        xlabel('time (second)')
        ylabel('Spike/sec')
        
    end
    sgtitle(['cluster:',num2str(id),'   preferred angle:',num2str(preferred_angle)],'FontSize',28)
    set(gca,'linewidth',1.5)
    ax = gca;
    ax.FontSize = 22;
    ax.FontName = "Arial";
    saveas(gcf,[outputDir,'\',num2str(bin_size),session_name,'_preferred raster_cluster-',num2str(id),'.png'])
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);
    saveas(gcf,[outputDir,'\',num2str(bin_size),session_name,'_preferred raster_cluster-',num2str(id),'.pdf'])
    close gcf
end

    
