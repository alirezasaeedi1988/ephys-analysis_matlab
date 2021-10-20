%%%% statistics
%% parameters
close all
clear
tic
data_folder    = 'D:\matwork\Data\neon_opto_laser\new';
%mouse_numbers  = [340 401 405 407];
mouse_numbers  = [417 419 420 462 463 471];
session_name   = 'neon';
session_num    = 5; %for neon
plot_raster    = 0;  %%% plot raster for preferred angle?
plot_power     = 0;  %%% plot power spectrum
stim_angle     = [0 45 90 135 180 225 270 315];
stim_type      = [1 2];
neon_type      = [1 2 3 4];
Ncond          = length(stim_type)*length(stim_angle)*length(neon_type); 
angle_loc      = reshape(1:Ncond,length(stim_angle),length(stim_type)*length(neon_type));
area_names     = {'V1';'LM';'LI';'LL'};
ave_act        = [];
ave_act_sweep  = [];
ave_act_feed   = [];
ave_act2       = [];
phase_angle    = [];
first_comp     = [];
zero_comp      = [];
f1f0           = [];
f2f1f0         = [];
f2f1           = [];
latency4       = [];
latency2       = [];
latency5p      = [];
preferred_ang  = [];
exited         = [];
cid            = [];
OSI            = [];
DSI            = [];
ave_latency    = [];
pv_simplicity  = [];
h_simplicity   = [];
S_S            = [];
Depths         = [];
TPlatency      = [];
PTratio        = [];
F1_dominant    = [];
RF_size        = []; 
RF_loc         = [];
RF_label       = [];
RF_size_area   = [];
Intrsct_ratio  = [];
RF_size_area_ellips   = [];
Intrsct_ratio_ellips  = [];
RF_RSS         = [];
Cell_area      = [];
detecting_chan = [];
H_values       = [];
P_values       = [];
tag            = cell(2,1);
dummyAllneuron = 0;
grand_psth_ldg  = [];
grand_psth_ncs  = [];
grand_psth_square= [];
opto_psth_ldg   = [];
opto_psth_ncs   = [];

for mouse_number=1:length(mouse_numbers)
    
    %% getting forders for analysis
    %%% Kilosort foldels
    mouse     = ['M',num2str(mouse_numbers(mouse_number))];
    KS_dirs   = cell(1,1);
    folders   = dir(data_folder);
    j=1;
    for i=1:length(folders)
        if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
            KS_dirs{j,1}=fullfile(data_folder,folders(i).name);
            j=j+1;
        end
    end
    
    
    
    
    %%
    for day=1:length(KS_dirs) %% loop for different days
        
        %% loading
        inputfile_stm           = [KS_dirs{day},KS_dirs{day}(end-15:end),'_stimulus_time.mat'];
        load(inputfile_stm)
        load([KS_dirs{day},'\',session_name,'_spike_data.mat'])
        load([KS_dirs{day},'\spike_Feature.mat'])
        %% asigning  previously calculated variables in each day
        tf             = stimulus_TemporalFreq{session_num};          % temporal frequency (cycle/s)
        cycles         = stimulus_cycles{session_num};
        sourc_folder   = fullfile(KS_dirs{day},session_name);
        spikes         = spike_data.spikes;
        pdfs           = spike_data.pdfs;
        RFs            = spike_data.RFs;
        dx             = spike_data.dx;
        smoothed_pdfs  = spike_data.smoothed_pdf;
        xp             = spike_data.xp;
        stim_len       = spike_data.stim_len;
        %surround_sup  = spike_data.surround_suppression;
        evoked_cids    = cell2mat(spike_data.evoked_cids);
        evoked_ind     = ismember(sp.cids,evoked_cids);
        depth          = sp.clusDepths(evoked_ind);
        RF_cell        = f_RF_select(KS_dirs{day},depth,0,'RF');

        %% concatenating
        Depths             = vertcat(Depths,depth);
        TPlatency          = vertcat(TPlatency,sp.clusTPlatency(evoked_ind)./30);
        PTratio            = vertcat(PTratio,sp.clusTPratio(evoked_ind));
        cid                = vertcat(cid,cell2mat(spike_data.evoked_cids)'); %#ok<*AGROW>
        ave_latency        = vertcat(ave_latency,cell2mat(spike_data.average_latency)');
        num_loc            = regexp(KS_dirs{day},'_+[0-9]');
        tarikh             = ['D',KS_dirs{day}(num_loc(3)+1:num_loc(3)+2),'/',KS_dirs{day}(num_loc(2)+1:num_loc(2)+2)];
        RF_size            = vertcat(RF_size,RF_cell.RF_size);
        RF_loc             = vertcat(RF_loc,RF_cell.RF_loc);
        RF_label           = vertcat(RF_label,RF_cell.RF_label);
        RF_size_area       = vertcat(RF_size_area,RF_cell.RF_area);
        RF_size_area_ellips= vertcat(RF_size_area_ellips,RF_cell.RF_area_ellips);
        Cell_area          = vertcat(Cell_area,RF_cell.area);  %% area: 1 = V1, 2 = LM , 3 = LI , 4 = LL
        Intrsct_ratio      = vertcat(Intrsct_ratio,RF_cell.intrsct_ratio);
        Intrsct_ratio_ellips= vertcat(Intrsct_ratio_ellips,RF_cell.intrsct_ratio_ellips);
        RF_RSS             = vertcat(RF_RSS,RF_cell.RF_RSS);
        [~, onset_ind]     = min(abs(xp));
        [~,offset_ind]     = min(abs(xp-stim_len));
        [~,halfset_ind]    = min(abs(xp-stim_len/2));
        Nc                 = length(pdfs);       %number of clusters
        %% preallocating some intermediate variables
        ave_activity_cond  = zeros(Nc,8);
        ave_sweep_cond     = zeros(Nc,8);
        ave_feed_cond      = zeros(Nc,8);
        ave_activity_cond2 = zeros(Nc,8);
        phase_angle_cond   = zeros(Nc,8);
        first_comp_cond    = nan(Nc,8);
        zero_comp_cond     = nan(Nc,8);
        f1f0_cond          = nan(Nc,8);
        f2f1f0_cond        = nan(Nc,8);
        f2f1_cond          = nan(Nc,8);
        latency_cond4      = nan(Nc,8);
        latency_cond2      = nan(Nc,8);
        latency_cond_5p    = nan(Nc,8);

        response_type      = zeros(Nc,1);
        dominant_f1        = zeros(Nc,1);
        h_values           = nan(Nc,8);
        p_values           = nan(Nc,8);
        %% some params
        transient_time = 150; %in ms to exclude the transient component
        transient_nbins= transient_time/(dx*1000);
        bin_size=10;
        FS=1/(dx*bin_size); % freq. of sampling in PSTH
        
        %%single neuron calculation
        for neuron=1:Nc        %%%loop over clusters in one day
            if plot_power
                neuron_name   = [mouse,'_D',KS_dirs{day}(num_loc(3)+1:num_loc(3)+2),'-',KS_dirs{day}(num_loc(2)+1:num_loc(2)+2),'_ID',num2str(spike_data.evoked_cids{neuron})];
                PS=figure('visible','off');
                cond_names = {'NCS', 'LDG','NCS+ring','NCS+squar','NCS+circle','LDG+circle','NCS+light','LDG+light'};
                cond_loc_plot = [1 2 7 8 6 5 3 4];
            end
            dummyAllneuron      = dummyAllneuron+1;
            tag{dummyAllneuron} = [mouse,tarikh,'id',num2str(spike_data.evoked_cids{neuron})];
            spike        = spikes{neuron};
            Pdf          = pdfs{neuron};
            Pdf2         = pdfs{neuron};
            RF           = RFs{neuron};
            smoothed_pdf = smoothed_pdfs{neuron};
            %% preferred_angle and OSI DSI
            cond_preferred_angle = nan(1,8);
            cond_preferred_ind   = nan(1,8);
            osi                  = nan(1,8);
            dsi                  = nan(1,8);
            psth_locs            = nan(1,8);
            for ip=1:8
                Report=0;
                [cond_preferred_angle(ip),cond_preferred_ind(ip),osi(ip),dsi(ip)]=...
                    f_find_preferred6_cond(spike_data,stim_angle,angle_loc,onset_ind,offset_ind,ip,KS_dirs{day},neuron,Report);
                psth_locs(ip)=angle_loc(cond_preferred_ind(ip),ip);
            end
            psth_locs(7) = angle_loc(cond_preferred_ind(1),7);
            psth_locs(8) = angle_loc(cond_preferred_ind(2),8);
            cond_preferred_angle(7) = cond_preferred_angle(1);
            cond_preferred_angle(8) = cond_preferred_angle(2);
            OSI = [OSI;osi];
            DSI = [DSI;dsi];
            preferred_ang=[preferred_ang;cond_preferred_angle];
            %% pre_stimulus_baseline
            pre_stimulus_baseline=mean(mean(Pdf(psth_locs,1:onset_ind-1)));
            if pre_stimulus_baseline
                Pdf=(Pdf-pre_stimulus_baseline);%./pre_stimulus_baseline; %% normalizing psth by pre_stimulus_baseline
            end
            pre_stimulus_baseline_smooth=mean(mean(smoothed_pdf(psth_locs,1:onset_ind-1)));
            if pre_stimulus_baseline_smooth
                smoothed_pdf=(smoothed_pdf-pre_stimulus_baseline_smooth);%./pre_stimulus_baseline_smooth; %% normalizing psth by pre_stimulus_baseline
            end
            %% detecting channel
            [~,cell_chan]=min(abs(sp.ycoords-depth(neuron)));
            detecting_chan=[detecting_chan,cell_chan];
            %% collect all psths
            grand_psth_ncs    = [grand_psth_ncs;Pdf(psth_locs(1),:)];
            grand_psth_ldg    = [grand_psth_ldg;Pdf(psth_locs(2),:)];
            grand_psth_square = [grand_psth_square;Pdf(psth_locs(4),:)];
            opto_psth_ldg     = [opto_psth_ldg;Pdf(psth_locs(8),:)];
            opto_psth_ncs     = [opto_psth_ncs;Pdf(psth_locs(7),:)];
            %% conditional calculations like ave_activity, phase and latency
            for cond=1:length(psth_locs) %%8 condition (1:neon 2:phy 3:ring block (ctrl) 4:squar block 5:circle 6:phy+circle, ,'NCS+light','LDG+light) of preferred angle 
                ave_activity_cond(neuron,cond)  = mean(Pdf(psth_locs(cond),onset_ind:offset_ind));
                ave_sweep_cond(neuron,cond)  = mean(Pdf(psth_locs(cond),onset_ind:onset_ind+transient_nbins));
                ave_feed_cond(neuron,cond)  = mean(Pdf(psth_locs(cond),onset_ind+transient_nbins:offset_ind));                
                ave_activity_cond2(neuron,cond) = mean(Pdf2(psth_locs(cond),onset_ind:offset_ind));
                cond_smoothed_pdf = smoothed_pdf(psth_locs(cond),onset_ind+transient_nbins-1:offset_ind+transient_nbins);
                s1 =length(cond_smoothed_pdf); m  = s1 - mod(s1, bin_size);
                cond_smoothed_pdf=mean(reshape(cond_smoothed_pdf(1:m),bin_size,[]));
                %% statistical test for evoked cells
                test_vector=spike(spike(:,3)==psth_locs(cond),1:2);
                sum_pre_activity=[];sum_post_activity=[];
                trials=unique(test_vector(:,2));
                if isempty(trials)
                    p_values(neuron,cond)=nan;h_values(neuron,cond)=nan;
                else
                    for tr=1:length(trials)
                        trial=trials(tr);
                        sum_pre_activity  = [sum_pre_activity;sum(test_vector(test_vector(:,2)==trial,1)<0)];
                        sum_post_activity = [sum_post_activity;sum(test_vector(test_vector(:,2)==trial,1)<0.6 & test_vector(test_vector(:,2)==trial,1)>0)/2];
                    end
                    [p_values(neuron,cond),h_values(neuron,cond)]=ranksum(sum_pre_activity,sum_post_activity,'alpha', 0.001);
                end
                %% phase  and modulation calculation
                L=length(cond_smoothed_pdf);
                freq=FS*(0:(L/2))/L;%%0:FS/length(x):FS-1;
                [~,f1_ind]    = min(abs(freq-tf));
                [~,f2_ind]    = min(abs(freq-tf*2));
                f             = fft(cond_smoothed_pdf)/L;
                phase_angle_cond(neuron,cond)= angle(f(f1_ind));
                f1            = abs(f(f1_ind))*2;
                f2            = abs(f(f2_ind))*2;
                f0=f(1);
                first_comp_cond(neuron,cond) = f1;
                zero_comp_cond(neuron,cond)  = f0;
                f1f0_cond(neuron,cond)       = (f1-f0)/(f0+f1);
                f2f1f0_cond(neuron,cond)     = (f2+f1-f0)/(f0+f2+f1);
                f2f1_cond(neuron,cond)       = (f2-f1)/(f2+f1);
                f(1)=0; f=abs(f);f=f(1:floor(L/2)+1);f(2:end-1)=2*f(2:end-1);
                if cond==2
                    ff=fft(cond_smoothed_pdf-mean(cond_smoothed_pdf));ff=abs(ff/L);ff=ff(1:floor(L/2)+1);ff(2:end-1)=2*ff(2:end-1);
                    [~ ,domominat_loc]=max(ff);
                    if domominat_loc==f1_ind
                        dominant_f1(neuron)=1;
                    end
                end
                %% power spectrum plot
                if plot_power
                    subplot(4,2,cond_loc_plot(cond))
                    plot(freq,f)
                    hold on
                    plot(f1_ind-1,f1,'r*')
                    plot((f2_ind-1),f2,'g*')
                    title([cond_names{cond},' f1f0=',num2str(f1f0_cond(neuron,cond)),',f1f2=',num2str(f2f1_cond(neuron,cond))])
                    xlim([0,25])
                    xlabel('frequency (Hz)')
                    ylabel('power')
                end

                %% conditional response latency
                cond_RF = RF(psth_locs(cond),:);
                n=4; %number of merging bins
                cond_RF4=cond_RF;
                s1=length(cond_RF4);
                m  = s1 - mod(s1, n);
                cond_RF4=reshape(cond_RF4(1:m),n,[]);
                cond_RF4=sum(cond_RF4);
                teta    = f_resp_latency(cond_RF4,floor(onset_ind/n),0,floor(offset_ind/(n)));
                if ~isempty(teta)
                    latency_cond4(neuron,cond)=teta*n;
                end
                cond_RF2 = movsum(cond_RF,3);
                teta    = f_resp_latency(cond_RF2,onset_ind,0,offset_ind);
                if ~isempty(teta)
                    latency_cond2(neuron,cond)=teta;
                end
                teta=f_resp_latency_5p(smoothed_pdf(psth_locs(cond),:),onset_ind);
                latency_cond_5p(neuron,cond)=teta;
            end
            if plot_power
                saveas(PS,[data_folder,'/PS_',neuron_name,'.png'])
                close(PS)
            end

            if mean(ave_activity_cond(neuron,:))>=0
                response_type(neuron)=1; %%%"Excited";
            end
            %% plot raster for preferred angle
            if plot_raster
                f_neonopto_plot_pref_raster(KS_dirs{day},cond_preferred_angle,psth_locs,spike_data,neuron,...
                    latency_cond_5p(neuron,:),latency_cond2(neuron,:),bin_size+10,area_names{RF_cell.area(neuron)},...
                    session_name,f1f0_cond(neuron,:),cell_chan,p_values(neuron,:)) %#ok<*UNRCH>
            end
                           
            %% simplicity T-test
%             tag{dummy_all_neuron}
            phy_smoothed_pdf   = (abs(smoothed_pdf(psth_locs(2),onset_ind+transient_nbins-1:offset_ind)));
            xp_phy             = xp(onset_ind+transient_nbins-1:offset_ind);
            [~,h_simpl]        = f_simplicity_test(phy_smoothed_pdf, xp_phy, spike, psth_locs);
            h_simplicity       = [h_simplicity;h_simpl];

        end
        
        ave_act      = vertcat(ave_act,ave_activity_cond);     %ave_act is the mean value of reduced responce (pre-stim reduced)
        ave_act_sweep= vertcat(ave_act_sweep,ave_sweep_cond);
        ave_act_feed = vertcat(ave_act_feed,ave_feed_cond);
        ave_act2     = vertcat(ave_act2,ave_activity_cond2);   %ave_act2 is the mean value of raw Firing rate
        exited       = vertcat(exited,response_type);
        phase_angle  = vertcat(phase_angle,phase_angle_cond);
        first_comp   = vertcat(first_comp,first_comp_cond);
        zero_comp    = vertcat(zero_comp,zero_comp_cond);
        f1f0         = vertcat(f1f0,f1f0_cond);
        f2f1f0       = vertcat(f2f1f0,f2f1f0_cond);
        f2f1         = vertcat(f2f1,f2f1_cond);
        latency4     = vertcat(latency4,latency_cond4);
        latency2     = vertcat(latency2,latency_cond2);
        latency5p     = vertcat(latency5p,latency_cond_5p);
        
        F1_dominant  = vertcat(F1_dominant,dominant_f1);

        H_values     = vertcat(H_values,p_values);
        P_values     = vertcat(P_values,p_values);

        
        %%  surround_suppression
        SR_day=NaN(length(evoked_cids),1);
%         for isr=1:length(evoked_cids)
%             id_loc=ismember(surround_sup(1,:),evoked_cids(isr));
%             if sum(id_loc)~=0
%                 SR_day(isr)=surround_sup(2,id_loc);
%             end
%         end
        SR_day(SR_day<-1)=nan;
        S_S   = vertcat(S_S,SR_day);

    end
end

phase_diff               = rad2deg(angdiff(phase_angle(:,1),phase_angle(:,2)));
twopi_phase_diff         = rad2deg(wrapTo2Pi((phase_angle(:,1)-phase_angle(:,2))));
norm_rel_amp_f1          = (first_comp(:,1)-first_comp(:,2))./(first_comp(:,1)+first_comp(:,2));
ave_IGR2                 = ((ave_act2(:,1))-(ave_act2(:,2)))./((ave_act2(:,1))+(ave_act2(:,2)));   %ave_act2 is the mean value of raw Firing rate
ave_IGR1                 = (abs(ave_act(:,1))-abs(ave_act(:,2)))./(abs(ave_act(:,1))+abs(ave_act(:,2))); %ave_act is the mean value of reduced responce (pre-stim reduced)
ave_IGR1_sweep           = (abs(ave_act_sweep(:,1))-abs(ave_act_sweep(:,2)))./(abs(ave_act_sweep(:,1))+abs(ave_act_sweep(:,2))); %ave_act is the mean value of reduced responce (pre-stim reduced)
ave_IGR1_feed            = (abs(ave_act_feed(:,1))-abs(ave_act_feed(:,2)))./(abs(ave_act_feed(:,1))+abs(ave_act_feed(:,2))); %ave_act is the mean value of reduced responce (pre-stim reduced)


%intersection_point       = f_intersect_gauss2_for_hist_v2(TPlatency,1000,.4,.5,1); %inputs(vect,Nbin,aa,bb,report) /sugest an in terval for the intersection [aa,bb]
intersection_point       = 0.5021;
Icell_ind                = TPlatency<=intersection_point;
how_fast                 = 80; %(ms)
fast_phy                 = latency4(:,2) <=how_fast;
fast_neon                = latency4(:,1) <=how_fast;
fast_phy_neon            = logical(fast_neon.*fast_phy);
shifted_ind              = abs(phase_diff) >=120;

%% being evoked in different condition based on response latency
latency2_copy=latency2;
latency2(latency2>300)=nan;
% evoke_phy_ind  = ~isnan(latency2(:,2));
% evoke_neon_ind = ~isnan(latency2(:,1));
% both_evoked    = evoke_neon_ind & evoke_phy_ind;
% only_phy       = evoke_phy_ind & ~evoke_neon_ind;
% only_neon      = ~evoke_phy_ind & evoke_neon_ind;
% 
% excited_idx      = logical(exited);
% evoked_idx       = true(length(exited),1);
% LDG_excited_idx  = evoke_phy_ind & ave_act(:,2)>=0;
% both_excited_idx = both_evoked & ave_act(:,2)>=0;
% com_ldg_excited  = LDG_excited_idx & f1f0(:,2)<=0;
% 
% evoked_any    = (sum(~isnan(latency2)')~=0)';
% evoked_V1_ind = Cell_area==1 & ~isnan(latency2);
% evoked_lm_ind = Cell_area==2 & ~isnan(latency2);
% evoked_li_ind = Cell_area==3 & ~isnan(latency2);
% evoked_ll_ind = Cell_area==4 & ~isnan(latency2);
%% being evoked in different condition based on ranksum test
pv=0.005;
evoke_phy_ind  = P_values(:,2)<pv;
evoke_neon_ind = P_values(:,1)<pv;
both_evoked    = evoke_neon_ind & evoke_phy_ind;
only_phy       = evoke_phy_ind & ~evoke_neon_ind;
only_ncs      = ~evoke_phy_ind & evoke_neon_ind;

excited_idx      = logical(exited);
LDG_excited_idx  = evoke_phy_ind & ave_act(:,2)>=0;
both_excited_idx = both_evoked & ave_act(:,2)>=0;
com_ldg_excited  = LDG_excited_idx & f1f0(:,2)<=0;

evoked_any    = (sum(P_values(:,1:6)<pv,2)~=0);
evoked_V1_ind = Cell_area==1 & P_values<pv;
evoked_lm_ind = Cell_area==2 & P_values<pv;
evoked_li_ind = Cell_area==3 & P_values<pv;
evoked_ll_ind = Cell_area==4 & P_values<pv;
%% percentage of evoked 

total_v1_evoked = sum (evoked_any & Cell_area==1);
total_lm_evoked = sum (evoked_any & Cell_area==2);
total_li_evoked = sum (evoked_any & Cell_area==3);
total_ll_evoked = sum (evoked_any & Cell_area==4);

NCS_V1_evoked = sum(evoked_V1_ind(:,1))/total_v1_evoked;
NCS_lm_evoked = sum(evoked_lm_ind(:,1))/total_lm_evoked;
NCS_li_evoked = sum(evoked_li_ind(:,1))/total_li_evoked;
NCS_ll_evoked = sum(evoked_ll_ind(:,1))/total_ll_evoked;
NCS_evoked_percent=[NCS_V1_evoked NCS_lm_evoked NCS_li_evoked NCS_ll_evoked]*100;

ldg_V1_evoked = sum(evoked_V1_ind(:,2))/total_v1_evoked;
ldg_lm_evoked = sum(evoked_lm_ind(:,2))/total_lm_evoked;
ldg_li_evoked = sum(evoked_li_ind(:,2))/total_li_evoked;
ldg_ll_evoked = sum(evoked_ll_ind(:,2))/total_ll_evoked;
LDG_evoked_percent=[ldg_V1_evoked ldg_lm_evoked ldg_li_evoked ldg_ll_evoked]*100;

sqr_V1_evoked = sum(evoked_V1_ind(:,4))/total_v1_evoked;
sqr_lm_evoked = sum(evoked_lm_ind(:,4))/total_lm_evoked;
sqr_li_evoked = sum(evoked_li_ind(:,4))/total_li_evoked;
sqr_ll_evoked = sum(evoked_ll_ind(:,4))/total_ll_evoked;
SQR_evoked_percent=[sqr_V1_evoked sqr_lm_evoked sqr_li_evoked sqr_ll_evoked]*100;

crcl_V1_evoked = sum(evoked_V1_ind(:,5))/total_v1_evoked;
crcl_lm_evoked = sum(evoked_lm_ind(:,5))/total_lm_evoked;
crcl_li_evoked = sum(evoked_li_ind(:,5))/total_li_evoked;
crcl_ll_evoked = sum(evoked_ll_ind(:,5))/total_ll_evoked;
CRCL_evoked_percent = [crcl_V1_evoked crcl_lm_evoked crcl_li_evoked crcl_ll_evoked]*100;

nl_V1_evoked = sum(evoked_V1_ind(:,7)& ave_act(:,7)>0)/total_v1_evoked;
nl_lm_evoked = sum(evoked_lm_ind(:,7)& ave_act(:,7)>0)/total_lm_evoked;
nl_li_evoked = sum(evoked_li_ind(:,7)& ave_act(:,7)>0)/total_li_evoked;
nl_ll_evoked = sum(evoked_ll_ind(:,7)& ave_act(:,7)>0)/total_ll_evoked;
nl_evoked_percent = [nl_V1_evoked nl_lm_evoked nl_li_evoked nl_ll_evoked]*100;

pl_V1_evoked = sum(evoked_V1_ind(:,8)& ave_act(:,8)>0)/total_v1_evoked;
pl_lm_evoked = sum(evoked_lm_ind(:,8)& ave_act(:,8)>0)/total_lm_evoked;
pl_li_evoked = sum(evoked_li_ind(:,8)& ave_act(:,8)>0)/total_li_evoked;
pl_ll_evoked = sum(evoked_ll_ind(:,8)& ave_act(:,8)>0)/total_ll_evoked;
pl_evoked_percent = [pl_V1_evoked pl_lm_evoked pl_li_evoked pl_ll_evoked]*100;


evoked_percent = [NCS_evoked_percent; LDG_evoked_percent;SQR_evoked_percent;CRCL_evoked_percent;nl_evoked_percent;pl_evoked_percent];

%% edge selectives  
edge_selective = f2f1>=0 & f2f1f0>=0 & ~isnan(latency4);
sum(edge_selective)
neon_edge=edge_selective(:,1) & exited;

%% saving

save([data_folder,'/neon_opto_stat_data.mat'])

%% 

toc
 

