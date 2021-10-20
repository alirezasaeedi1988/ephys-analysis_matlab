%%%% statistics
%% parameters
clear
tic
data_folder    = 'D:\matwork\data\higher_vision';
mouse_numbers  = [164 185 335 336 337 407];
session_name   = 'neon';
session_num    = 2; %for neon
plot_raster    = 0;  %%% plot raster for preferred angle?
plot_power     = 0;  %%% plot power spectrum
stim_angle     = [0 45 90 135 180 225 270 315];
stim_type      = [1 2];
neon_type      = [1 2 3];
Ncond          = length(stim_type)*length(stim_angle)*length(neon_type); 
angle_loc      = reshape(1:Ncond,length(stim_angle),length(stim_type)*length(neon_type));
area_names     = {'V1';'LM';'LI';'LL'};
ave_act        = [];
ave_act2       = [];
phase_angle    = [];
first_comp     = [];
zero_comp      = [];
f1f0           = [];
f2f1f0         = [];
f2f1           = [];
latency4       = [];
latency2       = [];
latency_phy    = [];
pref_ang       = [];
pref_ang_phy   = [];
pref_ang_neon  = [];
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
Cell_area      = [];
tag            = cell(2,1);
dummyAllneuron = 0;



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
        tf             =  stimulus_TemporalFreq{session_num};          % temporal frequency (cycle/s)
        cycles         =  stimulus_cycles{session_num};
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
        Cell_area          = vertcat(Cell_area,RF_cell.area);  %% area: 1 = V1, 2 = LM , 3 = LI , 4 = LL
        [~, onset_ind]     = min(abs(xp));
        [~,offset_ind]     = min(abs(xp-stim_len));
        [~,halfset_ind]    = min(abs(xp-stim_len/2));
        Nc                 = length(pdfs);       %number of clusters
        ave_activity_cond  = zeros(Nc,6);
        ave_activity_cond2 = zeros(Nc,6);
        phase_angle_cond   = zeros(Nc,6);
        first_comp_cond    = nan(Nc,6);
        zero_comp_cond     = nan(Nc,6);
        f1f0_cond          = nan(Nc,6);
        f2f1f0_cond        = nan(Nc,6);
        f2f1_cond          = nan(Nc,6);
        latency_cond4      = nan(Nc,6);
        latency4_phy       = nan(Nc,1);
        latency_cond2      = nan(Nc,6);
        response_type      = zeros(Nc,1);
        dominant_f1        = zeros(Nc,1);


        for neuron=1:Nc        %%%loop over clusters
            dummyAllneuron      = dummyAllneuron+1;
            tag{dummyAllneuron} = [mouse,tarikh,'id',num2str(spike_data.evoked_cids{neuron})];
            spike        = spikes{neuron};
            Pdf          = pdfs{neuron};
            Pdf2         = pdfs{neuron};
            RF           = RFs{neuron};
            smoothed_pdf = smoothed_pdfs{neuron};
            %%% preferred_angle and OSI DSI
            [preferred_angle,preferred_ind]=f_find_preferred6(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,'all');
            pref_ang=[pref_ang;preferred_angle];
            %psth_locs=angle_loc(preferred_ind,:);
            [preferred_angle_phy,~]=f_find_preferred6(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,'phy');
            pref_ang_phy=[pref_ang_phy;preferred_angle_phy];
            [preferred_angle_neon,~]=f_find_preferred6(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,'neon');
            pref_ang_neon=[pref_ang_neon;preferred_angle_neon];
            cond_preferred_angle = nan(1,6);
            cond_preferred_ind   = nan(1,6);
            osi                  = nan(1,6);
            dsi                  = nan(1,6);
            psth_locs            = nan(1,6);
            for ip=1:6
                [cond_preferred_angle(ip),cond_preferred_ind(ip),osi(ip),dsi(ip)]=f_find_preferred6_cond(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,ip);
                psth_locs(ip)=angle_loc(cond_preferred_ind(ip),ip);
            end
            OSI = [OSI;osi];
            DSI = [DSI;dsi];
            preferred_ang=[preferred_ang;cond_preferred_angle];
            %% ave_activity, phase and latency
            transient_time=150; %in ms to exclude the transient component
            transient_nbins=transient_time/(dx*1000);
            bin_size=10;
            FS=1/(dx*bin_size); % freq. of sampling in PSTH
            pre_stimulus_baseline=mean(mean(Pdf(psth_locs,1:onset_ind-1)));
            if pre_stimulus_baseline
                Pdf=(Pdf-pre_stimulus_baseline);%./pre_stimulus_baseline; %% normalizing psth by pre_stimulus_baseline
            end
            pre_stimulus_baseline_smooth=mean(mean(smoothed_pdf(psth_locs,1:onset_ind-1)));
            if pre_stimulus_baseline_smooth
                smoothed_pdf=(smoothed_pdf-pre_stimulus_baseline_smooth);%./pre_stimulus_baseline_smooth; %% normalizing psth by pre_stimulus_baseline
            end
            if plot_power
                neuron_name   = [mouse,'_D',KS_dirs{day}(num_loc(3)+1:num_loc(3)+2),'-',KS_dirs{day}(num_loc(2)+1:num_loc(2)+2),'_ID',num2str(spike_data.evoked_cids{neuron})];
                PS=figure('visible','off');
                cond_names = {'NCS', 'LDG','NCS+ring','NCS+squar','NCS+circle','LDG+circle'};
                cond_loc_plot = [1 2 6 5 3 4];
            end
 
            for cond=1:length(psth_locs) %%6 condition (1:neon 2:phy 3:ring block (ctrl) 4:squar block 5:circle 6:phy+circle) of preferred angle 
                ave_activity_cond(neuron,cond) = mean(Pdf(psth_locs(cond),onset_ind:offset_ind));
                ave_activity_cond2(neuron,cond) = mean(Pdf2(psth_locs(cond),onset_ind:offset_ind));
                cond_smoothed_pdf=smoothed_pdf(psth_locs(cond),onset_ind+transient_nbins-1:offset_ind+transient_nbins);
                s1=length(cond_smoothed_pdf);
                m  = s1 - mod(s1, bin_size);
                cond_smoothed_pdf=mean(reshape(cond_smoothed_pdf(1:m),bin_size,[]));
                
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
                    subplot(3,2,cond_loc_plot(cond))
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
                teta    = f_resp_latency(cond_RF2,floor(onset_ind),0,floor(offset_ind));
                if ~isempty(teta)
                    latency_cond2(neuron,cond)=teta;
                end
            end
            if plot_power
                saveas(PS,[data_folder,'/PS_',neuron_name,'.png'])
                close(PS)
            end
            %% phy response latency
            phy_RF =sum(RF(psth_locs(2:3),:));
            n=4;%number of merging bins
            s1=length(phy_RF);
            m  = s1 - mod(s1, n);
            phy_RF=reshape(phy_RF(1:m),n,[]);
            phy_RF=sum(phy_RF);
            teta    = f_resp_latency(phy_RF,floor(onset_ind/n),0,floor(offset_ind/(n*2)));
            if ~isempty(teta)
                latency4_phy(neuron)=teta*n;
            end
            if mean(ave_activity_cond(neuron,:))>=0
                response_type(neuron)=1; %%%"Excited";
            end
            %% plot raster for preferred angle
            if plot_raster
                f_neon6_plot_pref_raster(KS_dirs{day},preferred_angle,psth_locs,spike_data,neuron,latency_cond4(neuron,:),latency_cond2(neuron,:),bin_size,area_names{RF_cell.area(neuron)},'neon1',f1f0_cond(neuron,:)) %#ok<*UNRCH>
            end
                           
                
            %% simplicity T-test
%             tag{dummy_all_neuron}
            phy_smoothed_pdf   = (abs(smoothed_pdf(psth_locs(2),onset_ind+transient_nbins-1:offset_ind)));
            xp_phy             = xp(onset_ind+transient_nbins-1:offset_ind);
            [~,h_simpl] = f_simplicity_test(phy_smoothed_pdf, xp_phy, spike, psth_locs);
            h_simplicity       = [h_simplicity;h_simpl];

        end
        
        ave_act      = vertcat(ave_act,ave_activity_cond);     %ave_act is the mean value of reduced responce (pre-stim reduced)
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
        latency_phy  = vertcat(latency_phy,latency4_phy);
        F1_dominant  = vertcat(F1_dominant,dominant_f1);

        
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


%intersection_point       = f_intersect_gauss2_for_hist_v2(TPlatency,1000,.4,.5,1); %inputs(vect,Nbin,aa,bb,report) /sugest an in terval for the intersection [aa,bb]
intersection_point       = 0.5021;
Icell_ind                = TPlatency<=intersection_point;
how_fast                 = 80; %(ms)
fast_phy                 = latency4(:,2) <=how_fast;
fast_neon                = latency4(:,1) <=how_fast;
fast_phy_neon            = logical(fast_neon.*fast_phy);
shifted_ind              = abs(phase_diff) >=120;

%% neurons type based on being evoked in different condition
evoke_phy_ind  = ~isnan(latency2(:,2));
evoke_neon_ind = ~isnan(latency2(:,1));
both_evoked    = evoke_neon_ind & evoke_phy_ind;
only_phy       = evoke_phy_ind & ~evoke_neon_ind;
only_neon      = ~evoke_phy_ind & evoke_neon_ind;

excited_idx      = logical(exited);
evoked_idx       = true(length(exited),1);
LDG_excited_idx  = evoke_phy_ind & ave_act(:,2)>=0;
both_excited_idx = both_evoked & ave_act(:,2)>=0;
com_ldg_excited  = LDG_excited_idx & f1f0(:,2)<=0;

evoked_any    = (sum(~isnan(latency2)')~=0)';
evoked_V1_ind = Cell_area==1 & ~isnan(latency2);
evoked_lm_ind = Cell_area==2 & ~isnan(latency2);
evoked_li_ind = Cell_area==3 & ~isnan(latency2);
evoked_ll_ind = Cell_area==4 & ~isnan(latency2);

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

evoked_percent = [NCS_evoked_percent; LDG_evoked_percent;SQR_evoked_percent;CRCL_evoked_percent];

%% edge selectives  
edge_selective = f2f1>=0 & f2f1f0>=0 & ~isnan(latency4);
sum(edge_selective)
neon_edge=edge_selective(:,1) & exited;

%% saving

stat_data.ave_act        = ave_act;
stat_data.phase_angle    = phase_angle;
stat_data.latency        = latency4;
stat_data.f1f0           = f1f0;
stat_data.f2f1f0         = f2f1f0;
stat_data.f2f1           = f2f1;
stat_data.pref_ang       = pref_ang;
stat_data.cid            = cid;
stat_data.tag            = tag;
stat_data.ave_FR         = mean(ave_act,2);
stat_data.exited         = exited;
stat_data.ave_IGR2       = ave_IGR2;
stat_data.ave_IGR1       = ave_IGR1;
stat_data.ave_latency    = ave_latency./1000; 
stat_data.fast_phy_neon  = fast_phy_neon;
stat_data.Icell          = Icell_ind;
stat_data.shifted        = shifted_ind;
stat_data.log_phase_diff             = log(abs(phase_diff));
stat_data.phase_diff                 = phase_diff;
stat_data.twopi_phase_diff           = twopi_phase_diff;
stat_data.h_simplicity               = h_simplicity;
stat_data.norm_rel_amp_f1            = norm_rel_amp_f1;
stat_data.surround_suppression       = S_S;
stat_data.TPlatency                  = TPlatency;
stat_data.Depth                      = Depths;
stat_data.PTratio                    = PTratio;
stat_data.F1_dominant                = logical(F1_dominant);

save([data_folder,'/neon_stat_data.mat'],'stat_data')

%% 


toc

