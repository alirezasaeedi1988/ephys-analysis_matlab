function f_extract_spike_data_poisson_test(myKsDir)
%%% extracting spike times and...
%%% loading files

inputfile_stm = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)
spike_times   = readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster = readNPY([myKsDir ,'\spike_clusters.npy']);
[cids, cgs]   = readClusterGroupsCSV([myKsDir ,'\cluster_group.tsv']);

if length(stimulus_point)==3 %#ok<*USENS>
    session_numbers=2:3;
    session_names={'RF';'size';'neon'};
elseif length(stimulus_point)==4 
    session_numbers=2:3;
    session_names={'RF';'size';'neon';'MAE'};
elseif length(stimulus_point)==5
    session_numbers=2:3;
    session_names={'RF';'size';'neon';'size2';'MAE'};
elseif length(stimulus_point)==6
    session_numbers=2:5;
    session_names={'RF';'size';'contextual';'size2';'neon';'MAE'};
end

 for session_number=3%session_numbers
    
    session_name=session_names{session_number};

    fs  =30000;                                          %sampling frequency
    trial_schedule=stimulus_schedule{session_number};
    Nt=min([length(trial_schedule),length(stimulus_time{session_number})-2]);          %#ok<*IDISVAR>
    pst=0.3;                                             % pre-stimulus time sec
    dx= 1e-3;                                            % bin size
    stim_len=mean(stimulus_length{session_number})/fs;
    %%%
    Spt=floor(stimulus_point{session_number}-(pst*fs));      %starting point of trials in the data
    nsample=floor(mean(diff(Spt(2:end-1))));                          %number of samples for each trials
    Sps=stimulus_point{session_number};                      %starting point of stimulus
    
    %%%

    if strcmp(session_name,'neon')
        pst=0.3;                                                 % pre-stimulus time sec
        Spt=floor(stimulus_point{session_number}-(pst*fs));      %starting point of trials in the data
        nsample=floor(mean(diff(Spt)));                          %number of samples for each trials
        trial_schedule(Nt+1:end,:)=[];
        if min(size(trial_schedule))==2
            realcolor_loc=trial_schedule(:,1)==2;
            trial_schedule(realcolor_loc,2)=trial_schedule(realcolor_loc,2).*10;
            trial_schedule=trial_schedule(:,2);
        else
            phys1_loc=(trial_schedule(:,1)==1 & trial_schedule(:,2)==2);
            phys2_loc=(trial_schedule(:,1)==2 & trial_schedule(:,2)==2);
            ctrl_loc=(trial_schedule(:,1)==2 & trial_schedule(:,2)==1);
            
            trial_schedule(phys1_loc,3)=trial_schedule(phys1_loc,3).*10;
            trial_schedule(phys2_loc,3)=trial_schedule(phys2_loc,3).*100;
            trial_schedule(ctrl_loc,3)=trial_schedule(ctrl_loc,3).*1000;
            
            trial_schedule=trial_schedule(:,3);
        end
    elseif strcmp(session_name,'contextual')
        number_of_discarding_trials=length(trial_schedule)-length(Sps);
        trial_schedule(1:number_of_discarding_trials)=[]; 
    end
    condLabels = unique(trial_schedule);
    %%%
    xp=-pst:dx:(nsample/fs)-pst;
    [~, onset_ind] = min(abs(xp));
    [~,offset_ind] = min(abs(xp-stim_len));
    transient_time=150; %in ms to exclude the transient component
    transient_nbins=floor(transient_time/(dx*1000));
    %%%
    good_cid= cgs==2;
    good_cid=cids(good_cid);
    spikes=cell(1,1);
    pdfs=cell(1,1);
    RFs=cell(1,1);    
    smoothed_pdf=cell(1,1);
    
    evoked_cids=cell(1,1);
    average_latency=cell(1,2);
    trial_per_cond=cell(1,length(good_cid));
   
    e=1;
    for neuron=1:length(good_cid) %%%loop over clusters, to extract spikes and pdf ...
        RF=zeros(length(condLabels),length(xp));
        hist0=zeros(length(condLabels),length(xp));
        smoothed_hist=zeros(length(condLabels),length(xp));
        cid= spike_cluster==good_cid(neuron);
        Cluster=spike_times(cid);
        Cluster=Cluster';
        Cluster = double(Cluster);
        con_Ntrial=zeros(1,length(condLabels));
        spike=[];
        for cond=1:length(condLabels)%%% loop over all condition for selected neuron
            ind_trGroup= trial_schedule(1:Nt)==condLabels(cond);
            cond_spt=Spt(ind_trGroup);
            cond_sps=Sps(ind_trGroup);
            con_Ntrial(cond)=length(cond_spt);
            Trial=cell(con_Ntrial(cond),1);
            
            for trial_number=1:con_Ntrial(cond) %% loop over all trails in a condition for alignment spike times
                ind_tr= cond_spt(trial_number)<=Cluster & Cluster<cond_spt(trial_number)+nsample;
                Tri=Cluster(ind_tr);
                if ~isempty(Tri)
                    %%% align
                    Tri=(Tri-cond_sps(trial_number))./(fs);                
                    spike=[spike;Tri' , ones(length(Tri),1)*trial_number , ones(length(Tri),1)*cond]; %#ok<AGROW>
                    Trial{trial_number}=Tri;

                end
            end
            %%%PSTH calculation
            P=zeros(1,length(xp));
            for trial_number=1:con_Ntrial(cond) %% loop over all trails in a condition for psth
                if ~isempty(Trial{trial_number})
                    for jj=1:length(Trial{trial_number})
                        n=floor(((Trial{trial_number}(jj)-xp(1))/dx))+1;
                        P(n)=P(n)+1;
                    end
                end
            end
            PP=P./(con_Ntrial(cond)*dx); %% normalizing psth with number of trials and bin size
            RF(cond,:)=P;  %relative frequency
            hist0(cond,:)=PP;
            smoothed_hist(cond,:)=smoothdata(PP,'gaussian',15);%%% In 'smoothdata' the SD is fixed to 1/5th of window width
            
        end
        %%% poisson test for evoked and non-evoked neuron
        y    =  sum(RF);
        if strcmp(session_name,'size') || strcmp(session_name,'size2')
             y    = movsum(y,2);
        end
        teta=f_resp_latency(y,onset_ind,transient_nbins,offset_ind);
        if ~isempty(teta)
            average_latency{e}=teta;
            spikes{e}=spike;
            trial_per_cond{e}=con_Ntrial;
            pdfs{e}=hist0;
            RFs{e}=RF;
            smoothed_pdf{e}=smoothed_hist;
            evoked_cids{e}=good_cid(neuron);
            
            e=e+1;
        end

    end
    if strcmp(session_name,'neon')
        [sr, sr30,pref_size]=f_surround_suppression(myKsDir);
        spike_data.surround_suppression=sr;
        spike_data.surround_suppression30=sr30;
        spike_data.pref_size=pref_size;
    end
     
    spike_data.spikes          = spikes;
    spike_data.RFs             = RFs;
    spike_data.pdfs            = pdfs;
    spike_data.smoothed_pdf    = smoothed_pdf;
    spike_data.pst             = pst;
    spike_data.dx              = dx;
    spike_data.xp              = xp;
    spike_data.cluster_id      = good_cid;
    spike_data.trial_per_cond  = trial_per_cond;
    spike_data.stim_len        = stim_len;
    spike_data.evoked_cids     = evoked_cids;
    spike_data.average_latency = average_latency;

    save([myKsDir,'\',session_name,'_spike_data.mat'],'spike_data')
    
end