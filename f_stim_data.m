function f_stim_data(root,myKsDir)
%% compute stimulus data for merged data

% root                 = 'Y:\MPI_EphysData_2019\M146\2019_02_23';


outputfile           = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat']; % will be created		
Nsessions            = 1;      %%%number of sessions that needs to merged
fs                   = 30000;   %sampling frequency
file_name            ='100_ADC5.continuous';
% file_name            ='116_ADC4.continuous';

   
if ~exist(myKsDir, 'dir')
    mkdir(myKsDir)
end
ITI                  =30;  %inter trial interval in frame
%%%%%
ver_thr              = .7; %vertical threshold (as a factor of std)
hor_thr              =(ITI*fs)/(60*2);
sessions_strct=dir(fullfile(root));
sessions_dirindx=zeros(1,Nsessions);
jj=1;
for k=1:length(sessions_strct)
    if length(sessions_strct(k).name)==19 && strcmp(sessions_strct(k).name(1:3),'201')
        %keyboard
        sessions_dirindx(jj)= k;
        jj=jj+1;
    end
end
%%%%%%%%%%%%%%%%%
stimulus_time         = cell(1,Nsessions);
stimulus_point        = cell(1,Nsessions);
stimulus_length       = cell(1,Nsessions);
stimulus_schedule     = cell(1,Nsessions);
stimulus_cycles       = cell(1,Nsessions);
stimulus_TemporalFreq = cell(1,Nsessions);
stimulus_angles       = cell(1,Nsessions);
stimulus_sizes        = cell(1,Nsessions);
data                  = cell(1,Nsessions);
iii=1;
end_of_session=0; %dummy for final time in each session
for dirindx=sessions_dirindx  % go through different sessions
    if iii==1
        ITI       = .09;
        ver_thr   = 0.9; %vertical threshold (as a factor of std)
        hor_thr   = (ITI*fs);
    end
    stim_angles_in_degrees = 0;
    stim_sizes_in_degrees  = 0;
    input_root=[root,'\',sessions_strct(dirindx).name];
%     [d,~, info] = load_open_ephys_data([input_root,'\',file_name]);
    d = load_openEphys(input_root,file_name);
    d=double(d);
    d=smooth(d,7);            %smooth data
    d=(d-mean(d))/std(d);     %Z_score data
    data{iii}=d;
    ind=find(d>ver_thr);             %apply threshold for value of data as a factor of std
    dif=diff(ind);
    dif(dif<hor_thr)=0;
    [~,locs] = findpeaks(dif);
    Sts=ind(locs+1);Sts=Sts';    %starting time of stimulus
    Sts=[ind(1),Sts]; %#ok<AGROW>
    stimulus_point{iii}=Sts+end_of_session;
    stimulus_time{iii}=(Sts+end_of_session)./fs;
    Fts=ind(locs); Fts=Fts';                      %finishing time of trial
    stimulus_length{iii}=Fts-Sts(1:end-1);          % stimulus length
    % x=1:length(d);  sts1 = x(diff(d) > 7);

    end_of_session=length(d)+end_of_session;
    
    ff=dir(fullfile([root,'\',sessions_strct(dirindx).name]));
    for k=1:length(ff)
        if length(ff(k).name)>=10 && strcmp(ff(k).name(end-3:end),'.mat')
            matfile=ff(k).name;
        end
    end
    load([root,'\',sessions_strct(dirindx).name,'\',matfile])
    stimulus_schedule{iii}=trial_schedule;
    if exist('temporal_period','var')
        stimulus_cycles{iii}=stim_duration/temporal_period;
        stimulus_TemporalFreq{iii}=(60/temporal_period);
        clear temporal_period
    elseif exist('numFrames','var')
        stimulus_cycles{iii}=stim_duration/numFrames;
        stimulus_TemporalFreq{iii}=(60/numFrames);
        clear numFrames

    else
        stimulus_cycles{iii}=0;
        stimulus_TemporalFreq{iii}=0;
    end
    stimulus_angles{iii} = stim_angles_in_degrees;
    stimulus_sizes{iii}  = stim_sizes_in_degrees;
    if exist('distance_from_screen','var')
        SCR.dist_in_cm=distance_from_screen;
        SCR.width_in_cm = screen_width_in_cm;
        SCR.height_in_cm = screen_height_in_cm ;
        stim_shift=stim_shift_factor(1:2);
    end
    iii=iii+1;
end


save(outputfile,'stimulus_cycles','stimulus_point','stimulus_time',...
    'stimulus_length','stimulus_schedule','stimulus_TemporalFreq',...
    'stimulus_sizes','stimulus_angles','SCR','stim_shift')
                                                                                                             