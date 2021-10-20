function f_stim_data_blckrck(root,myKsDir,Nchan)
%% compute stimulus data for merged data recorded by black rock

addpath(genpath('D:\matwork\NPMK'))

% root                 = 'y:\MPI_EphysData_2020\higher_vision\M164_left\M164_2020_02_12'; %  where NSX files are	
% myKsDir              = 'D:\matwork\Data\higher_vision\M164_2020_02_12';   % destination     
outputfile           = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat']; % will be created	

fs                   = 30000;   %sampling frequency
% Nchan= 64;
% sessions_mat         = {'reti','neon','Lsiz','Usiz','2ret'};
sessions = dir(fullfile(root, 'da*.ns6'));
Nsessions           = length(sessions);

if Nsessions==3 %#ok<*USENS>
    sessions_mat={'reti';'size';'neon'};
elseif Nsessions==4 
    if Nchan==32
        sessions_mat={'reti';'size';'neon_c';'neon_f'}; %% for v1 recording
    elseif Nchan==64
        sessions_mat = {'reti','neon_f','Lsiz','Usiz'}; %% for lateral recording
    end
    
elseif Nsessions==5
    if Nchan==32
        sessions_mat = {'reti','size','neon_c','neon_f','2ret'}; %% for v1 recording
    elseif Nchan==64
        sessions_mat = {'reti','neon_f','Lsiz','Usiz','2ret'}; %% for lateral recording
    end
end









if ~exist(myKsDir, 'dir')
    mkdir(myKsDir)
end

ITI                  = 0.2;  %inter trial interval in second
ver_thr              = 0.5; %vertical threshold (as a factor of std)
hor_thr              = (ITI*fs);

%%%%%%%%%%%%%%%%%
stimulus_time         = cell(1,Nsessions);
stimulus_point        = cell(1,Nsessions);
stimulus_length       = cell(1,Nsessions);
stimulus_schedule     = cell(1,Nsessions);
stimulus_cycles       = cell(1,Nsessions);
stimulus_TemporalFreq = cell(1,Nsessions);
stimulus_angles       = cell(1,Nsessions);
stimulus_sizes        = cell(1,Nsessions);
% data                  = cell(1,Nsessions);
root_logfiles=dir(fullfile(root, '*.mat'));

end_of_session=0; %dummy for final time in each session
for session=1:Nsessions
    if session==1
        ITI       = .09;
        ver_thr   = 0.9; %vertical threshold (as a factor of std)
        hor_thr   = (ITI*fs);
    end
    NSX_fname       = fullfile(root,sessions(session).name);
    if exist(NSX_fname,'file')
        openNSx('report','read',NSX_fname);%,'t:3:10^4');
        d   = NS6.Data(Nchan+1,:);
        d=double(d);
        d=smooth(d,7);            %smooth data
        d=(d-mean(d))/std(d);     %Z_score data
        ind=find(d>ver_thr);             %apply threshold for value of data as a factor of std
        dif=diff(ind);
        dif(dif<hor_thr)=0;
        [~,locs] = findpeaks(dif);
        Sts=ind(locs+1);Sts=Sts';    %starting time of stimulus
        Sts=[ind(1),Sts]; %#ok<AGROW>
        stimulus_point{session}=Sts+end_of_session;
        stimulus_time{session}=(Sts+end_of_session)./fs;
        Fts=ind(locs); Fts=Fts';                      %finishing time of trial
        stimulus_length{session}=Fts-Sts(1:end-1);          % stimulus length
        
        end_of_session=length(d)+end_of_session;
                    
        for k=1:length(root_logfiles)
            if strcmp(root_logfiles(k).name(1:length(sessions_mat{session})),sessions_mat{session}) 
                matfile_name=root_logfiles(k).name
            end
        end
        load(fullfile(root,matfile_name))
        if exist('trial_schedule','var')
            stimulus_schedule{session}=trial_schedule;
        end
        clear trial_schedule
        if exist('temporal_period','var')
            stimulus_cycles{session}=stim_duration/temporal_period;
            stimulus_TemporalFreq{session}=60/temporal_period;
            clear temporal_period
        elseif exist('numFrames','var')
            stimulus_cycles{session}=stim_duration/numFrames;
            stimulus_TemporalFreq{session}=60/numFrames;
            clear numFrames
            
        else
            stimulus_cycles{session}=0;
            stimulus_TemporalFreq{session}=0;
        end
    end
    if exist('stim_angles_in_degrees','var')
        stimulus_angles{session} = stim_angles_in_degrees;
    end
    if exist('stim_sizes_in_degrees','var')
        stimulus_sizes{session}  = stim_sizes_in_degrees;
    end
    clear stim_sizes_in_degrees stim_angles_in_degrees
end


save(outputfile,'stimulus_cycles','stimulus_point','stimulus_time',...
    'stimulus_length','stimulus_schedule','stimulus_TemporalFreq','stimulus_sizes','stimulus_angles')
                                                                                                             