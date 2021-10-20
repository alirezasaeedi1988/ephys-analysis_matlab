function f_RF_single_unit_contextual(myKsDir,session_number,session_name,Plot)


% myKsDir     = 'D:\matwork\Data\higher_vision\M164_2020_02_12';
% session_number = 1;
% session_name='RF';
SD=2; %masking treshold
outputDir           = fullfile(myKsDir,[session_name,'_SD',num2str(SD)]);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end


%% load Stimulus data
stimTimeName = dir(fullfile(myKsDir, '*time.mat'));
if ~isempty(stimTimeName)
    stimTimeName = fullfile(myKsDir, stimTimeName(1).name); load(stimTimeName)
else
    error('there is no Stimulus-time in the directory');
end


%% laod kilosort results
spike_times    = readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster  = readNPY([myKsDir ,'\spike_clusters.npy']);
% [cids, cgs]    = readClusterGroupsCSV([myKsDir ,'\cluster_group.tsv']);

%% laod contextual data
neonDataName = dir(fullfile(myKsDir, 'contextual*.mat'));
if ~isempty(neonDataName)
    neonDataName = fullfile(myKsDir, neonDataName(1).name); load(neonDataName)
else
    error('there is no contextual_data in the directory');
end
evoked_cids = cell2mat(spike_data.evoked_cids);




%% other params
res_scale      = 4;
min_lim_size   = 5;   % degree      
max_lim_size   = 25;  % degree

Stim_onset_ind = floor(stimulus_point{session_number});      
trial_schedule = stimulus_schedule{session_number};
stim_size_deg  = stimulus_sizes{3}; %stim size for contextual
monitor_dist   = SCR.dist_in_cm;  %cm
monitor_X      = SCR.width_in_cm;%cm
monitor_Y      = SCR.height_in_cm; %cm
stim_loc_X     = (.5+stim_shift(1))* monitor_X; %cm
stim_loc_Y     = (.5+stim_shift(2))* monitor_Y; %cm
cm_per_degree  = tand(1/2) * monitor_dist; %% this is not exact but it is how we defined it for stimulus peresentation
stim_radius    = stim_size_deg*cm_per_degree;
nsample        = floor((mean(diff(Stim_onset_ind(2:end-1))))/2);      %number of samples for each trials
Nt             = min([length(trial_schedule),length(Stim_onset_ind)-2]);
grid_cols      = max(unique(trial_schedule(:,1)));
grid_rows      = max(unique(trial_schedule(:,2)));
col_width      = monitor_X/(grid_cols*res_scale);  % cm
row_width      = monitor_Y/(grid_rows*res_scale);  % cm
col_angle      = 2*atand(col_width/(2*monitor_dist));
row_angle      = 2*atand(row_width/(2*monitor_dist));



%% trial_counter
trial_counter  = zeros(grid_rows,grid_cols);
for trial_idx  = 1:numel(Stim_onset_ind)
    trial_counter(trial_schedule(trial_idx,2), trial_schedule(trial_idx,1)) = ...
        trial_counter(trial_schedule(trial_idx,2), trial_schedule(trial_idx,1)) ...
        + 1;
end


%% stimulus circle
th = 0:pi/50:2*pi;
x_circle    = stim_radius * cos(th) + stim_loc_X;
y_circle    = stim_radius * sin(th) + stim_loc_Y;
stim_circle = polyshape(x_circle,y_circle); 
stim_area = area(stim_circle);

%% define empty matrices
RF_size   = zeros(length(evoked_cids),1);
RSS       = zeros(length(evoked_cids),1); % goodness of fitting (1-(Residual_sum_of_squares[fit]/Residual_sum_of_squares[tot]))
RF_HWHM_X = zeros(length(evoked_cids),1);
RF_HWHM_Y = zeros(length(evoked_cids),1);
RF_loc    = zeros(length(evoked_cids),2);
X_Bound   = {};%zeros(length(evoked_cids),1);
Y_Bound   = {};%zeros(length(evoked_cids),1);
RF_intrsct_ratio= nan(length(evoked_cids),1);
RF_area   = nan(length(evoked_cids),1);
fs=30000;
 for neuron=1:length(evoked_cids) %% loop over clusters, spike time for selected cluster is extracted

     firnig_rate = zeros(grid_rows, grid_cols);
     cid         = spike_cluster==evoked_cids(neuron);
     Neuron_spike_times = spike_times(cid);
     Neuron_spike_times = Neuron_spike_times';
     Neuron_spike_times = double(Neuron_spike_times);
     for trial_number = 1:Nt
         Trial_spike_ind   = Stim_onset_ind(trial_number)<=Neuron_spike_times & ...
              Neuron_spike_times<Stim_onset_ind(trial_number)+nsample;
         Trial_spike_times =  Neuron_spike_times(Trial_spike_ind);                    
         firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))=...
             firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))+ ...
             length(Trial_spike_times)/((trial_counter(trial_schedule(trial_number,2), trial_schedule(trial_number,1)))*(nsample/fs));

     end
     
     S = imresize(firnig_rate,res_scale);
     if ~std(S(:))==0
         S=(S-mean(S(:)))/std(S(:));S(abs(S)<= SD)=0; S=abs(S);
         
         %% get the boundaries of RF
         BW = imbinarize(S);
         BW_filled = imfill(BW,'holes');
         boundaries2 = bwboundaries(BW_filled,4);
         allLengths = cellfun(@length,boundaries2);
         [~,maxLoc] = max(allLengths);
         Boundaries = boundaries2{maxLoc};
         Boundaries(:,2)    = Boundaries(:,2)*col_width; %cm
         Boundaries(:,1)    = Boundaries(:,1)*row_width; %cm
         X_Bound{neuron}    = Boundaries(:,2); %cm
         Y_Bound{neuron}    = Boundaries(:,1); %cm
         
         %% area intersection
         RF_poly=polyshape(X_Bound{neuron},Y_Bound{neuron} );
         RF_area(neuron) = area(RF_poly);  % cm^2
         
         
         polyout = intersect(RF_poly,stim_circle);
         area_int= area(polyout);
         
         RF_intrsct_ratio(neuron) = area_int/stim_area*100;
         
         %% 2D guassian fitting and finding the FR properties
         [A,RSS(neuron)]=f_Gauss2D_contextual(S,Plot,0,outputDir,evoked_cids(neuron),row_angle ,col_angle,col_width,row_width,Boundaries,stim_circle);
         %[Amp,x0,wx,y0,wy,theta]=f_Gauss2D(S,Plot,FitOrientation,outputDir,cell_id,row_angle ,col_angle);
         HWHM_X       = sqrt(2*log(2))*A(3);        % half width at half maximum
         HWHM_Y       = sqrt(2*log(2))*A(5);        % half width at half maximum
         HWHM_Y_angle = HWHM_Y*row_angle;
         HWHM_X_angle = HWHM_X*col_angle;
         RF_size(neuron)    = (HWHM_Y_angle+HWHM_X_angle)/2; %angle (degree)
         %RF_size(neuron)    = pi*HWHM_X*HWHM_Y;
         RF_loc(neuron,:)   = [A(2)*col_width, A(4)*row_width]; %cm
         RF_HWHM_X(neuron)  = HWHM_X*col_width; %cm
         RF_HWHM_Y(neuron)  = HWHM_Y*row_width; %cm
 
     end
 end
RF_size(RF_size<min_lim_size) = nan;
RF_size(RF_size>max_lim_size) = nan;
% good_rf = ~isnan(RF_size);
RF_loc_shifted      = RF_loc;
RF_loc_shifted(:,1) = RF_loc_shifted(:,1)- grid_cols;
RF_loc_shifted(:,2) = RF_loc_shifted(:,2)- grid_rows;

RF_data.intrsct_ratio= RF_intrsct_ratio;
RF_data.area         = RF_area;
RF_data.HWHM_X       = RF_HWHM_X; 
RF_data.HWHM_Y       = RF_HWHM_Y;
RF_data.X_Bound      = X_Bound;
RF_data.Y_Bound      = Y_Bound;
RF_data.RF_size      = RF_size;
RF_data.RF_loc       = RF_loc;
RF_data.RF_loc_shift = RF_loc_shifted;
RF_data.RSS          = RSS;   
RF_data.res_scale    = res_scale;
RF_data.session_name = session_name;
RF_data.evoked_cids  = evoked_cids;
RF_data.monitor_size = [monitor_X, monitor_Y]; %%cm
RF_data.stim_area    = stim_area;
RF_data.stim_circle  = stim_circle;
RF_data.stim_radius  = stim_radius;
RF_data.stim_loc_X   = stim_loc_X;
RF_data.stim_loc_Y   =stim_loc_Y;
save(fullfile(myKsDir,[session_name,'_data.mat']),'RF_data')

