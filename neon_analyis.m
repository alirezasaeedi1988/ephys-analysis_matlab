%%%% neon color analysis
%%% extracting stimulus data and spike data and ploting raster, psth and
%%% tuning curve...
clear
tic
%% parameters
%%
data_folder= 'D:\matwork\Data\NEON_V1';
data_folder2='Z:\';
mouse_numbers =[145 153 165 185 194 335];

addpath(genpath('D:\matwork\KiloSort'))
addpath(genpath('D:\matwork\npy-matlab'))
addpath(genpath('D:\matwork\spikes'))
addpath(genpath('D:\matwork\NPMK'))
Nchan=32;

for mouse_number=mouse_numbers
%% getting forders for analysis
%%
  %%% Kilosort folders
mouse=['M',num2str(mouse_number)]
KS_dirs=cell(2,1);
folders=dir(data_folder);
j=1;
for i=1:length(folders)
    if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
        KS_dirs{j}=fullfile(data_folder,folders(i).name);
        j=j+1;
    end 
end

  %%% OpenEphys folders

OE_dirs=cell(2,1);
folders=dir(fullfile(data_folder2,mouse));
% folders=dir(data_folder2);
j=1;
for i=1:length(folders)
    if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
        OE_dirs{j}=fullfile(data_folder2,mouse,folders(i).name);
%         OE_dirs{j}=fullfile(data_folder2,folders(i).name);
        j=j+1;
    end 
end

%% extracting stimulus data (this part will exrtact data for all the sessions)
%%

% for day=1:length(KS_dirs) %%  loop for different days
%    f_stim_data(OE_dirs{day},KS_dirs{day});
% %    f_stim_data_blckrck(OE_dirs{day},KS_dirs{day+7},Nchan);
% end

%% extracting spike data (this part will exrtact data for all the sessions)

% for day=1:length(KS_dirs) %%  loop for different days
% %     f_extract_spike_data(KS_dirs{day})
%     f_extract_spike_data_poisson_test(KS_dirs{day})
% end

%% extracting spike Feature (this part will calculate feature for all the sessions and all neourons)
% for day=1:length(KS_dirs) %%  loop for different days
%     f_spikeFeature(KS_dirs{day})
% end
% toc
%% plotting Master raster and psth (this will be applied on all the sessions)
%%
% 
% for day=1:length(KS_dirs) %% loop for different days
%     f_master_raster_PSTH(KS_dirs{day})
% end


%% plotting conditional raster, psth and tuning curve

% for day=1:length(KS_dirs) %% loop for different days
%     inputfile_stm           = [KS_dirs{day},KS_dirs{day}(end-15:end),'_stimulus_time.mat'];
%     load(inputfile_stm)
%     dim=[];
%     for session_number=1:length(stimulus_schedule)
%         dim=[dim,min(size(stimulus_schedule{session_number}))]; %#ok<*AGROW>
%     end
%     if max(dim)==2
% %         f_neon_conditional_raster_psth(KS_dirs{day})
%     elseif max(dim)==3
%         f_neon_conditional_raster_psth_ctrl(KS_dirs{day})
%         f_neon_tuning_curve(KS_dirs{day})
%     end
% end
%% plotting tuning curves
for day=1:length(KS_dirs) %% loop for different days
%     f_plot_size_tuning(KS_dirs{day})
    f_neon_tuning_curve(KS_dirs{day})
end
 
%% SS
%  for day=1:length(KS_dirs) %% loop for different days
%      f_SS_size_tuning(KS_dirs{day})
%  end
end
toc
