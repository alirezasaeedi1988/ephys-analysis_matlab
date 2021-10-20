%%%% neon color analysis
%%% extracting stimulus data and spike data and ploting raster, psth and
%%% tuning curve...
clear
tic
addpath(genpath('D:\matwork\KiloSort'))
addpath(genpath('D:\matwork\npy-matlab'))
%% parameters
%%
% data_folder  = 'D:\matwork\Data\neon_opto_laser\New';
% data_folder2 = 'Y:\MPI_EphysData_2021\ddot_laser';
% mouse_numbers = [417   419   420   462   463   471];
data_folder  = 'D:\matwork\data\higher_vision';
data_folder2 = 'Y:\MPI_EphysData_2020\higher_vision';
mouse_numbers = [164   185   335   336   337   405   407];        %%
Nchan=64;
addpath(genpath('D:\matwork\KiloSort'))
addpath(genpath('D:\matwork\npy-matlab'))
addpath(genpath('D:\matwork\spikes'))
addpath(genpath('D:\matwork\NPMK'))

for mouse_number=mouse_numbers
    %% getting forders for analysis
    %%
    %%% Kilosort folders
    mouse=['M',num2str(mouse_number)]
    KS_dirs=cell(1,1);
    folders=dir(data_folder);
    j=1;
    for i=1:length(folders)
        if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
            KS_dirs{j,1}=fullfile(data_folder,folders(i).name);
            j=j+1;
        end
    end
    
    %%% blachrock folders
    j=1;
    BL_dirs=cell(1,1);
    folders=dir(data_folder2);
    subfolders=[];
    for i=1:length(folders)
        if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
            subfolders=[subfolders;dir(fullfile(data_folder2,folders(i).name))];

        end
    end
    T = struct2table(subfolders); % convert the struct array to a table
    sortedT = sortrows(T, 'name'); % sort the table by 'name'
    sorted_subfolders = table2struct(sortedT);
    for ii=1:length(sorted_subfolders)
        if length(sorted_subfolders(ii).name)>4 && strcmp(sorted_subfolders(ii).name(1:4),mouse)
            BL_dirs{j,1}=fullfile(sorted_subfolders(ii).folder,sorted_subfolders(ii).name);
            j=j+1;
        end
    end
    
%     for i=1:length(folders)
%         if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
%             BL_dirs{j,1}=fullfile(data_folder2,folders(i).name);
%             j=j+1;
%         end
%     end
    
    %% extracting stimulus data (this part will exrtact data for all the sessions)
    %%
    
%     for day=9:length(KS_dirs) %%  loop for different days
%     %    f_stim_data(BL_dirs{day},KS_dirs{day});
%        f_stim_data_blckrck(BL_dirs{day},KS_dirs{day},Nchan);
%     end
%     
    %% RF for channel using MUA
%     for day=1:length(KS_dirs) %% loop for different days
%         TF=dir(fullfile(KS_dirs{day},'channel_RF_data.mat'));
%         if isempty(TF)
%             f_RF_cahnnel(KS_dirs{day},BL_dirs{day})
%         end
%     end
%     close all
    
        %% extracting spike Feature (this part will calculate feature for all the sessions and all neourons)
%     for day=1:length(KS_dirs) %%  loop for different days
%         TF=dir(fullfile(KS_dirs{day},'spike_Feature.mat'));
%         if isempty(TF)
%             f_spikeFeature(KS_dirs{day})
%         end
%     end
    %% extracting spike data (this part will exrtact data for all the sessions)
    
%     for day=1:length(KS_dirs) %%  loop for different days
%         f_extract_spike_data_higher_vis(KS_dirs{day},BL_dirs{day})
% 
%     end
%     close all
    
    tic
    for day=1:length(KS_dirs) %%  loop for different days
        Plot=0;
        f_RF_single_unit_neonFull(KS_dirs{day},BL_dirs{day},Plot)
    end
close all
 
    %% plotting Master raster and psth (this will be applied on all the sessions)
    %%
    %
    % for day=1:length(KS_dirs) %% loop for different days
    %     f_master_raster_PSTH(KS_dirs{day})
    % end


    %% plotting direction tuning curve
%     
%     for day=1:length(KS_dirs) %% loop for different days
% %         f_randomDot_tuning_curve(KS_dirs{day})
%         f_neon6_tuning_curve(KS_dirs{day})
%     end
    
    
    %% plotting size tuning
    % for day=1:length(KS_dirs) %% loop for different days
    %     f_plot_size_tuning(KS_dirs{day})
    % end
    
    
    

end
toc
