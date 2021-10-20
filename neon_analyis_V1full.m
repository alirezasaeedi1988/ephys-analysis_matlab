%%%% neon color fullscreen in V1
%%% extracting stimulus data and spike data and ploting raster, psth and
%%% tuning curve...
clear
tic
%% parameters
%%
data_folder= 'D:\matwork\Data\NEON_V1';
data_folder2='Y:\MPI_EphysData_2020\v1';
mouse_numbers =[ 185 335];

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
if mouse_number==185
 KS_dirs(6)=[];
 KS_dirs(3)=[];
end
  %%% BL folders

BL_dirs=cell(2,1);
% folders=dir(fullfile(data_folder2,mouse));
folders=dir(data_folder2);
j=1;
for i=1:length(folders)
    if length(folders(i).name)>4 && strcmp(folders(i).name(1:4),mouse)
        BL_dirs{j}=fullfile(data_folder2,folders(i).name);
%         OE_dirs{j}=fullfile(data_folder2,folders(i).name);
        j=j+1;
    end 
end

%% extracting stimulus data (this part will exrtact data for all the sessions)
%%
% 
% for day=1:length(KS_dirs) %%  loop for different days
%    f_stim_data_blckrck(OE_dirs{day},KS_dirs{day+7},Nchan);
% end

%% extracting spike data (this part will exrtact data for all the sessions)

for day=1:length(KS_dirs) %%  loop for different days
    f_extract_spike_data_V1full(KS_dirs{day},BL_dirs{day})
end


%% plotting tuning curves
for day=1:length(KS_dirs) %% loop for different days
%     f_plot_size_tuning(KS_dirs{day})
    f_neon6_tuning_curve_V1full(KS_dirs{day})
end
 

end
toc
