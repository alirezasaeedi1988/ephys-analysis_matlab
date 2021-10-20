function RF_cell=f_RF_select(myKsDir,depth,report,exp)

% myKsDir='D:\matwork\data\higher_vision\M164_2020_02_12';
%input depth
% RSS_margin = 5;
fname=dir(fullfile(myKsDir,[exp,'*.mat']));

if length(fname)==2
   RF_2    = load(fullfile(myKsDir,fname(1).name));
   RF_1    = load(fullfile(myKsDir,fname(2).name));
   Ncell   = length(RF_1.RF_data.RF_size);
   RF_loc  = zeros(Ncell,2);
   RF_size = zeros(Ncell,1);
   RF_area = zeros(Ncell,1);   
   RF_intrsct_ellips = zeros(Ncell,1);
   RF_area_ellips    = zeros(Ncell,1);   
   RF_intrsct = zeros(Ncell,1);
   RF_HWHM_X  = zeros(Ncell,1);
   RF_HWHM_Y  = zeros(Ncell,1);
   RF_RSS  = zeros(Ncell,1);
   RF_label   = strings(Ncell,1);
   %% rejecting low RSS (goodness of fit) based on Zocolan paper
   RF_1.RF_data.RF_size(RF_1.RF_data.RSS<0.5)= nan;
   RF_2.RF_data.RF_size(RF_2.RF_data.RSS<0.5)= nan;
   %% selection based on nan rejecting
   bothnan = isnan(RF_1.RF_data.RF_size) & isnan(RF_2.RF_data.RF_size);
   nan1    = isnan(RF_1.RF_data.RF_size) & ~isnan(RF_2.RF_data.RF_size);
   nan2    = ~isnan(RF_1.RF_data.RF_size) & isnan(RF_2.RF_data.RF_size);
   %
   RF_loc(bothnan,:)  = nan;
   RF_size(bothnan)   = nan;
   RF_HWHM_X(bothnan) = nan;
   RF_HWHM_Y(bothnan) = nan;
   RF_label(bothnan)  = "NaN";
   RF_intrsct(bothnan)= nan;
   RF_area(bothnan)   = nan;
   RF_intrsct_ellips(bothnan)= nan;
   RF_area_ellips(bothnan)   = nan;
   RF_RSS(bothnan)    = nan;
   %
   RF_loc(nan1,:)  = RF_2.RF_data.RF_loc(nan1,:);
   RF_size(nan1)   = RF_2.RF_data.RF_size(nan1);
   RF_HWHM_X(nan1) = RF_2.RF_data.HWHM_X(nan1);
   RF_HWHM_Y(nan1) = RF_2.RF_data.HWHM_Y(nan1); 
   RF_intrsct(nan1)= RF_2.RF_data.intrsct_ratio(nan1);
   RF_area(nan1)   = RF_2.RF_data.area(nan1);
   RF_intrsct_ellips(nan1)= RF_2.RF_data.intrsct_ratio_ellips(nan1);
   RF_area_ellips(nan1)   = RF_2.RF_data.area_ellips(nan1);
   RF_RSS(nan1)    = RF_2.RF_data.RSS(nan1);
   RF_label(nan1)  = convertCharsToStrings(RF_2.RF_data.session_name);
   %
   RF_loc(nan2,:)  = RF_1.RF_data.RF_loc(nan2,:);
   RF_size(nan2)   = RF_1.RF_data.RF_size(nan2);
   RF_HWHM_X(nan2) = RF_1.RF_data.HWHM_X(nan2);
   RF_HWHM_Y(nan2) = RF_1.RF_data.HWHM_Y(nan2);
   RF_intrsct(nan2)= RF_1.RF_data.intrsct_ratio(nan2);
   RF_area(nan2)   = RF_1.RF_data.area(nan2);
   RF_intrsct_ellips(nan2)= RF_1.RF_data.intrsct_ratio_ellips(nan2);
   RF_area_ellips(nan2)   = RF_1.RF_data.area_ellips(nan2);   
   RF_RSS(nan2)    = RF_1.RF_data.RSS(nan2);
   RF_label(nan2)  = convertCharsToStrings(RF_1.RF_data.session_name);
   %% selection based on goodness of fit
   %    bestFit1 = RF_1.RF_data.RSS < RF_2.RF_data.RSS - RSS_margin & ~isnan(RF_1.RF_data.RF_size);
   %    bestFit2 = RF_2.RF_data.RSS < RF_1.RF_data.RSS - RSS_margin & ~isnan(RF_2.RF_data.RF_size);
   bestFit1 = RF_1.RF_data.RSS > RF_2.RF_data.RSS  & ~isnan(RF_1.RF_data.RF_size);
   bestFit2 = RF_2.RF_data.RSS > RF_1.RF_data.RSS  & ~isnan(RF_2.RF_data.RF_size);
   %
   RF_loc(bestFit1,:)  = RF_1.RF_data.RF_loc(bestFit1,:);
   RF_size(bestFit1)   = RF_1.RF_data.RF_size(bestFit1);
   RF_HWHM_X(bestFit1) = RF_1.RF_data.HWHM_X(bestFit1);
   RF_HWHM_Y(bestFit1) = RF_1.RF_data.HWHM_Y(bestFit1);   
   RF_intrsct(bestFit1)= RF_1.RF_data.intrsct_ratio(bestFit1);
   RF_area(bestFit1)   = RF_1.RF_data.area(bestFit1);   
   RF_intrsct_ellips(bestFit1)= RF_1.RF_data.intrsct_ratio_ellips(bestFit1);
   RF_area_ellips(bestFit1)   = RF_1.RF_data.area_ellips(bestFit1);   
   RF_RSS(bestFit1)   = RF_1.RF_data.RSS(bestFit1);   
   RF_label(bestFit1)  = convertCharsToStrings(RF_1.RF_data.session_name);
   %
   RF_loc(bestFit2,:)  = RF_2.RF_data.RF_loc(bestFit2,:);
   RF_size(bestFit2)   = RF_2.RF_data.RF_size(bestFit2);
   RF_HWHM_X(bestFit2) = RF_2.RF_data.HWHM_X(bestFit2);
   RF_HWHM_Y(bestFit2) = RF_2.RF_data.HWHM_Y(bestFit2); 
   RF_intrsct(bestFit2)= RF_2.RF_data.intrsct_ratio(bestFit2);
   RF_area(bestFit2)   = RF_2.RF_data.area(bestFit2);  
   RF_intrsct_ellips(bestFit2)= RF_2.RF_data.intrsct_ratio_ellips(bestFit2);
   RF_area_ellips(bestFit2)   = RF_2.RF_data.area_ellips(bestFit2); 
   RF_RSS(bestFit2)   = RF_2.RF_data.RSS(bestFit2);      
   RF_label(bestFit2)  = convertCharsToStrings(RF_2.RF_data.session_name);
   %% selection based on horizontal location, more centered one
   RF1_centered  = abs(RF_1.RF_data.RF_loc_shift(:,1)) < abs(RF_2.RF_data.RF_loc_shift(:,1)) & RF_size ==0;
   RF2_centered  = abs(RF_2.RF_data.RF_loc_shift(:,1)) < abs(RF_1.RF_data.RF_loc_shift(:,1)) & RF_size ==0;
   %
   RF_loc(RF1_centered,:)  = RF_1.RF_data.RF_loc(RF1_centered,:);
   RF_size(RF1_centered)   = RF_1.RF_data.RF_size(RF1_centered);
   RF_HWHM_X(RF1_centered) = RF_1.RF_data.HWHM_X(RF1_centered);
   RF_HWHM_Y(RF1_centered) = RF_1.RF_data.HWHM_Y(RF1_centered);
   RF_intrsct(RF1_centered)= RF_1.RF_data.intrsct_ratio(RF1_centered);
   RF_area(RF1_centered)   = RF_1.RF_data.area(RF1_centered);
   RF_intrsct_ellips(RF1_centered)= RF_1.RF_data.intrsct_ratio_ellips(RF1_centered);
   RF_area_ellips(RF1_centered)   = RF_1.RF_data.area_ellips(RF1_centered);
   RF_RSS(RF1_centered)   = RF_1.RF_data.RSS(RF1_centered);
   RF_label(RF1_centered)  = convertCharsToStrings(RF_1.RF_data.session_name);
   %
   RF_loc(RF2_centered,:)  = RF_2.RF_data.RF_loc(RF2_centered,:);
   RF_size(RF2_centered)   = RF_2.RF_data.RF_size(RF2_centered);
   RF_HWHM_X(RF2_centered) = RF_2.RF_data.HWHM_X(RF2_centered);
   RF_HWHM_Y(RF2_centered) = RF_2.RF_data.HWHM_Y(RF2_centered);
   RF_intrsct(RF2_centered)= RF_2.RF_data.intrsct_ratio(RF2_centered);
   RF_area(RF2_centered)   = RF_2.RF_data.area(RF2_centered);
   RF_intrsct_ellips(RF2_centered)= RF_2.RF_data.intrsct_ratio_ellips(RF2_centered);
   RF_area_ellips(RF2_centered)   = RF_2.RF_data.area_ellips(RF2_centered);
   RF_RSS(RF2_centered)   = RF_2.RF_data.RSS(RF2_centered);
   RF_label(RF2_centered)  = convertCharsToStrings(RF_2.RF_data.session_name);
   
elseif length(fname)==1
    RF_1    = load(fullfile(myKsDir,fname(1).name));
    Ncell   = length(RF_1.RF_data.RF_size);
    RF_label   = strings(Ncell,1);
    %% rejecting high RSSs
    RF_1.RF_data.RF_size(RF_1.RF_data.RSS<0.5)= nan;
    RF_loc  = RF_1.RF_data.RF_loc;
    RF_size = RF_1.RF_data.RF_size;
    RF_HWHM_X  = RF_1.RF_data.HWHM_X;
    RF_HWHM_Y  = RF_1.RF_data.HWHM_Y;
    RF_intrsct = RF_1.RF_data.intrsct_ratio;
    RF_area    = RF_1.RF_data.area;
    RF_intrsct_ellips = RF_1.RF_data.intrsct_ratio_ellips;
    RF_area_ellips    = RF_1.RF_data.area_ellips;
    RF_RSS    = RF_1.RF_data.RSS;
    RF_label(:)= convertCharsToStrings(RF_1.RF_data.session_name);
    
else
    error('check RF_data file!')
end

%% report
if report
  [~,idx] = sort(depth);
  sorted_depth=depth(idx);
  sorted_RF_loc = RF_loc(idx,:);
  nan_idx=isnan(sorted_RF_loc(:,1));
  sorted_RF_loc(nan_idx,:)=[];
  sorted_depth(nan_idx)=[];
  % B = smoothdata(sorted_RF_loc(:,1),'gaussian');  % this was not suficient
  Outlier_idx = isoutlier(sorted_RF_loc(:,1),'movmedian',4);
  sorted_RF_loc(Outlier_idx,:)=[];
  sorted_depth(Outlier_idx)=[];
  BB = smooth(sorted_RF_loc(:,1),'rlowess');
  figure;plot(sorted_RF_loc(:,1),sorted_depth,'-b*')
  hold on;plot(BB,sorted_depth,'-rd')
  ylabel('neouron depth (\mum)')
  xlabel('horizontal location on screen')
  title('nasal   <<<<------------------------->>>>   temporal')
  legend('data','smoothed data')
  
  saveas(gcf,fullfile(myKsDir,'RF_Horizontal location.png'))
end

%% area selection



if strcmp(exp,'FRF') %full screen RF for V1 recording
    Visual_areas  = ones(Ncell,1);
else
    Visual_areas  = nan(Ncell,1); %% area: 1 = V1, 2 = LM , 3 = LI , 4 = LL
    load(fullfile(myKsDir,'channel_RF_data.mat'))
    reversals=RF_chan.reversals ;

switch length(reversals)
    case 1
        Visual_areas(depth<RF_chan.chanDepth(reversals(1)))= 1;
        Visual_areas(depth>=RF_chan.chanDepth(reversals(1))) = 2;
    case 2
        Visual_areas(depth<RF_chan.chanDepth(reversals(1)))= 1;        
        Visual_areas(depth>=RF_chan.chanDepth(reversals(1))& depth<=RF_chan.chanDepth(reversals(2)))=2;
        Visual_areas(depth>RF_chan.chanDepth(reversals(2))) = 3;
    case 3
        Visual_areas(depth<RF_chan.chanDepth(reversals(1)))= 1;        
        Visual_areas(depth>=RF_chan.chanDepth(reversals(1))& depth<=RF_chan.chanDepth(reversals(2)))=2;   
        Visual_areas(depth>RF_chan.chanDepth(reversals(2))& depth<=RF_chan.chanDepth(reversals(3)))=3;
        Visual_areas(depth>=RF_chan.chanDepth(reversals(3)))=4;
end
end
%% select the boundaries
for l=1:Ncell
    if strcmp(RF_label(l),"NaN")
        X_Bound{l}= nan;
        Y_Bound{l}= nan;
    else
        switch length(strfind(RF_label(l),'2'))
            case 0
                X_Bound{l} = RF_1.RF_data.X_Bound{l};
                Y_Bound{l} = RF_1.RF_data.Y_Bound{l};
            case 1
                X_Bound{l} = RF_2.RF_data.X_Bound{l};
                Y_Bound{l} = RF_2.RF_data.Y_Bound{l};
        end
    end
end

%% saving 

RF_cell.RF_size   = RF_size;
RF_cell.RF_loc    = RF_loc;
RF_cell.RF_label  = RF_label;
RF_cell.area      = Visual_areas;   % visual area
RF_cell.RF_RSS    = RF_RSS;
RF_cell.RF_HWHM_X = RF_HWHM_X;
RF_cell.RF_HWHM_Y = RF_HWHM_Y;
RF_cell.X_Bound   = X_Bound;
RF_cell.Y_Bound   = Y_Bound;
RF_cell.monitor_size  = RF_1.RF_data.monitor_size;
RF_cell.intrsct_ratio = RF_intrsct;
RF_cell.RF_area   = RF_area; % size area
RF_cell.RF_area_ellips   = RF_area_ellips; % size area
RF_cell.intrsct_ratio_ellips = RF_intrsct_ellips;

        save([myKsDir,'\selected_RF_data.mat'],'RF_cell')

