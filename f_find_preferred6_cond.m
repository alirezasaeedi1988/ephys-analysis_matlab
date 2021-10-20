function [preferred_angle,preferred_ind,osi,dsi]=...
    f_find_preferred6_cond(spike_data,stim_angle,angle_loc,onset_ind,offset_ind,cond,myKsDir,neuron,Report)
%% find_preferred,OSI and DSI
%% inputs  : spike_data, stim_angle, angle_loc,onset_ind,offset_ind,cond,myKsDir,neuron,cid,Report
%% outputs : preferred_angle,preferred_ind,osi,dsi
% cond= 0 for overall condition



id=spike_data.evoked_cids{neuron};
Pdf=spike_data.pdfs{neuron};

ave_activity_angle    = zeros(1,length(stim_angle));

pre=mean(mean(Pdf(:,1:onset_ind-1)));

for ang=1:length(stim_angle) %%% loop over angles to find the preferred angle
    switch cond 
        case 0 %zero for all
            p = mean(Pdf(angle_loc(ang,:),onset_ind:offset_ind));
        otherwise
            p = Pdf(angle_loc(ang,cond),onset_ind:offset_ind);
    end
    p = p-pre;
    ave_activity_angle(ang)= mean(p);
end
inv_activity_angle = ave_activity_angle.*sign(mean(ave_activity_angle)); %% to invert the activity for inhibited cells

[~,preferred_ind] = max(inv_activity_angle);
preferred_angle   = stim_angle(preferred_ind);

x_fit   =1:length(stim_angle);
y_fit   = circshift(inv_activity_angle,5-preferred_ind); %% shift the peak to the 5th element 
g       = @(A,X) A(1)*exp( -((X-5).^2/(2*A(2)^2)))+ A(3)*exp( -((X-5+4).^2/(2*A(4)^2)))+A(5); %% double gaussian with pre-defined mean valeus
A0      = [.5,1,.5,1, mean(y_fit)];  %% initial guess
[A,RSS] = lsqcurvefit(g,A0,x_fit,y_fit); %% fitting
Rsq     = 1-(RSS/(rssq(inv_activity_angle)^2));
new_y   = A(1)*exp( -((x_fit-5).^2/(2*A(2)^2)))+ A(3)*exp( -((x_fit-5+4).^2/(2*A(4)^2)))+A(5); %% generat data with the fitted function
new_y   = circshift(new_y,-(5-preferred_ind)); %% shift the peak to original location

% ffit = fit(stim_angle',ave_activity_angle','gauss2');
if Report
    session_name='neon';
    outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end
    figure('units','normalized','outerposition',[0 0 0.5 .8],'Visible','off');
    plot(stim_angle,inv_activity_angle,'-b*')
    hold on
    plot(stim_angle,new_y,'r')
    plot(stim_angle,ave_activity_angle,'-k*')
    legend('inverted data','gauss2','raw data')
    xticks(stim_angle)
    title(['neuron-',num2str(id),' cond=',num2str(cond),' RSQ=',num2str(Rsq)])
    saveas(gcf,[outputDir,'\tunning_fit_',num2str(id),'_',num2str(cond),'.png'])
    close gcf
end
%% OSI and DSI

% new_y   = abs(new_y);  %% instead of this shift by minimum
new_y   = new_y - min(new_y);
opp_ind = mod(preferred_ind+4,8);
if opp_ind ==0
    opp_ind=8;
end
dsi     = (new_y(preferred_ind)-new_y(opp_ind))/(new_y(preferred_ind)+new_y(opp_ind));
ort_ind = mod(preferred_ind+2,8);
if ort_ind ==0
    ort_ind=8;
end
osi     = (new_y(preferred_ind)-new_y(ort_ind))/(new_y(preferred_ind)+new_y(ort_ind));


%% interpolation for activities
% interpolated_angle = min(stim_angle):5:max(stim_angle);
%
% interp_act = interp1(stim_angle,ave_activity_angle,interpolated_angle,'spline');
% [~,inter_pref_ind]=max(interp_act);
% inter_pref_angle = interpolated_angle(inter_pref_ind);
%
% g = @(A,X) A(1)*exp( -((X-inter_pref_angle).^2/(2*A(2)^2)))+ A(3)*exp( -((X-inter_pref_angle+180).^2/(2*A(4)^2)))+A(5);
% A0=[1,1,1,1, mean(ave_activity_angle)];
% [A,RSS] = lsqnonlin(g,A0,interpolated_angle,interp_act);
% Rsq=1-(RSS/(rssq(interp_act)^2));
% new_y=A(1)*exp( -((interpolated_angle-inter_pref_angle).^2/(2*A(2)^2)))+ A(3)*exp( -((interpolated_angle-inter_pref_angle+180).^2/(2*A(4)^2)))+A(5);
%
% plot(interpolated_angle,interp_act);xticks(stim_angle)
% hold on
% plot(interpolated_angle,new_y);xticks(stim_angle)