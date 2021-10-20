function [preferred_angle,preferred_ind,osi,dsi,inter_pref_angle]=f_find_preferred(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,type)
%% find_preferred angle, OSI and DSI
%% inputs:Pdf,stim_angle,angle_loc,onset_ind,offset_ind
% type= 'all' , 'neon' , 'phy'

ave_activity_angle    = zeros(1,length(stim_angle));
interpolated_angle = min(stim_angle):15:max(stim_angle);

pre=mean(mean(Pdf(:,1:onset_ind-1)));

for ang=1:length(stim_angle) %%% loop over angles to find the preferred angle
    switch type
        case 'all'
            p=mean(Pdf(angle_loc(ang,1:2),:)); %-1 is to reject ctrl
        case 'neon'
            p=(Pdf(angle_loc(ang,1),:)); 
        case 'phy'
            p=(Pdf(angle_loc(ang,2),:)); %-1 is to reject ctrl
    end

    %             pre=mean(p(1:onset_ind));
    p=abs((p-pre)/pre);
    ave_activity_angle(ang)=mean(p(onset_ind:offset_ind));
end
[~,preferred_ind]=max(ave_activity_angle);
preferred_angle = stim_angle(preferred_ind);

%% OSI and DSI
ave_activity_angle=abs(ave_activity_angle);
opp_ind=mod(preferred_ind+4,8);
if opp_ind ==0
    opp_ind=8;
end
dsi=(ave_activity_angle(preferred_ind)-ave_activity_angle(opp_ind))/(ave_activity_angle(preferred_ind)+ave_activity_angle(opp_ind));
ort_ind=mod(preferred_ind+2,8);
if ort_ind ==0
    ort_ind=8;
end
osi=(ave_activity_angle(preferred_ind)-ave_activity_angle(ort_ind))/(ave_activity_angle(preferred_ind)+ave_activity_angle(ort_ind));

%% interpolation for activities
interp_act = interp1(stim_angle,ave_activity_angle,interpolated_angle,'spline');
[~,inter_pref_ind]=max(interp_act);
inter_pref_angle = interpolated_angle(inter_pref_ind);
