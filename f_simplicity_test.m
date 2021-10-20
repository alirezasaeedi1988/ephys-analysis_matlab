function [pv_simplicity,h_simplicity]=f_simplicity_test(phy_smoothed_pdf, xp_phy, spike, psth_locs)


%% inputs: phy_smoothed_pdf; xp_phy, spike, psth_locs
% mean_val=mean(phy_smoothed_pdf);
% phy_smoothed_pdf=(phy_smoothed_pdf-mean_val)/mean_val;
len=length(phy_smoothed_pdf);
% [~, m_loc]= max( phy_smoothed_pdf(1+floor(len/10):floor(len/2)));
% m_loc=m_loc+1+floor(len/10);
% window1=[xp_phy(m_loc-floor(len/10)),xp_phy(m_loc+floor(len/10))];
% window2=[xp_phy(m_loc+floor(len/5)),xp_phy(m_loc+floor(2*len/5))];
[~, m_loc]= max( phy_smoothed_pdf(floor(len/2):end-(1+floor(len/10))));
m_loc=m_loc+floor(len/2);
window1=[xp_phy(m_loc-floor(len/10)),xp_phy(m_loc+floor(len/10))];
window2=[xp_phy(m_loc-floor(2*len/5)),xp_phy(m_loc-floor(len/5))];
Ntrials=max(spike(:,2)-1);
test_vect=zeros(2,Ntrials);
% spik_time_loc=logical(double(spike(:,3)== psth_locs(2)) +double(spike(:,3)== psth_locs(3)));
spik_time_loc=spike(:,3)== psth_locs(2);
spikee=spike(spik_time_loc,1:2);
for i =1:Ntrials %%this loop is over trials
    spike_in_trial=spikee(spikee(:,2)==i,1);
    
    test_vect(1,i)=length(find(spike_in_trial<=window1(1) & spike_in_trial<=window1(2)));
    test_vect(2,i)=length(find(spike_in_trial<=window2(1) & spike_in_trial<=window2(2)));
    
end

[pv_simplicity,h_simplicity]= ranksum(test_vect(1,:),test_vect(2,:),'alpha',0.05);