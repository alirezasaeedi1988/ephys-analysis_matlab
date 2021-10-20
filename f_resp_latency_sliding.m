function teta=f_resp_latency_sliding(y,onset_ind,offset_ind,pre_stimulus_act,win_len,step_size)

% win_len=20;
% step_size=4
start_time=onset_ind-length(pre_stimulus_act);
end_time=start_time+win_len;
hv=zeros(1,ceil((offset_ind-(start_time+win_len))/step_size));
i=0;
 while end_time<offset_ind
    i=i+1;
    [~,hv(i)]=ranksum(pre_stimulus_act,y(start_time:end_time));
    start_time=start_time+step_size;
    end_time=end_time+step_size;
   
     
 end


N = 3; % Required number of consecutive numbers following a first one
t=find(hv); 
x = diff(t)==1;
f = find([false,x]~=[x,false]);
g = find(f(2:2:end)-f(1:2:end-1)>=N-1,1,'first');
first_t = t(f(2*g-1)); %% First t followed by >=N consecutive numbers
teta=first_t*step_size+win_len-length(pre_stimulus_act);