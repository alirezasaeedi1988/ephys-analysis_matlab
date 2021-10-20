function teta=f_resp_latency_5p(y,onset_ind)



y=abs(y(onset_ind:onset_ind*2));
[Max,loc]=max(y);
search_y=y(1:loc)./Max;
f=find(search_y<.05);
if ~isempty(f)
    teta=f(end);
else
    teta=nan;
end