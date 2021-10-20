function intersection_point=f_intersect_gauss2_for_hist_v2(vect,Nbin,aa,bb,report)
%%sugest an in terval for the intersection [aa,bb]
[y,xx] = hist(vect,Nbin);
y=y-1; y(y<0)=0;
%  y(xx>=1.3)=0;
idx= y==0;
xx=xx';
y=y';
y(idx)=nan;
 %% fitting
ff = fit(xx(~idx),y(~idx),'gauss2');
syms f(x,a,b,c)
f(x,a,b,c) = a*exp(-((x-b)/c)^2);
f1         = feval(f,xx,ff.a1,ff.b1,ff.c1);
f2         = feval(f,xx,ff.a2,ff.b2,ff.c2);
df1 = double(vpa(f1));
df2 = double(vpa(f2));

%% intersection
[~,ind_a] = min(abs(xx-aa));
[~,ind_b] = min(abs(xx-bb));

difff              =abs(df1-df2);
[~,index]          = min(difff(ind_a:ind_b));
index              = index+ind_a;
intersection_point = xx(index);

if report
    figure('units','normalized','outerposition',[0 0 .5 .8]);
    plot(ff,xx,y,'b*')
    hold on
    plot(xx,df1,'g','LineWidth',1)
    plot(xx,df2,'c','LineWidth',1)
    xline(xx(index),'-.',{'intersection point',num2str(intersection_point)},'LineWidth',2,'LabelVerticalAlignment','middle','LabelOrientation','horizontal');
end
% xlabel('trough-to-peak latency (ms)')
% ylabel('number of neurons')
% legend('Histogram data', 'Two-Term Gaussian fitted', 'First Gaussian','Second Gaussian','Intersection')
% saveas(gcf,[data_folder,'/neon_TPlatency_gauss2.png'])
% close gcf