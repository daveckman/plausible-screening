function draw_hist(X,varargin)
if nargin > 1
    bins = varargin{1};
else
    bins = linspace(min(X),max(X),sqrt(length(X)));
end
hist(X,bins)
h = findobj(gca,'Type','patch');
set(h,'facecolor',[0.3 0.8 0.8],'edgecolor','none')
h = findobj(gca,'Type','patch');
hold on
[f,xi] = ksdensity(X,min(X):((max(X)-min(X))/1000):max(X));
hold on
plot(xi,f*(bins(2)-bins(1))*length(X),'linewidth',2)
xlim([min(bins) max(bins)])
set(gca,'ytick',[])