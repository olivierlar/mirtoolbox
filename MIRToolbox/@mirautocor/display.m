function display(a,varargin)

if nargin>1 && isa(varargin{1},'mirdata')
    x = varargin{1};
    x = purgedata(x);
    display(x);
    p = get(a,'PeakPos');
    ph = get(a,'Phase');
    px = get(x,'Pos');
    dx = get(x,'Data');
    for i = 1:length(p)
        for j = 1:length(p{i})
            y = 0;
            maxx = max(max(dx{i}{j}));
            step = maxx/100;
            for k = 1:length(p{i}{j})
                for h = 1:length(p{i}{j}{k})
                    peak = p{i}{j}{k}(h);
                    phas = ph{i}{j}(peak,k);
                    pnts = phas:peak-1:size(px{i}{j},1);
                    y = y-1;
                    plot(px{i}{j}([1 end],k),y*step*[1 1],'r');
                    plot(px{i}{j}(pnts,k),y*step*ones(size(pnts)),'*r');
                end
            end
        end
    end
else
    display(mirdata(a),varargin{:});
end