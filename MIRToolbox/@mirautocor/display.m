function display(a,varargin)

display(mirdata(a),varargin{:});

ph = get(a,'Phase');
if ~isempty(ph)
    sr = get(a,'Sampling');
    x = get(a,'Input');
    x = purgedata(x);
    display(x);
    p = get(a,'PeakPosUnit');
    px = get(x,'Pos');
    dx = get(x,'Data');
    for i = 1:length(p)
        for j = 1:length(p{i})
            y = 0;
            b = [];
            maxx = max(max(dx{i}{j}));
            step = maxx/100;
            for k = 1:length(p{i}{j})
                for h = 1:length(p{i}{j}{k})
                    if ~isempty(b)
                        ok = 1;
                        kk = k;
                        hh = h;
                        for l = 1:length(b)
                            if px{i}{j}(1,kk) <= b(l)+.1
                                ok = 0;
                                break
                            end
                            if hh == length(p{i}{j}{kk})
                                kk = kk+1;
                                hh = 1;
                                if kk > length(p{i}{j})
                                    break
                                end
                            else
                                hh = hh+1;
                            end
                        end
                        if ok
                            y = 0;
                            b = [];
                        end
                    end
                    
                    lag = round(p{i}{j}{k}(h)*sr{i});
                    phas = ph{i}{j}{k}(h);
                    pnts = phas:lag:size(px{i}{j},1);
                    plot(px{i}{j}([1 end],k),-y*step*[1 1],'r');
                    plot(px{i}{j}(pnts,k),-y*step*ones(size(pnts)),'*r');
                    b = [b px{i}{j}(end,k)];
                    y = y+1;
                end
            end
        end
    end
end