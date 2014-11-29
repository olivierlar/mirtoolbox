function mov = display(k,filename)
% KEYSOM/DISPLAY display of key som

mov = [];
load keysomaudiodata; 
w = get(k,'Weight');
n = get(k,'Name');
if nargin>1
    mov = avifile(filename);
%elseif and(nargout>0,size(w)>0)
%    mov = moviein(size(w{1}{1},2)); %no longer needed as of MATLAB Release 11 (5.3).
end
for i = 1:length(w)
    wi = w{i};
    fig = figure;
    for j = 1:length(wi)
        wj = wi{j};
        for k = 1:size(wj,2)
            for l = 1:size(wj,3)
                h = pcolor(squeeze(wj(:,k,l,:)));
                shading interp
                axis([1,36,1,24]), view(2) , caxis([-1 1])
                axis off;
                hold on
                for m=1:24 
                    text(0.99*keyx(m)-1, 0.98*keyy(m)+1, keyN(m,:),...
                                            'FontSize',16,'FontName','Arial'); 
                end
                hold off
                set(gca,'PlotBoxAspectRatio',[1.5 1 1])
                colormap('jet')
                title('Self-organizing map projection of chromagram') 
                drawnow
                if nargin>1
                    mov = addframe(mov,gca);
                elseif nargout>0
                    colormap('jet')
                    if k == 1
                        mov = getframe;
                    else
                        mov(k) = getframe; 
                    end
                end
            end
        end
    end
    if isa(fig,'matlab.ui.Figure')
        fig = fig.Number;
    end
    disp(['The key som related to file ',n{i},...
                ' is displayed in Figure ',num2str(fig),'.']);
end
if nargin>1
    mov = close(mov);
    disp(['Data exported to file ',filename,'.']);
end
disp(' ');