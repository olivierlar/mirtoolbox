function r = miropt_pulseclarity(x)

vb = mirverbose;
mirverbose(0);
if not(isa(x,'miraudio'))
    a = miraudio('Design'); 
else
    a = x;
end
r = mirstruct;
tot = 0;
jj = 0;
kk = 0;
ll = 0;
mm = 0;
fea = {'Envelope','Pitch'};
for n = 1:length(fea)
    disp(fea{n});
    if strcmpi(fea{n},'Envelope')
        nc = {'Gammatone','Scheirer','Klapuri','Spectro'};
    %elseif strcmpi(fea{n},'Flux')
    %    nc = {'Flux'};
    else
        nc = {'Pitch'};
    end
    for l = 1:length(nc)
        disp(nc{l});
        if strcmpi(nc{l},'Pitch') %|| strcmpi(nc{l},'Flux')
            ft = {nc{l}};
        elseif strcmpi(nc{l},'Spectro')
            ft = {'Freq','Mel'};
        else
            ft = {'IIR','HalfHann'};
        end
        for j = 1:length(ft)
            disp(ft{j});
            ll = ll+1;
            if strcmpi(nc{l},'Pitch')
                r.tmp.onsets{ll} = mironsets(a,'Diff',0,'Pitch',...
                            'Sum',0,'Detect',0);
            %elseif strcmpi(nc{l},'Flux')
            %    r.tmp.onsets{ll} = mirspectrum(a,'Frame',.05,.5);
            elseif strcmpi(nc{l},'Spectro')
                r.tmp.onsets{ll} = mironsets(a,'Diff',0,...
                            'Method','Spectro',ft{j},...'UpSample',...
                            'Sum',0,'Detect',0,'Filterbank',0);
            elseif strcmpi(ft{j},'IIR')
                r.tmp.onsets{ll} = mironsets(a,'Diff',0,...
                            'FilterbankType',nc{l},'Filterbank',20,...
                            'Envelope','Sum',0,'Detect',0);
            else
                r.tmp.onsets{ll} = mironsets(a,'Diff',0,...
                            'FilterbankType',nc{l},'Filterbank',20,...
                            'FilterType','HalfHann','PreDecim',180,...
                            'Envelope','Sum',0,'Detect',0);
            end
            if strcmpi(nc{l},'Pitch')
                su = {'After'};
                llg = 1;    % No log at all
                order = 1;
                hwr = 0;
            %elseif strcmpi(nc{l},'Flux')
            %    su = {'After'};
            %    llg = 1;    % No log at all
            %    order = 1;
            %    hwr = [0 1];
            else
                if strcmpi(nc{l},'Spectro')
                    r.tmp.logonsets{ll} = ...
                        mironsets(r.tmp.onsets{ll},'Power','Detect',0,'Sum',0);
                else
                    r.tmp.logonsets{ll} = ...
                        mironsets(r.tmp.onsets{ll},'Log','Detect',0,'Sum',0);
                end
                llg = 1:2;
                hwr = [0 .8 1];
                if strcmpi(ft{j},'Freq')
                    su = {'Before'};
                elseif strcmpi(nc{l},'Scheirer') || strcmpi(nc{l},'Spectro')
                    su = {'Before','After'};
                else
                    su = {'Before','Adjacent','After'};
                end
            end
            for hw = 1:length(hwr)
                for lg = 1:length(llg)
                    if strcmpi(ft{j},'Freq') && (hwr(hw)<1)
                        order = [1 8 i];
                    else
                        order = [1 8];  % Differentiator order
                    end
                    for o = 1:length(order)
                        mm = mm+1;
                        if llg(lg) == 2
                            logonsets = r.tmp.logonsets{ll};
                        else
                            logonsets = r.tmp.onsets{ll};
                        end
                        if 0 %strcmpi(nc{l},'Flux')
                            if hwr(hw)
                                r.tmp.modifonsets{mm} = ...
                                    mirflux(logonsets,'Complex');
                            else
                                r.tmp.modifonsets{mm} = ...
                                    mirflux(logonsets);
                            end
                        else
                            if hwr(hw)
                                r.tmp.modifonsets{mm} = ...
                                    mironsets(logonsets,...
                                        'Complex',order(o) == i,...
                                        'HalfwaveDiff',abs(order(o)),...
                                        'Lambda',hwr(hw),...
                                        'Sum',0,'Detect',0);
                            else
                                r.tmp.modifonsets{mm} = ...
                                    mironsets(logonsets,...
                                        'Complex',order(o) == i,...
                                        'Diff',abs(order(o)),...
                                        'Sum',0,'Detect',0);
                            end
                        end
                        %mirsum(r.tmp.modifonsets{mm})
                        for m = 1:length(su)
                            jj = jj+1;
                            [r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}.hwr{hw}...
                                .res{1}.su{m}.max r.tmp.ac{jj}] = ...
                                mirpulseclarity(r.tmp.modifonsets{mm},...
                                     'Enhanced',0,'Resonance',0,...
                                     'Sum',su{m},'MaxAutocor');
                            res = {0,'ToiviainenSnyder','vanNoorden'};
                            for k = 1:length(res)
                                kk = kk+1;
                                
                                r.tmp.ac2{kk} = mirautocor(r.tmp.ac{jj},...
                                                    'Resonance',res{k});
                                %r.tmp.ac2{kk},drawnow
                                
                                if k > 1
                                    r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                        .hwr{hw}.res{k}.su{m}.max = ...
                                        mirpulseclarity(r.tmp.ac2{kk},...
                                            'MaxAutocor');
                                end
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %        .hwr{hw}.res{k}.su{m}.max
                                    
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.min = ...
                                    mirpulseclarity(r.tmp.ac2{kk},...
                                            'MinAutocor');
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.min
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.kurtosis = ...
                                    mirpulseclarity(r.tmp.ac2{kk},...
                                            'KurtosisAutocor','Total',1);
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.kurtosis
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.entropy0 = ...
                                    mirpulseclarity(r.tmp.ac2{kk},...
                                            'EntropyAutocor');
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.entropy0
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.interf0 = ...
                                    mirpulseclarity(r.tmp.ac2{kk},...
                                            'InterfAutocor','Total',Inf);
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.interf0
                                
                                r.tmp.ac3{kk} = mirautocor(r.tmp.ac2{kk},...
                                                    'Enhanced');
                                %r.tmp.ac3{kk},drawnow
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.entropy1 = ...
                                    mirpulseclarity(r.tmp.ac3{kk},...
                                            'EntropyAutocor');
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.entropy1
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.interf = ...
                                    mirpulseclarity(r.tmp.ac3{kk},...
                                            'InterfAutocor','Total',Inf);
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.interf
                                
                                r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                    .hwr{hw}.res{k}.su{m}.tempo = ...
                                    mirpulseclarity(r.tmp.ac3{kk},...
                                            'TempoAutocor');
                                %r.flog{n}.nc{l}.ft{j}.ord{o}.log{lg}...
                                %    .hwr{hw}.res{k}.su{m}.tempo
                                tot = tot+6;
                            end
                        end
                        if strcmpi(fea{n},'Envelope') && hw == 1 && lg == 1
                            r.tmp.sumdiffonsets{mm} = ...
                                mironsets(mirsum(...
                                    mironsets(r.tmp.modifonsets{mm})),...
                                          'Detect','Contrast',.3);
                            %r.tmp.sumdiffonsets{mm},drawnow
                            
                            r.flog{n}.nc{l}.ft{j}.ord{o}.attack = ...
                                mirpulseclarity(r.tmp.sumdiffonsets{mm},...
                                                    'AttackDiff');
                            %r.flog{n}.nc{l}.ft{j}.ord{o}.attack
                        else
                            r.tmp.sumdiffonsets{mm} = r.tmp.modifonsets{mm};
                        end
                    end
                end
            end
            
            if strcmpi(fea{n},'Envelope')
                r.tmp.sumonsets{ll} = mironsets(mirsum(r.tmp.onsets{ll}),...
                                                'Detect','Contrast',.3);
                %r.tmp.sumonsets{ll},drawnow
                
                r.fea{n}.nc{l}.ft{j}.xtrem = ...
                    mirpulseclarity(r.tmp.sumonsets{ll},'ExtremEnvelop');
                %r.fea{n}.nc{l}.ft{j}.xtrem
                
                r.fea{n}.nc{l}.ft{j}.attack = ...
                    mirpulseclarity(r.tmp.sumonsets{ll},'Attack','Diff');
                %r.fea{n}.nc{l}.ft{j}.attack
                
                r.fea{n}.nc{l}.ft{j}.attack2 = ...
                    mirpulseclarity(r.tmp.sumonsets{ll},'Attack','Gauss');
                %r.fea{n}.nc{l}.ft{j}.attack2
            else
                r.tmp.sumonsets{ll} = r.tmp.onsets{ll};
            end
        end
    end
end
r.artic = mirpulseclarity(a,'Articulation');

if not(isa(x,'miraudio'))
    r = mireval(r,x,'Single','pulseclarity_predict.txt');
end

mirverbose(vb);