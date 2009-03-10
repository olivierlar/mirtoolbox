function r = miropt_pulseclarity(x)

if not(isa(x,'miraudio'))
    a = miraudio('Design'); 
else
    a = x;
end
r = mirstruct;
%h = waitbar(0,'Preparing design...');
tot = 0;
jj = 0;
kk = 0;
fea = {'Envelope','SpectralFlux','Pitch'};
for n = 1:length(fea)
    if strcmp(fea{n},'Envelope')
        nc = 20;
    else
        nc = 1;
    end
    r.tmp.onsets{n} = mironsets(a,'Diff',0,'Filterbank',nc,...
                                            fea{n},'Sum',0,'Detect',0);
    if strcmp(fea{n},'Envelope')
        r.tmp.diffonsets{n} = mironsets(r.tmp.onsets{n},'Diff',1,...
            'Sum',0,'Detect',0); % 'SUM' AND 'DETECT' OPTIONS SHOULD BE AVOIDED HERE
    else
        r.tmp.diffonsets{n} = r.tmp.onsets{n};
    end
    if nc > 1
        su = {'Before','After'};
    else
        su = {'Before'};
    end
    for m = 1:length(su)
        jj = jj+1;
        [r.fea{n}.res{1}.su{m}.max r.tmp.ac{jj}] = ...
            mirpulseclarity(r.tmp.diffonsets{n},...
                 fea{n},'Enhanced',0,'Resonance',0,'Sum',su{m},...
                 'MaxAutocor');
             
        res = {0,'ToiviainenSnyder','vonNoorden'};
        for k = 1:length(res)
            kk = kk+1;
            r.tmp.ac2{kk} = mirautocor(r.tmp.ac{jj},'Resonance',res{k});
            if k > 1
                r.fea{n}.res{k}.su{m}.max = ...
                    mirpulseclarity(r.tmp.ac2{kk},'MaxAutocor');
            end
            r.fea{n}.res{k}.su{m}.min = ...
                mirpulseclarity(r.tmp.ac2{kk},'MinAutocor');
            r.fea{n}.res{k}.su{m}.kurtosis = ...
                mirpulseclarity(r.tmp.ac2{kk},'KurtosisAutocor','Total',1);
            r.fea{n}.res{k}.su{m}.entropy0 = ...
                mirpulseclarity(r.tmp.ac2{kk},'EntropyAutocor');
            r.fea{n}.res{k}.su{m}.interf0 = ...
                mirpulseclarity(r.tmp.ac2{kk},'InterfAutocor','Total',Inf);
            r.tmp.ac3{kk} = mirautocor(r.tmp.ac2{kk},'Enhanced');
            r.fea{n}.res{k}.su{m}.entropy1 = ...
                mirpulseclarity(r.tmp.ac3{kk},'EntropyAutocor');
            r.fea{n}.res{k}.su{m}.interf = ...
                mirpulseclarity(r.tmp.ac3{kk},'InterfAutocor','Total',Inf);
            r.fea{n}.res{k}.su{m}.tempo = ...
                mirpulseclarity(r.tmp.ac3{kk},'TempoAutocor');
            tot = tot+6;
        end
    end
    %waitbar((n-1)/length(fea)+...
    %        (k-1)/length(res)/length(fea),...
    %        ['Preparing design (',num2str(tot),' operations designed)...']);
    if strcmpi(fea{n},'Envelope')
        r.tmp.onsets2 = mirsum(mironsets(r.tmp.onsets{n},'Detect'));
        r.xtrem = mirpulseclarity(r.tmp.onsets2,'ExtremEnvelop');
        r.attack = mirpulseclarity(r.tmp.onsets2,'Attack','Diff');
        r.attack2 = mirpulseclarity(r.tmp.onsets2,'Attack','Gauss');
    end
end
r.artic = mirpulseclarity(a,'Articulation');
%close(h)
drawnow

if not(isa(x,'miraudio'))
    r = mireval(r,x,'Single');
end

mirexport('pulseclarity_predictions.txt',r);