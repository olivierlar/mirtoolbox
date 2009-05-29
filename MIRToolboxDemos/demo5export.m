function demo5export
% Example of use of the export function.
% For more examples, cf. mirfeatures.

s = mirspectrum('Folder');
if not(isempty(get(s,'Data')))
    ch = mirchromagram(s);
    k = mirkeystrength(ch);
    %h = mirhisto(k);
    e = mirenvelope('Folder');
    d = mirenvelope(e,'Diff');
    mirexport('NoFrameNoPeaks.txt',ch,k);

    ch = mirpeaks(ch,'Total',1);
    k = mirpeaks(k,'Total',1);
    %h = mirpeaks(h);
    e = mirpeaks(e,'Total',1);
    d = mirpeaks(d,'Total',1);
    mirexport('NoFrameButPeaks.txt',ch,k,e,d);
    clear ch k e d
    
    ac = mirautocor('Folder','Frame',.5,.2);
    s = mirspectrum('Folder','Frame',.5,.2);
    sf = mirflux(s);
    ms = mirspectrum(s,'Mel');
    m = mirmfcc(ms);
    ch = mirchromagram(s);
    k = mirkeystrength(ch);
    %h = mirhisto(k);
    e = mirenvelope('Folder','Frame',.5,.2);
    d = mirenvelope(e,'Diff');
    mirexport('FrameNoPeaks.txt',ms,m,ch,k);
    
    ac = mirpeaks(ac,'Total',1,'NoBegin');
    s = mirpeaks(s,'Total',1);
    sf = mirpeaks(sf,'Total',1);
    ms = mirpeaks(ms,'Total',1);
    m = mirpeaks(m,'Total',1,'Interpol',0);
    ch = mirpeaks(ch,'Total',1);
    k = mirpeaks(k,'Total',1);
    %h = mirpeaks(h)
    e = mirpeaks(e,'Total',1);
    d = mirpeaks(d,'Total',1);
    mirexport('FrameAndPeaks.txt',ac,s,sf,ms,m,ch,k,e,d);
end