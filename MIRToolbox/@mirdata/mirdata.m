function d = mirdata(orig,varargin)
%   d = mirdata(x) creates a MIR object.

if nargin > 0 && isa(orig,'mirdata')
    d.pos = orig.pos;
    d.data = orig.data;
    d.unit = orig.unit;
    d.framepos = orig.framepos;
    d.framerate = orig.framerate;
    d.framed = orig.framed;
    d.sr = orig.sr;
    d.length = orig.length;
    d.nbits = orig.nbits;
    d.name = orig.name;
    d.name2 = orig.name2;
    d.label = orig.label;
    d.channels = orig.channels;
    d.clusters = orig.clusters;
    d.multidata = orig.multidata;
    d.peak = orig.peak;
    d.onset = orig.onset;
    d.release = orig.release;
    d.track = orig.track;
    d.title = orig.title;
    d.abs = orig.abs;
    d.ord = orig.ord;
    d.interchunk = orig.interchunk;
    d.tmpidx = orig.tmpidx;
    d.acrosschunks = orig.acrosschunks;
    d.interpolable = orig.interpolable;
    d.tmpfile = orig.tmpfile;
    d.index = orig.index;
else
    d.pos = {};
    d.data = {};
    d.unit = '';
    d.framepos = {};
    d.framerate = {};
    d.framed = 0;
    d.sr = {};
    d.length = {};
    d.nbits = {};
    d.name = {};
    d.name2 = {};
    d.label = {};
    d.channels = [];
    d.clusters = {};
    d.multidata = [];
    d.peak.pos = {};
    d.peak.val = {};
    d.peak.precisepos = {};
    d.peak.preciseval = {};
    d.peak.mode = {};
    d.onset = {};
    d.release = {};
    d.track = {};
    d.title = 'Unspecified data';
    d.abs = 'Unspecified abscissa';
    d.ord = 'Unspecified ordinate';
    d.interchunk = [];
    d.tmpidx = 0;
    d.acrosschunks = [];
    d.interpolable = 1;  % If the abscissae axis is non-numeric (0), 
                         % then peak picking has to be done without interpolation.
    d.tmpfile = [];
    d.index = NaN;
end
d = class(d,'mirdata');
d = set(d,varargin{:});