function varargout = mircompute(algo,varargin)

l = length(varargin{1});    % number of audio files.
for i = 1:l
    % for each audio file
    v = varargin;
    for j = 1:length(v)
        % for each variable associated to that audio
        if isa(v{j},'mirdata')
            v{j} = get(v{j},'Data');
        end
        if iscell(v{j})
            v{j} = v{j}{i};
        end
        %if not(iscell(v{j}))
        %    v{j} = {v{j}};
        %end
    end
    % final result for that audio
    for k = 1:length(v{1})
        % for each segment in that audio file
        vk = v;
        for j = 1:length(vk)
            if iscell(vk{j})
                vk{j} = vk{j}{k};
            end
        end
        if nargout == 1
            res = algo(vk{:});
        else 
            [res res2] = algo(vk{:}); 
        end
        if iscell(res)
            lr = length(res);
            for j = 1:lr
                varargout{j}{i}{k} = res{j};
            end
        else
            varargout{1}{i}{k} = res;
            if nargout == 2
                varargout{2}{i}{k} = res2;
            end
        end
    end
end