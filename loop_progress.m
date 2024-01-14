function loop_progress(i,imax,varargin)

%% VERY SIMPLE MONITORING OF LOOP PROGRESS
%  assuming that each pass in a for loop requires 
%  a constant period of time

%% INPUT:
%  i = index task progression
%  imax = index of task completion
    
%% INSTRUCTION FOR USE
%
%  for i=1:imax
%      loop_progress(i,imax)
%      ... 
%      core job
%      ...
%  end
%  loop_progress(i,imax,'done') 

%  Note that loop_progress(i>imax,imax)
%  has the samme effect as
%  loop_progress(i,imax,'done') 

persistent t0 period

for k = 1:length(varargin)
    switch lower(varargin{k}) % varargin is a cell array
      case {'done'}
        i=imax+1;
      otherwise
        error('option %s not recognized',varargin{k});
    end
end

if (i>imax)
    for j=1:43
        fprintf('\b')
    end
    fprintf('100%%   | Time to completion: %s\n\n', ...
            time2str(0));
end

per = 100*i/imax;

t = clock();
if isempty(t0)
    t0 = t;
else
    if (i==2)
        period=etime(clock(),t0);
        fprintf('Loop period = %d [s]\n',period)
        fprintf('Expected loop duration = %d [s]\n',period*imax)
    end
    if (i>=2 && i<imax)
        for j=1:43
            fprintf('\b')
        end
        fprintf('%05.2f%%',per)  
        dur = period *(imax-i+1);
        fprintf(' | Time to completion: %s',time2str(dur)); 
    end
end

end % End of function loop_progress