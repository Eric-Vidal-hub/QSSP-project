function savepdf(file_name,varargin)

%% PDF OUTPUTS OF CURRENT GRAPH OF A DATA FILE
%
%  WARNING: In order to work, a figure must have been opened 
%           by the calling script !

%% INPUT
%  file_name = name of data file currently plotted,
%              featuring a file type extension delimited by period, e.g.: 
%              example.dat, example.pdf, etc.
%              Files without any extension can also be processed.

%% OUTPUT
%  Default: nothing is done. 

%% Recognized options in varargin 
%
% 'pdf'       Produce basic pdf output, e.g.: example.pdf
%
% 'latex'     Produce combined pdf & LaTeX files, e.g.:
%             example.tex and example-inc.pdf


%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

  pdf = false;   % If true, save plot to pdf file           
latex = false;   % If true, save plot to combined pdf & latex files

%% PARSE OPTIONAL ARGUMENT LIST

for k = 1:length(varargin);
    switch lower(varargin{k}) % varargin is a "cell array"
      case {'pdf'}
        pdf=true;
      case {'latex'}
        latex=true;
      case {'no'}
        pdf=false; latex=false;
      otherwise
        error(['error in function savepdf: ',...
               'option %s not recognized.\n'],varargin{k});
    end
end

%% CHECK IF A FIGURE IS AVAILABLE

h = findall(0,'type','figure');
if isempty(h)
    error(['function savepdf: can not produce pdf file ' ...
           'because no figure is open.'])
end

%% OPTIONS

if (pdf)
    pdf_file=changefiletype(file_name,'pdf');
    print(pdf_file,'-dpdf','-landscape','-fillpage');% Basic pdf output
    fprintf('Plot of window %s saved in %s\n',get(gcf,'name'),pdf_file);
end

if (latex)
    pdf_file=changefiletype(file_name,'pdf');
    tex_file=changefiletype(file_name,'tex');
    inc_file=sprintf('%s-inc.pdf',pdf_file(1:end-4));
    print(pdf_file,'-dpdflatex');% Combined pdf & LaTeX files
			         % .tex file to be input in LaTeX document
    fprintf('Plot of window %s saved in %s and %s (LaTeX)\n',...
            get(gcf,'name'),tex_file,inc_file);
end

end % End of function savepdf