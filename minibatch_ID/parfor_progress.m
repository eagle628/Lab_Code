function percent = parfor_progress(N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

narginchk(0, 1);

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar

if N > 0
    f = fopen('parfor_progress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    Ns = sprintf('%d', N);
    if nargout == 0
        disp(['  0% (', blanks(length(Ns)-1), '0/', Ns ') ' '[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    f = fopen('parfor_progress.txt', 'r');
    if f == -1
        return
    end
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    delete('parfor_progress.txt');
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        maxitr = sprintf('%d', progress(1));
        pr = sprintf('%d', (length(progress)-1));
        pr = [blanks(length(maxitr)-length(pr)) pr];
        ratio = [pr '/' maxitr];
        disp([repmat(char(8), 1, (w+9+2*length(maxitr)+5)), char(10), perc, ' (' ratio, ') [', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
%     fprintf('Test');
else
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        maxitr = sprintf('%d', progress(1));
        pr = sprintf('%d', (length(progress)-1));
        pr = [blanks(length(maxitr)-length(pr)) pr];
        ratio = [pr '/' maxitr];
        disp([repmat(char(8), 1, (w+9+2*length(maxitr)+5)), char(10), perc, ' (' ratio, ') [', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end