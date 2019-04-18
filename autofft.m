function [spectrum, freq] = autofft(xs, ts, fftset)
% AUTOFFT Evaluates a frequency spectrum using wFFT algorithm 
% v1.22a (build 26. 11. 2018)
%
% Lubos Smolik, University of West Bohemia (carlist{at}ntis.zcu.cz)
%
% [spectrum, freq] = autofft(xs, ts)
% [SPECTRUM, freq] = autofft(XS, ts)
% [...] = autofft(..., fftset)
%
% [spectrum, freq] = autofft(xs, ts) returns column vectors spectrum and
%       freq, which contain values of discrete autospectrum and their
%       respective frequencies. xs is a vector of samples and ts is a
%       vector of time stamps (s) or the sampling frequency (Hz).
% [SPECTRUM, freq] = autofft(XS, ts) returns array SPECTRUM which contains
%       column vectors of discrete autospectra and vector freq. XS is an
%       array consisting of column vectors of samples and ts is a vector of
%       time stamps (s) or the sampling frequency (Hz), which is valid for
%       column vectors in XS.
% [...] = autofft(..., fftset) performs Fast Fourier analysis in accordance
%       with the input specified in structured variable fftset.
%       
%       Construction of fftset:
%       fftset = struct('param', 'value', ...);
%
%       List of fftset parameters:
%        - 'nwin' - number of samples in each window  (default: length(xs))
%        - 'twin' - window length in (s)              (default: fs / nwin)
%           - Note: If both nwin and twin are specified, twin is ignored.
%        - 'overlap' - overlaping of two successive time windows in (%)
%                                                     (default: 50)
%        - 'lowpass' - frequency for low-pass filter  (default: fs/2)
%        - 'window' - time weighting window
%          - 'b' - Blackmann-Harris
%          - 'f' - flat-top
%          - 'h' - Hann
%          - 'm' - Hamming
%          - 'k' - Kaiser-Bessel, beta = 0.5
%          - 'kA.A' - Kaiser-Bessel, beta = A.A
%          - 'u' - uniform                            (default)
%        - 'averaging' - spectral averaging
%          - 'energy','rms' - energy averaging
%          - 'linear','lin' - linear averaging        (default)
%          - 'max','peak'   - maximum peak hold aveargaing
%          - 'min'          - minimum peak hold aveargaing
%          - 'none'         - no averaging (returns spectrum for all segments)
%        - 'jw' - j-omega weighting
%          - '1/jw2'        - double integration
%          - '1/jw'         - integrace
%          - '1'            - as input                (default)
%          - 'jw'           - derivation
%          - 'jw2'          - double derivation
%        - 'unit' - spectral unit
%          - 'pow'          - autospectrum             (default)
%          - 'rms'          - linear spectrum with rms magnitude
%          - 'pk'           - linear spectrum with 0-peak magnitude
%          - 'pp'           - linear spectrum with peak-peak magnitude
%          - 'psd'          - power spectral density
%          - 'rsd','rmspsd' - root mean square of power spectral density 
%
% Changelist
% v1.22a- Error occuring during peak hold averaging has been fixed.
% v1.22 - The Kaiser-Bessel window parameter (beta) can now be specified.
%       - Autospectrum is now properly square of RMS rather than 0-Pk.
%       - Relations for evaluation of PSD and RMSPSD now consider the noise
%         power bandwidth of the used time weighting window. 
% v1.21 - New function - low-pass filtering
%       - New types of averaging - no averaging     ('none')
%                                - energy averaging ('energy' or 'rms')
%                                - minimum value    ('min')
%       - Options for 'unit' and 'peak' has been merged (into 'unit').
%       - Relations for evaluation of PSD amd RMSPSD have been fixed.
%       - Dealing with a content at the Nyquist frequency has been fixed.
% v1.2  - Input parameters are now specified in a structured variable.
%       - v1.2 is not compatible with v1.12 and older versions!
% v1.12 - Input can now be an array.
% v1.11 - Performance optimization
%       - Handling of input vectors with the even number of samples has
%         been fixed.
% v1.1  - Handling of non-uniform time weighting windows has been fixed.
%       - Parameters nwin and overlap can now be skipped by user.
%
%% This code is published under BSD-2-Clause License
%
% Copyright (c) 2017-2018, Lubos Smolik
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%% nargin check
if nargin < 2
    error('Not enough input arguments.');
elseif nargin > 3
    error('Too many input arguments.');
end
%
%% Convert row vectors to column vectors if needed
if size(xs, 1) == 1         % samples
    xs = xs(:);                     
end
if size(ts(:), 1) == 1      % sampling frequency
	fs = ts;
else
    fs = 1 / (ts(2) - ts(1)); 
end
%
%% Specify default fftset
defset = struct('nwin', size(xs, 1), ...
                'twin', size(xs, 1) / fs, ...
                'overlap', 50, ...
                'lowpass', fs/2, ...
                'window', 'u', ...
                'averaging', 'lin', ...
                'jw', '1', ...
                'unit', 'pow');
deffields = fieldnames(defset);
%
%% Set analyser parameters
if nargin == 2  % use default fftset
    fftset = defset;
else            % use user-defined fftset  
    % Check whether there is user-defined 'nwin' or 'twin' parameter
    if isfield(fftset, 'nwin')
        fftset.twin = fftset.nwin / fs;
    elseif isfield(fftset, 'twin')
        fftset.nwin = round(fftset.twin * fs);
    end
    % Set unspecified parameters to default
    for i = 1:numel(deffields)        
        if isfield(fftset, deffields{i}) == 0
            fftset.(deffields{i}) = defset.(deffields{i});
        end
    end
end
% Generate frequency vector
freq = (fs * (0:(fftset.nwin/2)) / fftset.nwin)';
% Set allowed frequencies for the low-pass filtering
freq = freq(freq <= fftset.lowpass);
% Set handling of the last spectral line
if freq(end) == fs/2
    sh = 1; % Nyquist frequency is at the last spectral line
else
	sh = 0; % Nyquist frequency is not at the last spectral line
end
% Set number of overlaping samples
fftset.overlap = round(fftset.nwin * fftset.overlap / 100);
% Set indices for the signal segmentation
imax = floor((size(xs, 1)-fftset.overlap) / (fftset.nwin-fftset.overlap));
                                                % number of windows
ind = zeros(imax,2);                            % matrix of indices
ni = 1;                                         % pointer
for i = 1:imax                                  % cycle through windows
    ind(i,1) = ni; 
    ni = ni + fftset.nwin - 1;
    ind(i,2) = ni;
    ni = ni - fftset.overlap + 1;
end
% Generate the time weighting window
[fftset.window, fftset.noiseband] = timeweight(fftset.window, fftset.nwin, fs);
% Set constant for the jw weigthing
switch lower(fftset.jw)
    case '1/jw2'
        fftset.jw = - 1 ./ (4 * pi^2 * freq.^2);
    case '1/jw'
        fftset.jw = 1 ./ (2i * pi * freq);
    case 'jw'
        fftset.jw = 2i * pi * freq;
    case 'jw2'
        fftset.jw = - 4 * pi^2 * freq.^2;
    otherwise
        fftset.jw = 1;
end
%
%% Spectral analyser
spectrum = zeros(size(freq, 1), size(xs, 2));
%
for i = 1:size(xs, 2)
    % Preallocate an array for temporary spectra
    spectra = zeros(size(freq, 1), imax);
    % Frequency analysis of individual segments
    for j = 1:imax
        % Fast Fourier transformation of the weighted segment
        fftSamp = fft(fftset.window .* xs(ind(j,1):ind(j,2), i), ...
                      fftset.nwin) / fftset.nwin;
        % Application of the jw weigthing and the low-pass filtering
        spectra(:, j) = fftset.jw .* fftSamp(1:size(freq, 1));
        % Evaluation of spectral unit
        switch lower(fftset.unit)
            case 'rms'           % Linear spectrum with rms magnitude
                spectra(2:end-sh,j) = (2/sqrt(2))*abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = abs(spectra(end-sh+1:end,j))/sqrt(2);
            case 'pk'            % Linear spectrum with 0-peak magnitude
                spectra(2:end-sh,j) = 2 * abs(spectra(2:end-sh,j));
            case 'pp'            % Linear spectrum with peak-peak magnitude
                spectra(2:end-sh,j) = 4 * abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = 2 * abs(spectra(end-sh+1:end,j));
            case 'psd'           % Power spectral density
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = (1 / fftset.noiseband) * spectra(:,j) ...
                                       .* conj(spectra(:,j));
            case {'rsd','rmspsd'}% Root mean square of PSD 
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = sqrt((1 / fftset.noiseband) * ...
                                       spectra(:,j) .* conj(spectra(:,j)));                
            otherwise            % Autospectrum
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = spectra(:,j) .* conj(spectra(:,j));
        end
    end
    % Spectral averaging
    switch lower(fftset.averaging)
        case {'energy', 'rms'}   % Energy averaging
            spectrum(:, i) = sqrt(sum(spectra.^2, 2) ./ imax);
        case {'max', 'peak'}     % Maximum peak hold averaging
            spectrum(:, i) = max(spectra, [], 2);
        case 'min'              % Minimum peak hold averaging
            spectrum(:, i) = min(spectra, [], 2);
        case 'none'             % No averaging
            if size(xs, 2) == 1
                spectrum = spectra;
            else
                spectrum{i} = spectra;
            end
        otherwise               % Linear averaging
            spectrum(:, i) = mean(spectra, 2);
    end
    % Remove imaginary part (residue due to spectra .* conj(spectra))
    spectrum = real(spectrum);
end
% End of main fucntion
if nargout == 0
    plot(freq, spectrum)
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction weightsig generates the time weighting window
%
% Input
%   - sym - symbol for the time weighting window
%   - n   - length of the time weighting window (samples)
%   - fs  - sampling frequency (Hz) 
%
% Output
%   - window    - the time weighting window
%   - noiseband - the noise power bandwidth of the time weighting window
%  
function [window, noiseband] = timeweight(sym, n, fs)
    % Generate the specified time weighting window
    switch sym(1)
        case 'b'    % Blackmann-Harris
            window = blackmanharris(n);
        case 'f'    % flat-top
            window = flattopwin(n);
        case 'h'    % Hann
            window = hann(n);
        case 'k'    % Kaiser-Bessel
            if length(sym) == 1
                % window with default beta = 0.5
                window = kaiser(n, 0.5);
            else
                % window with user specified beta
                window = kaiser(n, str2double(sym(2:end)));
            end
        case 'm'    % Hamming
            window = hamming(n);
        otherwise   % uniform
            window = rectwin(n);
    end
    % Adjust magnitude of the generated time weighting window
    window = window / mean(window);
    % Calculate the noise power bandwidth in Hz
    noiseband = enbw(window, fs);
end
