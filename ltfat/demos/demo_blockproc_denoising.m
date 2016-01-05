function demo_blockproc_denoising(source,varargin) 
%DEMO_BLOCKPROC_DENOISING Variable coefficients thresholding
%   Usage: demo_blockproc_denoising('gspi.wav')
%
%   For additional help call DEMO_BLOCKPROC_DENOISING without arguments.
%
%   The present demo allows you to set the coefficient threshold during the
%   playback using the control panel.
% 
%
%   Url: http://ltfat.sourceforge.net/doc/demos/demo_blockproc_denoising.php

% Copyright (C) 2005-2014 Peter L. Soendergaard <soender@users.sourceforge.net>.
% This file is part of LTFAT version 1.4.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if demo_blockproc_header(mfilename,nargin)
   return;
end

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               {'Thr','Treshold',0,0.1,0,1000}
               });
    
% Buffer length
bufLen = 1024;
% Number of frequency channels
M = 1000;

% Setup blocktream
fs=block(source,varargin{:},'loadind',p,'L',bufLen);

% Window length in ms
winLenms = 20; %floor(fs*winLenms/1e3)
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',40,M);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');


flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   
  % Obtain parameters from the control panel
  gain = 10^(p.getParam('GdB')/20); % dB -> val
  thres = p.getParam('Thr');
  %bufLen = floor(p.getParam('bufLen'));

  % Read block of length bufLen
  [f,flag] = blockread();
  % Apply analysis frame
  c = blockana(Fa, f*gain); 
  % Plot
  % blockplot(fobj,F,c);
  % Apply thresholding
  c = thresh(c,thres,'soft');
  % Apply synthesis frame
  fhat = real(blocksyn(Fs, c, size(f,1)));
  % Play the block
  %fhat = f;
  blockplay(fhat);
end
blockdone(p);

