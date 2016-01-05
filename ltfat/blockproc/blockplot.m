function cola=blockplot(p,F,c,cola)
%BLOCKPLOT Plot block coefficients
%   Usage: blockplot(p,F,c);
%
%   Input parameters:
%         p     : JAVA object of the class net.sourceforge.ltfat.SpectFrame.
%         F     : Frame object.
%         c     : Block coefficients.
%         cola  : (Optional) overlap from previous block.
%
%   Output parameters:
%         cola  : Overlap to the next block.
%
%   BLOCKPLOT(p,F,c) appends the block coefficients c to the running 
%   coefficient plot in p. The coefficients must have been obtained by
%   c=blockana(F,...). The format of c is changed to a rectangular 
%   layout according to the type of F.
%
%   BLOCKPLOT(p,[],c) does the same, but expects c to be already
%   formatted matrix of real numbers. The matrix dimensions are not
%   restricted, but it will be shrinked or expanded to a vertical
%   strip in the sliding image.
%
%   cola=BLOCKPLOT(p,F,c,cola) does the same, but adds cola to the 
%   first respective coefficients in c and returns last coefficients from
%   c. This is only relevant for the sliced window blocking approach.
%
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockplot.php

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

complainif_notenoughargs(nargin,3,'BLOCKPLOT');
complainif_notvalidframeobj(F,'BLOCKPLOT');


if ~isempty(F)
    ctf = framecoef2tf(F,c(:,1));
    
    if size(c,2)>1
        error('%s: Only one channel input is supported.',upper(mfilename));
    end

    if strcmp(F.blockalg,'sliced')
       % DO the coefficient overlapping or cropping
       %ctf = ctf(:,floor(end*3/8):floor(end*5/8)+1);

       if nargin>3 
          olLen = ceil(size(ctf,2)/2);
          if isempty(cola)
             cola = zeros(size(ctf,1),olLen,class(ctf));
          end

          ctf(:,1:olLen) = ctf(:,1:olLen) + cola;
          cola = ctf(:,end+1-olLen:end);
          ctf = ctf(:,1:olLen);
       end
    end
    
    ctf = abs(ctf);
else
    if ~isreal(c)
        error('%s: Complex values are not supported',upper(mfilename));
    end
    ctf = c;
end


if isoctave
   % The JAVA 2D-array handling is row-major
   ctf = cast(ctf,'double').';
   javaMethod('append',p,ctf(:),size(ctf,2),size(ctf,1));
else
   % Matlab casts correctly
   ctf = cast(ctf,'single');
   javaMethod('append',p,ctf);
end



