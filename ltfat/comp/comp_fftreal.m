function f=comp_fftreal(f)
%COMP_FFTREAL  Compute an FFTREAL
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_fftreal.php

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
  
N=size(f,1);
N2=floor(N/2)+1;

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
f=fft(f,N,1);
f=f(1:N2,:);


