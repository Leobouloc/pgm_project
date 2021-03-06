function L=comp_framelength_tensor(F,Ls);
%COMP_FRAMELENGTH_TENSOR  Helper function for the Tensor frame
%  
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_framelength_tensor.php

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
error(['For the tensor product frame, call framelength for the ', ...
       'individual subframes']);
