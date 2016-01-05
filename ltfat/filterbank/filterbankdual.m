function gdout=filterbankdual(g,a,L);
%FILTERBANKDUAL  Dual filters
%   Usage:  gd=filterbankdual(g,a);
%           gd=filterbankdual(g,a,L);
%
%   FILTERBANKDUAL(g,a) computes the canonical dual filters of g for a
%   channel subsampling rate of a (hop-size).
%
%   The input and output format of the filters g are described in the
%   help of FILTERBANK.
%
%   FILTERBANKDUAL(g,a,L) computes canonical dual filters for a system
%   of length L. If L is not specified, the shortest possible
%   transform length is choosen.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the IFILTERBANK function.
%
%   See also: filterbank, ufilterbank, ifilterbank
%
%   Url: http://ltfat.sourceforge.net/doc/filterbank/filterbankdual.php

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

complainif_notenoughargs(nargin,2,'FILTERBANKDUAL');

if nargin<3
   L = [];
end

[g,info] = filterbankwin(g,a,L,'normal'); 
M=info.M;

Lsut = filterbanklength(L,a);
if (~isempty(L)) && (L~=filterbanklength(L,a))
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts. Next suitable L is %i; obtainable by ',...
           'filterbanklength(a,L).'],upper(mfilename),Lsut);
end;


% Prioritize painless over uniform algorithm
if info.isuniform && info.ispainless
    info.isuniform = 0;
end

if info.isuniform
    
  % Uniform filterbank, use polyphase representation
  if isempty(L)
      error('%s: You need to specify L.',upper(mfilename));
  end;

  a=a(1);
  
  % G1 is done this way just so that we can determine the data type.
  G1=comp_transferfunction(g{1},L);
  thisclass=assert_classname(G1);
  G=zeros(L,M,thisclass);
  G(:,1)=G1;
  for ii=2:M
    G(:,ii)=comp_transferfunction(g{ii},L);
  end;
  
  N=L/a;
  
  gd=zeros(N,M,thisclass);
  
  for w=0:N-1
    idx = mod(w-(0:a-1)*N,L)+1;
    H = G(idx,:);
    
    H=pinv(H)';
    
    gd(idx,:)=H;
  end;
  
  gd=ifft(gd)*a;
  
  if isreal(g)
    gd=real(gd);
  end;
  
  gdout=cell(1,M);
  for m=1:M
    gdout{m}=struct('h',cast(gd(:,m),thisclass),'offset',0);
  end;
  
else

    if info.ispainless
                
        if isempty(L)
            error('%s: You need to specify L.',upper(mfilename));
        end;
        
        gdout = comp_painlessfilterbank(g,info.a,L,'dual',0);

%        F=comp_filterbankresponse(g,info.a,L,0);
        
%         gdout=cell(1,M);
%         for m=1:M
%             thisgd=struct();
%             if isfield(g{m},'H')
%                H=circshift(comp_transferfunction(g{m},L)./F,-g{m}.foff);
%                thisgd.H=H(1:numel(g{m}.H));
%                thisgd.foff=g{m}.foff;
%                thisgd.realonly=0;
%                thisgd.delay=0;
%             elseif isfield(g{m},'h')
%                H=comp_transferfunction(g{m},L)./F; 
%                thisgd.h = ifft(H);
%                thisgd.offset = 0;
%             end
%             
%             gdout{m}=thisgd;
%         end;
        
    else
        error(['%s: The canonical dual frame of this system is not a ' ...
               'filterbank. You must either call an iterative ' ...
               'method to perform the desired inverstion or transform ',...
               'or transform the filterbank to uniform one. Please see ' ...
               'FRANAITER or FRSYNITER for the former and ',...
               'NONU2UFILTERBANK for the latter case.'],upper(mfilename));        

    end;
    
end;

