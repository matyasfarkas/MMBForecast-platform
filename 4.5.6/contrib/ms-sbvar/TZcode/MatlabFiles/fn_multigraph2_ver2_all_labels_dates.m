function scaleout = fn_multigraph2_ver2_all_labels_dates(dates,imf1,imf2,...
                         xlab,ylab,tstring,XTick,YTickIndx,scaleIndx,legend_inp, nrowg,ncolg)
%scaleout = fn_multigraph2_ver2_all_labels_dates(dates,imf1,imf2,...
%                         xlab,ylab,tstring,XTick,YTickIndx,scaleIndx,nrowg,ncolg)
%
%Stacking two sets of impulse responses in one graph.  See fn_multigraph2_ver2_all_labels.m.
%   dates:  time horizon in the graphics,
%   imf1, imf2: row: the same dimension as the dates (in the graphics),
%               column: "nrow "variables (row in the graphics),
%               3rd dim: across "ncol" different situations (column in the graphics, but may be overwritten by ncolg).
%               If the 3rd D is only 1, then we can specifiy nrowg and ncolg.
%   xlab:   x-axis labels on the top.  Use [] if no label is needed.
%   ylab:   y-axis labels on the left.  Use [] if no label is needed.
%   tstring:  string for time (e.g., month or quarter) -- x-axis labels at the bottom.  Use [] if no label is needed.
%   XTick: ticks on the x-axis with grids on.  If [], no x-ticks and no grids.
%   YTickIndx:  1: enable YTick; 0: disable.  Often set to zero to get a better picture.
%   scaleIndx:  1: enable scale along Y-axis; 0: disable
%   legend_inp: legend inputs.  If [], no legend.
%   nrowg: number of rows in the graph (valid when the 3rd D is one)
%   ncolg: number of columns in the graph (valid when the 3rd D is one)
%
%  See fn_multigraph2_ver2_all_labels.m, fn_multigraph2_ver2, imrgraph, imcerrgraph, imrerrgraph

%
% Copyright (C) 1997-2012 Tao Zha
%
% This free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% If you did not received a copy of the GNU General Public License
% with this software, see <http://www.gnu.org/licenses/>.
%

nstp = size(imf1,1);
nrow = size(imf1,2);
ncol = size(imf1,3);

if (length(dates) ~= nstp)
   error('.../fn_multigraph2_ver2_all_labels_dates: dates must have the dimension as the row number of imf1 and imf2!')
end

if (nargin < 11)
  nrowg = nrow;
  ncolg = ncol;
end
if (nrowg*ncolg < nrow*ncol)
   error('fn_multigraph2_ver2_all_labels_dates.m: nrowg*ncolg must be greater or equal to nrow*ncol')
end


temp1=zeros(ncol,1);
temp2=zeros(ncol,1);
maxval=zeros(nrow,1);
minval=zeros(nrow,1);
for i = 1:nrow
   for j = 1:ncol
      %jnk1=max(firsth(:,i,j));
      %jnk2=max(firstl(:,i,j));
      %jnk3=max(firsth1(:,i,j));
      jnk4=max(imf2(:,i,j));
      jnk5=max(imf1(:,i,j));

      temp1(j)=max([jnk4 jnk5]);
      %
      %jnk1=min(firstl(:,i,j));
      %jnk2=min(firsth(:,i,j));
      %jnk3=min(firstl1(:,i,j));
      jnk4=min(imf2(:,i,j));
      jnk5=min(imf1(:,i,j));

      temp2(j)=min([jnk4 jnk5]);
	end
   maxval(i)=max(temp1);
   minval(i)=min(temp2);
end

scaleout = [maxval(:) minval(:)];

%--------------
%  Column j: Shock 1 to N; Row i: Responses to
%-------------
%figure
rowlabel = 1;
for i = 1:nrow      % column: from top to bottom
   columnlabel = 1;

   if minval(i)<0
      if maxval(i)<=0
         yt=[minval(i) 0];
		else
         yt=[minval(i) 0 maxval(i)];
		end
   else % (minval(i) >=0)
      if maxval(i) > 0
         yt=[0 maxval(i)];
		else % (identically zero responses)
			yt=[-1 0 1];
		end
	end


   scale=[dates(1) dates(end) minval(i) maxval(i)];
   for j = 1:ncol        % row: from left to right
      k1=(i-1)*ncol+j;
      subplot(nrowg,ncolg,k1)
      plot(dates,imf1(:,i,j),'k-',dates,imf2(:,i,j),'k--')
                %t,[firstl(:,i,j) firsth(:,i,j)],':');
      grid;
      if (~isempty(legend_inp))
         legend(legend_inp)
      end
      if scaleIndx
         axis(scale);
      end
      %set(gca,'YLim',[minval(i) maxval(i)])
      %
		set(gca,'XTick',XTick)
      if YTickIndx
         set(gca,'YTick',yt)
      end

      %  if i<nrow
      %    set(gca,'XTickLabelMode','manual','XTickLabel',[])
      %  end
      %set(gca,'XTickLabel',' ');
      if (scaleIndx) && (j>1)
         set(gca,'YTickLabel',' ');
      end
      %  if rowlabel == 1
      %     %title(['x' num2str(j)])
      %     %title(eval(['x' num2str(j)]))
      %     title(char(xlab(j)))
      %  end
      if (~isempty(xlab))
         title(char(xlab(i)))
      end
      %
      %  if columnlabel == 1
      %     %ylabel(['x' num2str(i)])
      %     %ylabel(eval(['x' num2str(i)]))
      %     ylabel(char(ylab(i)))
      %  end
      if (~isempty(ylab))
         ylabel(char(ylab(i)))
      end
      if (i==nrow) && (~isempty(tstring))
         xlabel(tstring)
      end
      columnlabel = 0;
   end
   rowlabel = 0;
end
