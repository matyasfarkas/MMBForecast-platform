function disp_identification(pdraws, idemodel, idemoments, name, advanced)

% Copyright (C) 2008-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global options_

if nargin < 5 || isempty(advanced)
    advanced=0;
end

[SampleSize, npar] = size(pdraws);
% jok = 0;
% jokP = 0;
% jokJ = 0;
% jokPJ = 0;

% for j=1:npar,
%     %     if any(idemodel.ind(j,:)==0),
%     %         pno = 100*length(find(idemodel.ind(j,:)==0))/SampleSize;
%     %         disp(['Parameter ',name{j},' is not identified in the model for ',num2str(pno),'% of MC runs!' ])
%     %         disp(' ')
%     %     end
%     %     if any(idemoments.ind(j,:)==0),
%     %         pno = 100*length(find(idemoments.ind(j,:)==0))/SampleSize;
%     %         disp(['Parameter ',name{j},' is not identified by J moments for ',num2str(pno),'% of MC runs!' ])
%     %         disp(' ')
%     %     end
%     if any(idemodel.ind(j,:)==1),
%         iok = find(idemodel.ind(j,:)==1);
%         jok = jok+1;
%         kok(jok) = j;
%         mmin(jok,1) = min(idemodel.Mco(j,iok));
%         mmean(jok,1) = mean(idemodel.Mco(j,iok));
%         mmax(jok,1) = max(idemodel.Mco(j,iok));
%         [ipmax, jpmax] = find(abs(squeeze(idemodel.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
%         if ~isempty(ipmax)
%             jokP = jokP+1;
%             kokP(jokP) = j;
%             ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
%             [N,X]=hist(ipmax,[1:npar]);
%             jpM(jokP)={find(N)};
%             NPM(jokP)={N(find(N))./SampleSize.*100};
%             pmeanM(jokP)={mean(squeeze(idemodel.Pco(j,find(N),iok))')};
%             pminM(jokP)={min(squeeze(idemodel.Pco(j,find(N),iok))')};
%             pmaxM(jokP)={max(squeeze(idemodel.Pco(j,find(N),iok))')};
%         end
%     end
%     if any(idemoments.ind(j,:)==1),
%         iok = find(idemoments.ind(j,:)==1);
%         jokJ = jokJ+1;
%         kokJ(jokJ) = j;
%         mminJ(jokJ,1) = min(idemoments.Mco(j,iok));
%         mmeanJ(jokJ,1) = mean(idemoments.Mco(j,iok));
%         mmaxJ(jokJ,1) = max(idemoments.Mco(j,iok));
%         [ipmax, jpmax] = find(abs(squeeze(idemoments.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
%         if ~isempty(ipmax)
%             jokPJ = jokPJ+1;
%             kokPJ(jokPJ) = j;
%             ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
%             [N,X]=hist(ipmax,[1:npar]);
%             jpJ(jokPJ)={find(N)};
%             NPJ(jokPJ)={N(find(N))./SampleSize.*100};
%             pmeanJ(jokPJ)={mean(squeeze(idemoments.Pco(j,find(N),iok))')};
%             pminJ(jokPJ)={min(squeeze(idemoments.Pco(j,find(N),iok))')};
%             pmaxJ(jokPJ)={max(squeeze(idemoments.Pco(j,find(N),iok))')};
%         end
%     end
% end

disp(['  ']),


no_warning_message_display=1;

if any(idemodel.ino) || any(any(idemodel.ind0==0)) || any(any(idemodel.jweak_pair))
    no_warning_message_display=0;
    disp('WARNING !!!')
    if SampleSize>1
        disp(['The rank of H (model) is deficient for ', num2str(length(find(idemodel.ino))),' out of ',int2str(SampleSize),' MC runs!'  ]),
    else
        disp(['The rank of H (model) is deficient!'  ]),
    end
    skipline()
    for j=1:npar
        if any(idemodel.ind0(:,j)==0)
            pno = 100*length(find(idemodel.ind0(:,j)==0))/SampleSize;
            if SampleSize>1
                disp(['    ',name{j},' is not identified in the model for ',num2str(pno),'% of MC runs!' ])
            else
                disp(['    ',name{j},' is not identified in the model!' ])
            end
            disp(['    [dJ/d(',name{j},')=0 for all tau elements in the model solution!]' ])
        end
    end
    npairs=size(idemodel.jweak_pair,2);
    jmap_pair=dyn_unvech(1:npairs);
    jstore=[];
    skipline()
    for j=1:npairs
        iweak = length(find(idemodel.jweak_pair(:,j)));
        if iweak
            [jx,jy]=find(jmap_pair==j);
            jstore=[jstore jx(1) jy(1)];
            if SampleSize > 1
                disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear (with tol = 1.e-10) for ',num2str((iweak)/SampleSize*100),'% of MC runs!' ])
            else
                disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear (with tol = 1.e-10) !' ])
            end
        end

    end
    skipline()
    for j=1:npar
        iweak = length(find(idemodel.jweak(:,j)));
        if iweak && ~ismember(j,jstore)
            %         disp('WARNING !!!')
            %         disp(['Model derivatives of parameter ',name{j},' are multi-collinear (with tol = 1.e-10) for ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
            if SampleSize>1
                disp([name{j},' is collinear w.r.t. all other params ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
            else
                disp([name{j},' is collinear w.r.t. all other params!' ])
            end
        end
    end
    %         if npar>(j+1),
    %             [ipair, jpair] = find(squeeze(idemodel.Pco(j,j+1:end,:))'>(1-1.e-10));
    %         else
    %             [ipair, jpair] = find(squeeze(idemodel.Pco(j,j+1:end,:))>(1-1.e-10));
    %         end
    %         if ~isempty(jpair),
    %             for jx=j+1:npar,
    %                 ixp = find(jx==(jpair+j));
    %                 if ~isempty(ixp)
    %                     if SampleSize > 1,
    %                         disp(['    [',name{j},',',name{jx},'] are PAIRWISE collinear (with tol = 1.e-10) for ',num2str(length(ixp)/SampleSize*100),'% of MC runs!' ])
    %                     else
    %                         disp(['    [',name{j},',',name{jx},'] are PAIRWISE collinear (with tol = 1.e-10)!' ])
    %                     end
    %                 end
    %             end
    %         end
end

if no_warning_message_display
    disp(['All parameters are identified in the model (rank of H).' ]),
    skipline()
end

no_warning_message_display = 1;

if any(idemoments.ino) || any(any(idemoments.ind0==0)) || any(any(idemoments.jweak_pair))
    no_warning_message_display = 0;
    skipline()
    disp('WARNING !!!')
    if SampleSize > 1
        disp(['The rank of J (moments) is deficient for ', num2str(length(find(idemoments.ino))),' out of ',int2str(SampleSize),' MC runs!'  ]),
    else
        disp(['The rank of J (moments) is deficient!'  ]),
    end
    %     disp('WARNING !!!')
    %     disp(['The rank of J (moments) is deficient for ', num2str(length(find(idemoments.ino))/SampleSize*100),'% of MC runs!'  ]),
    %     indno=[];
    %     for j=1:SampleSize, indno=[indno;idemoments.indno{j}]; end
    %     freqno = mean(indno)*100;
    %     ifreq=find(freqno);
    %     disp('MOMENT RANK FAILURE DUE TO COLLINEARITY OF PARAMETERS:');
    skipline()
    for j=1:npar
        if any(idemoments.ind0(:,j)==0)
            pno = 100*length(find(idemoments.ind0(:,j)==0))/SampleSize;
            if SampleSize > 1
                disp(['    ',name{j},' is not identified by J moments for ',num2str(pno),'% of MC runs!' ])
            else
                disp(['    ',name{j},' is not identified by J moments!' ])
            end
            disp(['    [dJ/d(',name{j},')=0 for all J moments!]' ])
        end
    end
    skipline()
    npairs=size(idemoments.jweak_pair,2);
    jmap_pair=dyn_unvech(1:npairs);
    jstore=[];
    for j=1:npairs
        iweak = length(find(idemoments.jweak_pair(:,j)));
        if iweak
            [jx,jy]=find(jmap_pair==j);
            jstore=[jstore'  jx(1) jy(1)]';
            if SampleSize > 1
                disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear (with tol = 1.e-10) for ',num2str((iweak)/SampleSize*100),'% of MC runs!' ])
            else
                disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear (with tol = 1.e-10) !' ])
            end
        end

    end
    skipline()
    for j=1:npar
        iweak = length(find(idemoments.jweak(:,j)));
        if iweak && ~ismember(j,jstore)
            %             disp('WARNING !!!')
            %             disp(['Moment derivatives of parameter ',name{j},' are multi-collinear (with tol = 1.e-10) for ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
            if SampleSize > 1
                disp([name{j},' is collinear w.r.t. all other params ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
            else
                disp([name{j},' is collinear w.r.t. all other params!' ])
            end
        end
    end
    %             if npar>(j+1),
    %                 [ipair, jpair] = find(squeeze(idemoments.Pco(j,j+1:end,:))'>(1-1.e-10));
    %             else
    %                 [ipair, jpair] = find(squeeze(idemoments.Pco(j,j+1:end,:))>(1-1.e-10));
    %             end
    %             if ~isempty(jpair),
    %                 for jx=j+1:npar,
    %                     ixp = find(jx==(jpair+j));
    %                     if ~isempty(ixp)
    %                         if SampleSize > 1
    %                             disp(['    [',name{j},',',name{jx},'] are PAIRWISE collinear (with tol = 1.e-10) for ',num2str(length(ixp)/SampleSize*100),'% of MC runs!' ])
    %                         else
    %                             disp(['    [',name{j},',',name{jx},'] are PAIRWISE collinear (with tol = 1.e-10) !' ])
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %     end
end
if no_warning_message_display
    skipline()
    disp(['All parameters are identified by J moments (rank of J)' ]),
    skipline()
end

% if ~ options_.noprint && advanced,
%     disp('Press KEY to continue with identification analysis')
%     pause;
%     dyntable('Multi collinearity in the model:',char('param','min','mean','max'), ...
%              char(name(kok)),[mmin, mmean, mmax],10,10,6);
%     disp(' ')
%     dyntable('Multi collinearity for moments in J:',char('param','min','mean','max'), ...
%              char(name(kokJ)),[mminJ, mmeanJ, mmaxJ],10,10,6);
%     disp(' ')
% end


% if advanced && (~options_.noprint),
%     for j=1:length(kokP),
%         dyntable([name{kokP(j)},' pairwise correlations in the model'],char(' ','min','mean','max'), ...
%                  char(name{jpM{j}}),[pminM{j}' pmeanM{j}' pmaxM{j}'],10,10,3);
%     end
%
%     for j=1:length(kokPJ),
%         dyntable([name{kokPJ(j)},' pairwise correlations in J moments'],char(' ','min','mean','max'), ...
%                  char(name{jpJ{j}}),[pminJ{j}' pmeanJ{j}' pmaxJ{j}'],10,10,3);
%     end
% end
% disp(' ')

% identificaton patterns
if SampleSize==1 && advanced
    skipline()
    disp('Press ENTER to print advanced diagnostics'), pause(5),
    for  j=1:size(idemoments.cosnJ,2)
        pax=NaN(npar,npar);
        fprintf('\n')
        disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
        fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
        for i=1:npar
            namx='';
            for in=1:j
                dumpindx = idemoments.pars{i,j}(in);
                if isnan(dumpindx)
                    namx=[namx ' ' sprintf('%-15s','--')];
                else
                    namx=[namx ' ' sprintf('%-15s',name{dumpindx})];
                    pax(i,dumpindx)=idemoments.cosnJ(i,j);
                end
            end
            fprintf('%-15s [%s] %14.7f\n',name{i},namx,idemoments.cosnJ(i,j))
        end
    end
end
