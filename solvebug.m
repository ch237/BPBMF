% clc;

% load wordnetresult50;%wn50
% U=U./repmat(sum(U,2),1,size(U,2));
% zetair=multinomialpara(U,idselect,lambda);
eshowid=enameid('/m/03zb6t');%/m/042xh
rshowid=rnameid('/sports/sports_team/roster./sports/sports_team_roster/position');%/influence/influence_node/influenced_by;
% eshowid=4449;
% simscore=zeros(para.R,para.N);
% for r=1:para.R
%     simscore(r,:)=U(eshowid,:)*lambda{r}*U';
% end
simscore=U(eshowid,:)*lambda{rshowid}*U';

relatedEntN=20;

fprintf('--------------------------------------------------------------\t\n');
[score tailshowid]=sort(simscore,'descend');
fprintf('%s---%s---\t\n',eidname(num2str(eshowid)),ridname(num2str(rshowid)));
for i=1:relatedEntN
    fprintf(' %s\t\n',eidname(num2str(tailshowid(i))));
end
fprintf('-----------------------in test-------------------------------\t\n');
for i=1:length(idtest{1})
    if idtest{1}(i)==eshowid && idtest{3}(i)==rshowid
       fprintf('%s\t\n',eidname(num2str(idtest{2}(i))));
    end
end
fprintf('-----------------------in train-------------------------------\t\n');
for i=1:length(idtrain{1})
    if idtrain{1}(i)==eshowid && idtrain{3}(i)==rshowid
       fprintf('%s\t\n',eidname(num2str(idtrain{2}(i))));
    end
end


% rshowid=10;
% for rshowid=1:para.R
%     fprintf('--------------------------------------------------------------\t\n');
%     [score tailshowid]=sort(simscore(rshowid,:),'descend');
%     fprintf('%s---%s---\t\n',eidname(num2str(eshowid)),ridname(num2str(rshowid)));
%     for i=1:relatedEntN
%         fprintf(' %s\t\n',eidname(num2str(tailshowid(i))));
%     end
%     fprintf('--------------------------------------------------------------\t\n')

% end