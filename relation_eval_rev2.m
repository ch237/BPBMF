
function [auc_test auc_pr zetai] = relation_eval_rev2(xi,id,U,G,eta)
if iscell(id)
else
    id_row=size(id,1);
    id=mat2cell(id,id_row,ones(1,3));
end
Nnon0=length(id{1,1});
K=size(G{1},1);
[M R]=size(eta);

% compute all the lambda's
for r=1:R
    lambda{r} = zeros(K);
    for m=1:M
        lambda{r} = lambda{r} + (eta(m,r)*G{m});
    end
end

Np=length(id{1,1});
zetai=zeros(Np,1);
% zetair=zeros(Np,K^2);
for i=1:Np
%     tempmat=(U(id{1,1}(i),:)'*U(id{1,2}(i),:)).*lambda{id{1,3}(i)};
%     zetair(i,:)=reshape(tempmat,[1,K^2]);
    zetai(i)=U(id{1,1}(i),:)*lambda{id{1,3}(i)}*U(id{1,2}(i),:)';
%     if zetai(i)==0
%         zetai(i)=M;
%     end
%     if all(zetair(i,:))==1
%         zetair(i,:)=1/K^2;
% %         zetair(i,:)=1e-12;
%     end
%     if all(~zetair(i,:))==1
%         zetair(i,:)=1e-12;
% %         zetair(i,:)=1/(K^2*M);
%     end
end
% zetai=sum(zetair,2);
prob=1-exp(-zetai);
auc_test = compute_AUC(xi,prob,ones(1,length(xi)));

[prec, tpr, fpr, thresh] = prec_rec(prob, xi);
auc_pr = trapz([0;tpr],[1;prec]);

% auc_pr= auc_pr(1-prob', xi');
