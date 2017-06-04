function [mae rmse mse auc_test auc_test2 zetai]=relation_eval(xi,id_temp,U,lambda)
% K=length(U);
if iscell(id_temp)
    id=id_temp;
else    
id=cell(1,3);
    for k=1:3
        id{1,k}=id_temp(:,k);
    end
end 


R=size(lambda{1},1);
Np=length(id{1,1});
zetai=zeros(Np,1);
% zetair=zeros(Np,R^2);
for i=1:Np
%     tempmat=(U(id{1,1}(i),:)'*U(id{1,2}(i),:)).*lambda{id{1,3}(i)};
%     zetair(i,:)=reshape(tempmat,[1,R^2]);
%     if all(zetair(i,:))==1
%         zetair(i,:)=1/R^2;
% %         zetair(i,:)=1e-12;
%     end
    zetai(i)=U(id{1,1}(i),:)*lambda{id{1,3}(i)}*U(id{1,2}(i),:)';
end
% zetai=sum(zetair,2);
% llike = sum(xi.*log(1-exp(-zetai)) - (1-xi).*zetai);
rmse=norm(xi-(1-exp(-zetai)),'fro')/(sqrt(size(xi,1) * size(xi,2)));
mae=sum(abs(xi-(1-exp(-zetai))))/Np;
mse=rmse^2;
auc_test = compute_AUC(xi,1-exp(-zetai),ones(1,length(xi)));

%[prec, tpr, fpr, thresh] = prec_rec(1-exp(-zetai), xi,  'numThresh',3000);
[prec, tpr, fpr, thresh] = prec_rec(1-exp(-zetai), xi);
% auc_test2=0;
auc_test2 = trapz([0;tpr],[1;prec]);
% acc=sum(round(1-exp(-zetai))==xi)/length(xi);
% acc=sum(((1-exp(-zetai))>0.00004)==xi)/length(xi);
% acc=sum((xi-(1-exp(-zetai))).^2)/length(xi);
%     xi_pred = double((1-exp(-zetai)) >= 0.01);
%     figure(2),hist(1-exp(-zetai),50);
%     drawnow
% zetai=1-exp(-zetai);
% [precision,recall,T,AUC] = perfcurve(xi,zetai,1,'XCrit','prec');
% plot(recall,precision);
% xlabel('Recall');
% ylabel('Precision'); 
% drawnow;
% figure(100),plot(1-exp(-zetai))