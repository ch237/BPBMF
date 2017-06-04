% function [U lambda]=RelationOnlineGibbs(xi,id,xi_tensor,para)
    rng(0);
    a=1e-2;
    U=sampleDirMat(a*ones(1,para.N),para.K);
    U=U';
    
    Train=length(idtrain{1});
    xitrain=ones(Train,1);
    
    
    gamconst1=1e-0; 
    gamconst2=1e-0; 
    gamma0=gamconst1*ones(para.R,1);
    c0=gamconst2*ones(para.R,1);
    beta=gamconst2*ones(para.R,1);
    e0=gamconst1;
    f0=gamconst1;
    rk=cell(1,para.R);
    
    for r=1:para.R
        rk{r}=randg(gamma0(r)/para.K,[para.K,1])./c0(r);%gamrnd(gamma0(r)/para.K,1/c0(r),para.K,1);
    end
    
    epsilon=zeros(1,para.R);
    for r=1:para.R
        epsilon(r)=randg(e0)/f0;%gamrnd(e0,1/f0);
    end
    rkrkt=cell(1,para.R);
    for r=1:para.R
        rkrkt{r}=rk{r}*rk{r}';
        for i=1:para.K
            rkrkt{r}(i,i)=epsilon(r)*rk{r}(i);
        end
    end
    
    lambda=cell(1,para.R);%zeros(para.K,para,K,para.R);
    for r=1:para.R
          lambda{r}=randg(rkrkt{r})./beta(r);
    end 
    
    Np=floor(Train*para.batchfrac);
    idselect=cell(1,3);
    t0=0;
    
    step=1;
    sampleN=(para.itermax-para.burnin)/step;
    iterperpass=ceil(1/para.batchfrac);
    sN=1;
for iter=1:para.itermax   
%     iter
    if para.batchfrac==1
        gam_t=1;
    else
        gam_t = (iter+t0)^(-0.5);
    end
    idid=randperm(Train,Np);
    xiselect = xitrain(idid);
    for k=1:3
        idselect{k} = idtrain{1, k}(idid);      
    end
    zetair=multinomialpara(U,idselect,lambda);
    zetai=sum(zetair,2);
    if iter==1
        zetai=ones(size(zetai));
    end
%     zetai=zetai./median(zetai);
    xi_count=truncated_Poisson_rnd_1(xiselect',zetai);
%     multrnd_histc(xi_count,zetair)
    xi_k1k2=mnrnd(xi_count,zetair./repmat(sum(zetair,2),1,para.K^2));
    [xi_k,x_kk]=calculate_xik(xi_k1k2,idselect,para.N,para.R);
    if iter==1%iter<=3*iterperpass%mod(iter,iterperpass)==1%iter==1
%         if iter==1
%         else
%             xi_k=xi_k_old+xi_k;
%             for r=1:para.R
%                 x_kk{r}=x_kk_old{r}+x_kk{r};
%             end
%             
%         end
% %         xi_k=xi_k;%./para.batchfrac;
        xi_k=xi_k./para.batchfrac;
        for r=1:para.R
% %             x_kk{r}=x_kk{r};%./para.batchfrac;
            x_kk{r}=x_kk{r}./para.batchfrac;
        end
    else
%         xi_k=xi_k_old+xi_k;%(1-gam_t)*xi_k_old+gam_t*xi_k./para.batchfrac;
        xi_k=(1-gam_t)*xi_k_old+gam_t*xi_k;%./para.batchfrac;
        for r=1:para.R
%             x_kk{r}=x_kk_old{r}+x_kk{r};
            x_kk{r}=(1-gam_t)*x_kk_old{r}+gam_t*x_kk{r};%./para.batchfrac;
        end
    end
    U_temp=zeros(size(U));
    lambda_temp=lambda;
    rk_temp=rk;
    for r=1:para.R
        lambda_temp{r}=zeros(size(lambda{r}));
        rk_temp{r}=zeros(size(rk{r}));
    end
    c0_temp=zeros(size(c0));
    epsilon_temp=zeros(size(epsilon));
    for s=1:sN
        
        
    for k=1:para.K
%         U(:,k)=(a+a+xi_k(:,k))./sum(a+a+xi_k(:,k));  
        U(:,k)=sampleDirMat(a+xi_k(:,k)',1)';
    end
    U_temp=U_temp+U./sN;
    %%
    theta_kk=U(idselect{1},:)'*(U(idselect{2},:))./para.batchfrac^2;
    
    for r=1:para.R
        rkrkt{r}=rk{r}*rk{r}';
        for k=1:para.K            
            rkrkt{r}(k,k)=epsilon(r)*rk{r}(k);
        end  
    end
    for r=1:para.R
        lambda{r}=randg(rkrkt{r}+x_kk{r})./(beta(r)+theta_kk);
        lambda_temp{r}= lambda_temp{r}+lambda{r}./sN;
    end
    
    L_KK=zeros(para.K,para.K);
    temp_p_tilde_k=zeros(para.K,1);
%     p_kk_prime_one_minus = beta./(beta+theta_kk);
    for r=1:para.R
        p_kk_prime_one_minus = beta(r)./(beta(r)+theta_kk);
        for k=randperm(para.K)%1:para.K
            R_KK=rk{r}';
            R_KK(k)=epsilon(r);
            L_KK(k,:) = CRT_sum_mex_matrix(sparse(round(x_kk{r}(k,:))),rk{r}(k)*R_KK);
            temp_p_tilde_k(k) = -sum(R_KK.*log(max(p_kk_prime_one_minus(k,:), realmin)));
%             rk{r}(k)=gamrnd(gamma0(r)/para.K+sum(L_KK(k,:)),1/(c0(r)+temp_p_tilde_k(k)));
            rk{r}(k)=randg(gamma0(r)/para.K+sum(L_KK(k,:)))./(c0(r)+temp_p_tilde_k(k));
            rk_temp{r}(k)=rk_temp{r}(k)+rk{r}(k)./sN;
        end

        ell = sum(CRT_sum_mex_matrix(sparse(round(diag(x_kk{r}))'),epsilon(r)*rk{r}'));
        epsilon(r) = randg(ell+e0)/(f0-sum(rk{r}.*log(max(diag(p_kk_prime_one_minus), realmin))));
        epsilon(r)=max(epsilon(r),1e-6);
        epsilon_temp(r)=epsilon_temp(r)+epsilon(r)./sN;

%         ell_tilde = CRT_sum_mex(round(sum(L_KK,2)),gamma0(r)/para.K);
%         sum_p_tilde_k_one_minus = -sum(log(c0(r)./(c0(r)+temp_p_tilde_k) ));
%         gamma0new = randg(e0 + ell_tilde)./(f0 + 1/para.K*sum_p_tilde_k_one_minus);
%         AcceptProb1 = CalAcceptProb1(rk{r},c0(r),gamma0(r),gamma0new,ell_tilde,1/para.K*sum_p_tilde_k_one_minus,para.K);
%         if AcceptProb1>rand(1)
% %             AcceptProb1
%             gamma0(r)=gamma0new;
%         end 
%         c0(r)=randg(gamconst1+gamma0(r))./(gamconst2+sum(rk{r}));
% %         c0_temp(r)=c0_temp(r)+c0(r)/sN;
% % 
%         R_KK = rk{r}*rk{r}';
%         R_KK(sparse(1:para.K,1:para.K,true)) = epsilon(r)*rk{r};        
%         beta(r)=randg(gamconst1+sum(sum(R_KK)))./(gamconst2+sum(sum(lambda{r})));
%         beta(r)=randg(gamconst1+para.K^2)./(gamconst2+sum(sum(lambda{r})));
    end
    end
%     c0=c0_temp;
    epsilon=epsilon_temp;
    U=U_temp;
    lambda=lambda_temp;
    rk=rk_temp;
    
    
    xi_k_old=xi_k;
    x_kk_old=x_kk;
    theta_kk_old=theta_kk;
    
    
    [mae rmse mse auc_test auc_pr_test zetaitest]=relation_eval(xitest,idtest,U,lambda);
    rmsevec(iter)=rmse;maevec(iter)=mae;msevec(iter)=mse;auc(iter)=auc_test;auc_pr(iter)=auc_pr_test;
    fprintf('iter= %d; mae=%f, mse=%f, rmse=%f, auc=%f, auc_pr=%f,maxxicount=%f;\n', iter, mae, mse, rmse, auc_test,auc_pr_test,max(xi_count));
    figure(1),plot(auc(1:iter));
    xlabel('Iterations');
    ylabel('AUC');

    
%     figure(2),bar(zetai)
%     figure(2),imagesc(U)
    drawnow;
    
    if iter>para.burnin & mod(iter-para.burnin,step)==0
        if iter==(para.burnin+step)
            zetaifinal=zetaitest/sampleN;
        else
            zetaifinal=zetaifinal+zetaitest/sampleN;
        end
    end
    
    
%     figure(100),bar(epsilon);drawnow;
end 

auc_sam = compute_AUC(xitest,zetaifinal,ones(1,length(xitest)));
fprintf('Final test AUC=%f\n', auc_sam);
    
% end