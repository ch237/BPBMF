% function [U lambda]=RelationOnlineGibbs(xi,id,xi_tensor,para)
    a=1e-2;
    U=sampleDirMat(a*ones(1,para.N),para.K);
    U=U';
    
    Nnon0=length(xi);
    Train=floor(para.TrainFrac*Nnon0);
    idall=randperm(Nnon0);
    idtest=cell(1,length(id));
    idtrain=cell(1,length(id));
    for k=1:length(id)
        idtest{k} = id{k}(idall(Train+1:end)); 
        idtrain{k} = id{k}(idall(1:Train)); 
    end

    zeroN=round(0.1*nnz(~xi_tensor));%ceil((Nnon0-Train)*numel(xi_tensor)/length(xi));
    zeroID=cell(1,3);     
    [zeroID{1},zeroID{2},zeroID{3}] = ind2sub(size(xi_tensor),find(xi_tensor== 0));
    zeroIDid=randperm(length(zeroID{1}));
    zeroIDid=zeroIDid(1:zeroN);
    for k=1:3
        idtest{k}=[idtest{k};zeroID{k}(zeroIDid)];
    end
    xitest=[ones(Nnon0-Train,1);zeros(zeroN,1)];

    gamma0=1e-0*ones(para.R,1);
    c0=1e-0*ones(para.R,1);%gamrnd(1,1);
    beta=1e-0*ones(para.R,1);%1;%gamrnd(1,1);
    rk=cell(1,para.R);
    
    for r=1:para.R
        rk{r}=gamrnd(gamma0(r)/para.K,1/c0(r),para.K,1);
    end
    e0=1;
    f0=1;
    epsilon=zeros(1,para.R);
    for r=1:para.R
        epsilon(r)=gamrnd(e0,1/f0);
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
%         lambda{r}=gamrnd(rkrkt{r},1/beta);
        lambda{r}=rkrkt{r};
    end 
    
    Np=floor(Train*para.batchfrac);
    idselect=cell(1,3);
    t0=0;
    
    step=5;
    sampleN=(para.itermax-para.burnin)/step;
    
for iter=1:para.itermax   
    tic
%     iter
    gam_t = (iter+t0)^(-0.5);
    idid=randperm(Train,Np);
    idid=idall(idid);
    xiselect = xi(idid);
    for k=1:3
        idselect{k} = id{1, k}(idid);      
    end
    
    zetair = multinomialpara(U,idselect,lambda);
    zetai=sum(zetair,2);
    Ezetai=zetai./(1-exp(-zetai)+eps);
    Eratio=Ezetai./(zetai+eps);
    xi_k1k2=zetair.*repmat(Eratio,1,para.K^2);
    [xi_k,x_kk]=calculate_xik(xi_k1k2,idselect,para.N,para.R);
    
    
    if iter==1
        xi_k=xi_k./para.batchfrac;
        for r=1:para.R
            x_kk{r}=x_kk{r}./para.batchfrac;
        end

    else       
        xi_k=(1-gam_t)*xi_k_old+gam_t*xi_k./para.batchfrac;
        for r=1:para.R
            x_kk{r}=(1-gam_t)*x_kk_old{r}+gam_t*x_kk{r}./para.batchfrac;
        end
    end
    
    if iter==1
        dirpara=a+xi_k;
    end
    idunique=unique([idselect{1};idselect{2}]);
    for k=1:para.K
        dirpara(idunique,k)=a+xi_k(idunique,k);
        U(:,k)=dirpara(:,k)./sum(dirpara(:,k));%sampleDirMat(a+xi_k(:,k)',1)';    
    end

    %%
    if iter==1
        theta_kk=U(idselect{1},:)'*U(idselect{2},:)./para.batchfrac;
    else
        theta_kk = (1-gam_t)*theta_kk_old + gam_t*U(idselect{1},:)'*U(idselect{2},:)./para.batchfrac;
    end

    
    
    L_KK=zeros(para.K,para.K);
    temp_p_tilde_k=zeros(para.K,1);
    for r=1:para.R
        p_kk_prime_one_minus = beta(r)./(beta(r)+theta_kk);
        for k=randperm(para.K)%1:para.K
            R_KK=rk{r}';
            R_KK(k)=epsilon(r);
            L_KK(k,:) = CRT_sum_mex_matrix(sparse(round(x_kk{r}(k,:))),rk{r}(k)*R_KK);
            temp_p_tilde_k(k) = -sum(R_KK.*log(max(p_kk_prime_one_minus(k,:), realmin)));
%             rk{r}(k)=gamrnd(gamma0(r)/para.K+sum(L_KK(k,:)),1/(c0(r)+temp_p_tilde_k(k)));
            rk{r}(k)=(gamma0(r)/para.K+sum(L_KK(k,:)))./(c0(r)+temp_p_tilde_k(k));
        end

        ell = sum(CRT_sum_mex_matrix(sparse(round(diag(x_kk{r}))'),epsilon(r)*rk{r}'));
%         epsilon(r) = randg(ell+e0)/(f0-sum(rk{r}.*log(max(diag(p_kk_prime_one_minus), realmin))));
        epsilon(r) = (ell+e0)/(f0-sum(rk{r}.*log(max(diag(p_kk_prime_one_minus), realmin))));
        epsilon(r)=max(epsilon(r),1e-6);
    end
    
    for r=1:para.R
        rkrkt{r}=rk{r}*rk{r}';
        for k=1:para.K            
            rkrkt{r}(k,k)=epsilon(r)*rk{r}(k);
        end  
    end
    for r=1:para.R
        lambda{r}=(rkrkt{r}+x_kk{r})./(beta(r)+theta_kk);
%         lambda{r}=gamrnd(rkrkt{r}+x_kk{r},1./(beta(r)+theta_kk));
    end

    
    xi_k_old=xi_k;
    x_kk_old=x_kk;
    theta_kk_old=theta_kk;
    time_periter=toc
    
    [mae rmse mse auc_test auc_pr_test zetaitest]=relation_eval(xitest,idtest,U,lambda);
%     llikevec(iter)=llike;
    rmsevec(iter)=rmse;maevec(iter)=mae;msevec(iter)=mse;auc(iter)=auc_test;auc_pr(iter)=auc_pr_test;
    fprintf('iter= %d; mae=%f, mse=%f, rmse=%f, auc=%f, auc_pr=%f\n', iter, mae, mse, rmse, auc_test,auc_pr_test);
    
    figure(1),plot(auc(1:iter));
    xlabel('Iterations');
    ylabel('AUC');
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