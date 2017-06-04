% function [U lambda]=RelationGibbs(xi,id,xi_tensor,para)
    a=1e-2;
    U=sampleDirMat(a*ones(1,para.N),para.K);
    U=U';
%     rng(0);


    Train=length(idtrain{1});
    xitrain=ones(Train,1);
    
    gam_a=1e0;
    gam_b=1e0;
    gamma0=gam_a*ones(para.M,1);
    c0=gam_b*ones(para.M,1);%gamrnd(1,1);
    beta=gam_b*ones(para.M,1);%1;%gamrnd(1,1);
    h0=gam_a*1e-0;%h0=1e-2;q0=1e-3;
    q0=gam_b*1e-0;
    
    rk=cell(1,para.M);
    for m=1:para.M
        rk{m}=gamrnd(gamma0(m)/para.K,1/c0(m),para.K,1);
    end
    e0=gam_a;
    f0=gam_b;

    epsilon=zeros(1,para.M);
    for m=1:para.M
        epsilon(m)=gamrnd(e0,1/f0);
    end
    rkrkt=cell(1,para.M);
    for m=1:para.M
        rkrkt{m}=rk{m}*rk{m}';
        for i=1:para.K
            rkrkt{m}(i,i)=epsilon(m)*rk{m}(i);
        end
    end
    G=cell(1,para.M);
    for m=1:para.M
        G{m}=randg(rkrkt{m})./beta(m);
    end
    eta=randg(h0,[para.M,para.R])./q0;    
    
step=5;
sampleN=(para.itermax-para.burnin)/step;
L_KK=cell(1,para.M);

for m=1:para.M
    L_KK{m}=zeros(para.K,para.K);
end



Np=floor(Train*para.batchfrac);
idselect=cell(1,3);
t0=10;

step=5;
sampleN=(para.itermax-para.burnin)/step;
    

tic;

for iter=1:para.itermax  
    
    gam_t =(iter+t0)^(-0.5);
    
    idid=randperm(Train,Np);
    xiselect = xitrain(idid);
    for k=1:3
        idselect{k} = idtrain{1, k}(idid);      
    end
    
    zetair = multinomialpara_rev(U,idselect,G,eta);
    zetai=sum(zetair,2);
    Ezetai=zetai./(1-exp(-zetai)+eps);
    Eratio=Ezetai./(zetai+eps);
%     xi_count=truncated_Poisson_rnd_1(xiselect',zetai);
%     xi_k1k2m=mnrnd(xi_count,zetair./repmat(sum(zetair,2),1,para.K^2*para.M));
    xi_k1k2m=zetair.*repmat(Eratio,1,para.K^2*para.M);
    [xi_k,x_kkm,x_rm]=calculate_xik_rev(xi_k1k2m,idselect,para.N,para.R,para.M,para.K);
    
    if iter==1
        xi_k=xi_k./para.batchfrac;
        x_kkm=x_kkm./para.batchfrac;
        for r=1:para.R
            x_rm{r}=x_rm{r}./para.batchfrac;
        end
    else
        xi_k=(1-gam_t)*xi_k_old+gam_t*xi_k./para.batchfrac;
        x_kkm=(1-gam_t)*x_kkm_old+gam_t*x_kkm./para.batchfrac;
        for r=1:para.R
            x_rm{r}=(1-gam_t)*x_rm_old{r}+gam_t*x_rm{r}./para.batchfrac;
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
    
    if iter==1
        theta_kk=U(idselect{1},:)'*U(idselect{2},:)./para.batchfrac;%%%%%%%%%%%%%%%%%%%%%
    else
        theta_kk = (1-gam_t)*theta_kk_old + gam_t*U(idselect{1},:)'*U(idselect{2},:)./para.batchfrac;
    end
    
    
    
    eta_sum=sum(eta,2);    
%     L_KK=zeros(para.K,para.K);
    temp_p_tilde_k=zeros(para.K,1);
    for m=1:para.M       
        p_kk_prime_one_minus = beta(m)./(beta(m)+theta_kk.*eta_sum(m));
        for k=1:para.K%randperm(para.K)%1:para.K
            R_KK=rk{m}';
            R_KK(k)=epsilon(m);
            L_KK{m}(k,:) = CRT_sum_mex_matrix(sparse(round(reshape(x_kkm(k,:,m),[1,para.K]))),rk{m}(k)*R_KK);
            temp_p_tilde_k(k) = -sum(R_KK.*log(max(p_kk_prime_one_minus(k,:), realmin)));
            rk{m}(k)=(gamma0(m)/para.K+sum(L_KK{m}(k,:)))./(c0(m)+temp_p_tilde_k(k));
        end
        
        ell = sum(CRT_sum_mex_matrix(sparse(round(diag(x_kkm(:,:,m))')),epsilon(m)*rk{m}'));
%         epsilon(m) = randg(ell+e0)/(f0-sum(rk{m}.*log(max(diag(p_kk_prime_one_minus), realmin))));
        epsilon(m) = (ell+e0)/(f0-sum(rk{m}.*log(max(diag(p_kk_prime_one_minus), realmin))));
        epsilon(m)=max(epsilon(m),1e-6);
    end
    
    for m=1:para.M
        rkrkt{m}=rk{m}*rk{m}';
        for k=1:para.K            
            rkrkt{m}(k,k)=epsilon(m)*rk{m}(k);
        end  
    end
    
    for m=1:para.M
        theta_G_sum=sum(theta_kk(:).*G{m}(:));
        for r=1:para.R
%             eta(m,r)=randg(h0+x_rm{r}(m))/(q0+theta_G_sum);
            eta(m,r)=(h0+x_rm{r}(m))/(q0+theta_G_sum);
        end
    end    
    
%     eta_sum=sum(eta,2);
       
    for m=1:para.M
%         G{m}=randg(rkrkt{m}+reshape(x_kkm(:,:,m),[para.K,para.K]))./(eta_sum(m)+beta(m));
        G{m}=(rkrkt{m}+reshape(x_kkm(:,:,m),[para.K,para.K]))./(eta_sum(m)+beta(m));
    end

    xi_k_old=xi_k;
    x_kkm_old=x_kkm;
    x_rm_old=x_rm;
    theta_kk_old=theta_kk;
    
    
    if iter==1 
            time_trace(iter) = toc;
        else
            time_trace(iter) = time_trace(iter-1) + toc;
    end
    
    
    stepEva=5;
    if mod(iter,stepEva)==1
        [auc_test,auc_pr_test]=relation_eval_rev2(xitest,idtest,U,G,eta);
        auc(iter)=auc_test;
        aucpr(iter)=auc_pr_test;
        fprintf('iter= %d;auc=%f;auc_pr=%f;\n',iter,auc_test,auc_pr_test); 
    else
        fprintf('iter= %d;\n',iter); 
    end
    
    plot(1:stepEva:iter,auc(1:stepEva:iter));
    xlabel('Iterations');
    ylabel('AUC');
    drawnow;
    
    tic;
    
%     if iter>para.burnin & mod(iter-para.burnin,step)==0
%             if iter==(para.burnin+step)
%                 zetaifinal=zetaitest/sampleN;
%             else
%                 zetaifinal=zetaifinal+zetaitest/sampleN;
%             end
%     end
end 
% auc_sam = compute_AUC(xitest,zetaifinal,ones(1,length(xitest)));
% fprintf('Final test AUC=%f\n', auc_sam);