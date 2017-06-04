function [xi_k,x_kkm,x_rm]=calculate_xik_rev(xi_k1k2m,id,N,R,M,K)
%     if iscell(id)
%     else
%         id_row=size(id,1);
%         id=mat2cell(id,id_row,ones(1,3));
%     end
    
%     K=sqrt(length(xi_k1k2m(1,:)));
    xi_k=zeros(N,K);
    x_rm=cell(1,R);
    x_kkm=zeros(K,K,M);
    
    for r=1:R
        x_rm{r}=zeros(1,M);
    end
    
    for i=1:length(id{1})
%         if id{1}~=id{2}
            xi_k1k2mtensor=reshape(xi_k1k2m(i,:),[K, K, M]);
            xi_k1k2mat=sum(xi_k1k2mtensor,3);
            xi_rm_vec=reshape(sum(sum(xi_k1k2mtensor,1),2),[1,M]);
            xi_k(id{1}(i),:)=xi_k(id{1}(i),:)+sum(xi_k1k2mat,2)';
%             xi_k(id{2}(i),:)=xi_k(id{2}(i),:)+sum(xi_k1k2mtensor,1);
            x_kkm=x_kkm+xi_k1k2mtensor;
            x_rm{id{3}(i)}=x_rm{id{3}(i)}+xi_rm_vec;
%         end
    end
end

