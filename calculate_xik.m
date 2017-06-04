function [xi_k,x_kk]=calculate_xik(xi_k1k2,id,N,R)
    if iscell(id)
    else
        id_row=size(id,1);
        id=mat2cell(id,id_row,ones(1,3));
    end
    
    K=sqrt(length(xi_k1k2(1,:)));
    xi_k=zeros(N,K);
    x_kk=cell(1,R);
    
    for r=1:R
        x_kk{r}=zeros(K,K);
    end
    
    for i=1:length(id{1})
%         if id{1}~=id{2}
            xi_k1k2mat=reshape(xi_k1k2(i,:),[K, K]);
            xi_k(id{1}(i),:)=xi_k(id{1}(i),:)+sum(xi_k1k2mat,2)';
%             xi_k(id{2}(i),:)=xi_k(id{2}(i),:)+sum(xi_k1k2mat,1);
            x_kk{id{3}(i)}=x_kk{id{3}(i)}+xi_k1k2mat;
%         end
    end
end

