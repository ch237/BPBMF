function zetair = multinomialpara_rev(U,id,G,eta)
Nnon0=length(id{1,1});
K=size(G{1},1);
M=size(eta,1);
zetair=1e-12+zeros(Nnon0,K^2*M);%repmat(ones(1,K^2*M),Nnon0,1);
% zetair=repmat(ones(1,K^2*M),Nnon0,1);
tempmat=zeros(K,K,M);
for i=1:Nnon0
    utu=(U(id{1,1}(i),:)'*U(id{1,2}(i),:));
    for m=1:M
        tempmat(:,:,m)=utu.*(eta(m,id{3}(i)).*G{m});
    end
    zetair(i,:)=reshape(tempmat,[1,K^2*M]);
end



% function zetair = multinomialpara_rev(U,id,lambda)
% if iscell(id)
% else
%     id_row=size(id,1);
%     id=mat2cell(id,id_row,ones(1,3));
% end
% Nnon0=length(id{1,1});
% K=size(lambda{1},1);
% zetair=repmat(ones(1,K^2),Nnon0,1);
% % zetair_normal=zeros(size(zetair));
% for i=1:Nnon0
%     tempmat=(U(id{1,1}(i),:)'*U(id{1,2}(i),:)).*lambda{id{1,3}(i)};
%     zetair(i,:)=reshape(tempmat,[1,K^2]);
%     if all(~zetair(i,:))==1
%         zetair(i,:)=1e-12;%1e-1;%1/K^2;
%     end
% end