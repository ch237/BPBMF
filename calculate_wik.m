function wik=calculate_wik(U,id,lambda,N)
    if iscell(id)
    else
        id_row=size(id,1);
        id=mat2cell(id,id_row,ones(1,3));
    end
    
    K=size(U,2);
    wik=zeros(N,K);
    
    for i=1:length(id{1})
%         if id{1}~=id{2}
            wik(id{1}(i),:)=wik(id{1}(i),:)+U(id{2}(i),:)*lambda{id{3}(i)};
    end
end

