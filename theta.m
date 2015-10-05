function t = theta(A)

    n = length(A);
    blk{1,1} = 's';      blk{1,2} = n;
    e = nnz(A)/2;

    As = cell(1,1 + e);
    constraint_ind = 1;

    As{constraint_ind} = speye(n,n);
    constraint_ind = constraint_ind + 1;

    for i = 2:n
        for j = 1:i-1
            if(A(i,j)==1)
                As{constraint_ind} = sparse(j,i,1,n,n);
                constraint_ind = constraint_ind + 1;
            end
        end
    end

    At(1) = svec(blk(1,:),As,1);
    C{1} = -ones(n,n);
    b = [1;zeros(e,1)];
%     [obj,~] = sqlp(blk,At,C,b);
    [~,obj,~,~,~] = evalc('sqlp(blk,At,C,b)');
    t = -obj(1);
%     disp(X{1})
%     disp(['Independent number is ' num2str(-alpha)])
%     disp(['SDP relaxation theta_-delta is ' num2str(-obj(1))]);
    
end