%% Generate random graph

n = 5;
A = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 1 0 0 1 0];
% A = graph(n);
disp([num2str(nnz(A)/2) ' edges in the graph'])

%% Original Lovasz theta

% blk{1,1} = 's';      blk{1,2} = n;
% e = nnz(A)/2;
% 
% As = cell(1,1 + e);
% constraint_ind = 1;
% 
% As{constraint_ind} = speye(n,n);
% constraint_ind = constraint_ind + 1;
% 
% for i = 2:n
%     for j = 1:i-1
%         if(A(i,j)==1)
%             As{constraint_ind} = sparse(j,i,1,n,n);
%             constraint_ind = constraint_ind + 1;
%         end
%     end
% end
% 
% At(1) = svec(blk(1,:),As,1);
% C{1} = -ones(n,n);
% b = [1;zeros(e,1)];
% [obj,X] = sqlp(blk,At,C,b);
% % [~,obj,X,y,Z] = evalc('sqlp(blk,At,C,b)');
% disp(X{1})
% disp(['Independent number is ' num2str(-alpha)])
% disp(['SDP relaxation theta_-delta is ' num2str(-obj(1))]);


%% Equivalalent Lovasz Theta

% blk{1,1} = 's';      blk{1,2} = n+1;
% e = nnz(A)/2;
% 
% 
% As = cell(1,1 + n + e);
% constraint_ind = 1;
% 
% % X_00 = 1;
% As{constraint_ind} = sparse(1,1,1,n+1,n+1);
% constraint_ind = constraint_ind + 1;
% 
% s = sqrt(2);
% % X_0i = X_ii
% for i = 2:n+1
%     As{constraint_ind} = sparse([1 i],[i i],[-1/2 1],n+1,n+1); % no clue why divide by 2
%     constraint_ind = constraint_ind + 1;
% end
% 
% % X_ij = 0 for ij in E
% for i = 2:n
%     for j = 1:i-1
%         if(A(j,i)==1)
%             As{constraint_ind} = sparse(j+1,i+1,1,n+1,n+1);
%             constraint_ind = constraint_ind + 1;
%         end
%     end
% end
% 
% At(1) = svec(blk(1,:),As,1);
% C{1} = -speye(n+1);  C{1}(1) = 0;
% b = [1;zeros(n,1);zeros(e,1)];
% [obj,X] = sqlp(blk,At,C,b);
% % [~,obj,X,y,Z] = evalc('sqlp(blk,At,C,b)');
% disp(X{1})
% disp(['Independent number is ' num2str(-alpha)])
% disp(['SDP relaxation theta_-delta is ' num2str(-obj(1))]);

%% Lovasz Theta Minus

% n2 = nchoosek(n+1,2); % number of inequality constraints
% 
% blk{1,1} = 's';      blk{1,2} = n;
% blk{2,1} = 'l';      blk{2,2} = n2;
% 
% e = nnz(A)/2;
% 
% As = cell(1,1 + e + n2);
% constraint_ind = 1;
% 
% As{constraint_ind} = speye(n,n);
% constraint_ind = constraint_ind + 1;
% 
% % X_ij = 0 for ij in E
% for i = 2:n
%     for j = 1:i-1
%         if(A(i,j)==1)
%             As{constraint_ind} = sparse(j,i,1,n,n);
%             constraint_ind = constraint_ind + 1;
%         end
%     end
% end
% 
% % -X_ij <= 0
% for i = 1:n
%     for j = 1:i
%         As{constraint_ind} = sparse(j,i,-1,n,n);
%         constraint_ind = constraint_ind + 1;
%     end
% end
% 
% At(1) = svec(blk(1,:),As,1);
% At{2,1} = [sparse(n2,1+e), speye(n2)];
% C{1,1} = -ones(n,n);
% C{2,1} = zeros(n2,1);
% b = [1;zeros(e,1);zeros(n2,1)];
% % [obj,X] = sqlp(blk,At,C,b);
% [~,obj,X,y,Z] = evalc('sqlp(blk,At,C,b)');
% % disp(['Independent number is ' num2str(-alpha)])
% disp(['SDP relaxation theta_-delta is ' num2str(-obj(1))]);

%% Lovasz Theta Minus Delta

% n2 = nchoosek(n+1,2) + 2 * nchoosek(n,2);% + 3 * nchoosek(n,3); % number of inequality constraints
n2 = 0;
blk{1,1} = 's';      blk{1,2} = n;

e = nnz(A)/2;

As = cell(0);
constraint_ind = 1;

As{constraint_ind} = speye(n,n);
constraint_ind = constraint_ind + 1;

% X_ij = 0 for ij in E
for i = 2:n
    for j = 1:i-1
        if(A(i,j)==1)
            As{constraint_ind} = sparse(j,i,1,n,n);
            constraint_ind = constraint_ind + 1;
        end
    end
end

% -X_ij <= 0
for i = 1:n
    for j = 1:i
        As{constraint_ind} = sparse(j,i,-1,n,n);
        constraint_ind = constraint_ind + 1;
        n2 = n2 + 1;
    end
end

% X_ji - X_ii <= 0
for i = 1:n
    for j = 1:i - 1
        As{constraint_ind} = sparse([j j],[i j],[1/2 -1],n,n);
        constraint_ind = constraint_ind + 1;
        n2 = n2+1;
        As{constraint_ind} = sparse([j i],[i i],[1/2 -1],n,n);
        constraint_ind = constraint_ind + 1;
        n2 = n2+1;
    end
end

% X_jk + X_ik - Xij - X_kk <= 0
for k = 1:n
    for i = 2:n
        if (i==k)
            continue
        else
            for j = 1:i-1
                if(j==k)
                    continue
                else
                    As{constraint_ind} = sparse([min(j,k) min(i,k) j k],[max(j,k) max(i,k) i k],[1/2 1/2 -1/2 -1],n,n);
                    constraint_ind = constraint_ind + 1;
                    n2 = n2+1;
                end
            end
        end
    end
end
blk{2,1} = 'l';      blk{2,2} = n2;

At(1) = svec(blk(1,:),As,1);
At{2,1} = [sparse(n2,1+e), speye(n2)];
C{1,1} = -ones(n,n);
C{2,1} = zeros(n2,1);
b = [1;zeros(e,1);zeros(n2,1)];
% [obj,X] = sqlp(blk,At,C,b);
[~,obj,X,y,Z] = evalc('sqlp(blk,At,C,b)');
% disp(['Independent number is ' num2str(-alpha)])
disp(['SDP relaxation theta_-delta is ' num2str(-obj(1))]);