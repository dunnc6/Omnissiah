function [opt,iter] = SimplexVanderBei(A,b,c)
  % A port of Vanderbei's script from the course text (Chapter 5.5)
  % This code will not work if A,b,c are not the right sizes; it assumes the 
  % origin is feasible for the problem (so b should have non-negative entries).
  % The code will return a optimization flag and iteration count.  The flag 
  % output will be returned as 0 if an optimal solution was found and -1 if the 
  % problem was detected as unbounded.
  [m,n] = size(b);
  if n>1
    b = b';
  end
  [m,n] = size(c);
  if m>1
    c = c';
  end
  iter = 0;
  opt = 0;
  while max(c)>eps
    % pick largest coefficient
    [cj,col] = max(c);
    Acol = A(:,col);
    % select leaving variable
    if sum(Acol<-eps) ==0
      opt = -1; %unbounded
      break;
    end
    nums = b.*(Acol<-eps);
    dens = -Acol.*(Acol<-eps);
    [t,row] = min(nums./dens);
    Arow = A(row,:);
    
    a = A(row,col); %pivot element
    A = A-Acol*Arow/a;
    A(row,:) = -Arow/a;
    A(:,col) = Acol/a;
    A(row,col) = 1/a;
    brow = b(row);
    b = b-brow*Acol/a;
    b(row) = -brow/a;
    
    ccol = c(col);
    c = c-ccol*Arow/a;
    c(col) = ccol/a;
    
    iter = iter+1;
  end
end
