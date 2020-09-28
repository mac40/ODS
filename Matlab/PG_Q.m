function [x,it,fx,ttot,fh,timeVec] = PG_Q(Q,c,x,...
    verbosity,maxit,maxtime,eps,fstop,stopcr)
         
gamma=1e-4;

flagls=0;

%x = full(x);

[~,n] = size(Q);

if (maxit < Inf)
    fh=zeros(1,maxit);
    timeVec=zeros(1,maxit);
else
    fh=zeros(1,100*n);
    timeVec=zeros(1,100*n);
end

it=1;

proj_simplex_vector = @(y) max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);

tstart = tic;

while (it <= maxit && flagls==0)
    
    Qx = Q*x;
    xQx= x'*Qx;
    cx = c'*x;
    
    if (it==1)
        fx=0.5*xQx - cx;
        timeVec(it) = 0;  
    else
        %timeVec(it) = cputime - t;
        timeVec(it) = toc(tstart);
    end
    
    fh(it)=fx;
    
    % gradient evaluation
    g = Qx - c;
    
    if (timeVec(it) > maxtime)
        break;
    end
    
    % compute direction
    d = proj_simplex_vector(x-g)-x;
        
    % stopping criteria and test for termination
    switch stopcr
        case 1
            if (fx <= fstop)
                break;
            end
        case 2           
            [~,istar] = min(g);
            xstar = zeros(n,1);
            xstar(istar) = 1e0;
            if (g'*(xstar-x) >= -eps)
                break;
            end
        otherwise
            error('Unknown stopping criterion');
    end
        
    %Armijo search
    alpha=1e0;
    gd = g'*d;
    ref = gamma*gd;

    dQd = d'*Q*d;

    while (1)

        %fz = fx + alpha*(Qd'*x + c'*d) + 0.5*alpha*alpha*d'*Qd;
        %fz = fx + alpha*Qd'*(x+0.5*alpha*d) + alpha*c'*d;
        fz = fx + alpha*(gd + 0.5*alpha*dQd);

        if (fz <= fx + alpha*ref)
            z = x + alpha*d;
            break;
        else
            %disp(['alpha trial   = ' num2str(alpha) ' ' num2str(fz) ' ' num2str(fx)]);
            alpha=alpha*0.5;
        end

        if (alpha <= 1e-20)
            z=x;
            flagls=1;
            it = it-1;
            break;
        end

    end
    
    x = z;
    fx = fz;
          
    if (verbosity>0)
        disp(['-----------------** ' num2str(it) ' **------------------']);
        disp(['f(x)     = ' num2str(fx)]);
    end
    
    it = it+1;
        
end

%ttot = cputime - t;
ttot = toc(tstart);

if (it < size(fh,2))
    fh = fh(1:it);
    timeVec = timeVec(1:it);
end

end