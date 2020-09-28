function [x,it,fx,ttot,fh,timeVec] = FW_Q(Q,c,x,...
    verbosity,maxit,maxtime,eps,fstop,stopcr)

gamma=1e-4;

flagls=0;

[~,n] = size(Q);

if (maxit < Inf)
    fh=zeros(1,maxit);
    timeVec=zeros(1,maxit);
else
    fh=zeros(1,100*n);
    timeVec=zeros(1,100*n);
end

it=1;

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
    
    %solution of FW problem
    [~,istar]=min(g);
    xstar=zeros(n,1);
    xstar(istar)=1.0;

    %direction calculation
    d = xstar - x;
    gnr = g'*d;
        
    % stopping criteria and test for termination
    switch stopcr
        case 1
            if (fx <= fstop)
                break;
            end
        case 2
            if (gnr >= -eps)
                break;
            end
        otherwise
            error('Unknown stopping criterion');
    end
        
    %Armijo search
    alpha=1;
    ref = gamma*gnr;

    while (1)
        %z=x+alpha*d;
        %Smart computation of the o.f. at the trial point
        fz = 0.5*((1-alpha)^2*xQx +2*alpha*(1-alpha)*Qx(istar) + ...
            alpha^2*Q(istar,istar)) - ((1-alpha)*cx + alpha*c(istar));

        if (fz<=fx+alpha*ref)
            z = x + alpha*d;
            break;
        else
            %disp(['alpha trial   = ' num2str(alpha) ' ' num2str(fz) ' ' num2str(fx)]);
            alpha=alpha*0.5;
        end

        if (alpha <= 1e-20)
            z=x;
            fz=fx;
            flagls=1;
            it = it-1;
            break;
        end

    end

    x=z;
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

%x = full(x);

end