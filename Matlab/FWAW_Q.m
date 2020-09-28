function [x,it,fx,ttot,fh,timeVec] = FWAW_Q(Q,c,x,...
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
    
    %solution of Away step problem
    indcc=find(x>0e0);
    [~,istaraw]=max(g(indcc));

    xstar2=zeros(n,1);
    indaw=indcc(istaraw);
    xstar2(indaw)=1.0;

    %directions calculation
    dFW=xstar-x;
    dAS=x-xstar2;

    p1=g'*dFW;
    p2=g'*dAS;

    %keyboard

    %choice of the search direction
    if (p1<= p2)
        d=dFW;
        alpham=1.0;
        gnr = p1;
        cdir=1;
    else
        %display('AWAY-STEP chosen');
        d=dAS;
        alpham=x(indaw)/(1-x(indaw));
        gnr = p2;
        cdir=2;
    end
    
    % stopping criteria and test for termination
    switch stopcr
        case 1
            if (fx <= fstop)
                break;
            end
        case 2
            if (p1 >= -eps)
                break;
            end            
        otherwise
            error('Unknown stopping criterion');
    end
    
    %Armijo search
    alpha=alpham;
    ref = gamma*gnr;

    while (1)

        %Smart computation of the o.f. at the trial point
        if (cdir == 1)
            fz = 0.5*((1-alpha)^2*xQx +2*alpha*(1-alpha)*Qx(istar) + ...
                alpha^2*Q(istar,istar)) - ((1-alpha)*cx + alpha*c(istar));
        else
            fz = 0.5*((1+alpha)^2*xQx -2*alpha*(1+alpha)*Qx(indaw) + ...
                alpha^2*Q(indaw,indaw)) - ((1+alpha)*cx - alpha*c(indaw));
        end

        if (fz<=fx+alpha*ref)
            z=x+alpha * d;
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

%x = full(x);

end