function Main_MEB

    %number of samples
    m=2^5;
    
    %number of variables
    n=2^12;

    % Number of runs for each instance:
    nrun = 1;

    % Number of solvers:
    nsolvers = 3;

    maxit = 3000;
    maxtime = 100;

    frun_cell = cell(nrun,nsolvers);
    timeVectot_cell = cell(nrun,nsolvers);

    fstop = zeros(nrun,1);

    for krun = 1:nrun

        %%%%% Generation of the instance: %%%%%

        % seed changes at every run to generate starting point 
        rng(krun);

        x0 = zeros(n,1);
        x0(1) = 1e0;

        Q = randn(m,n);
        c = sum(Q.^2,1)';
        Q = 2e0*(Q'*Q);

        %----------------------------------------------------------

        count = 1;
        eps = 1e-6;
        stopcr = 2;


        disp('****************');
        disp('* AWAY STEP FW *');
        disp('****************');


        [xfwaw,iterfwaw,fxfwaw,tottimefwaw,frun_cell{krun,count},timeVectot_cell{krun,count}]=...
            FWAW_Q(Q,c,x0,0,maxit,maxtime,eps,fstop(krun),stopcr);


        % Print results:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'0.5*xQX - cx = %10.3e\n',fxfwaw);
        fprintf(1,'Number of non-zero components of x = %d\n',...
            sum((abs(xfwaw)>=0.0001)));
        fprintf(1,'Number of iterations = %d\n',...
            iterfwaw);
        fprintf(1,'CPU time = %10.3e\n', tottimefwaw);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %----------------------------------------------------------

        count = count + 1;

        disp('*****************');
        disp('*  FW STANDARD  *');
        disp('*****************');


        [xfw,iterfw,fxfw,tottimefw,frun_cell{krun,count},timeVectot_cell{krun,count}]=...
            FW_Q(Q,c,x0,0,maxit,maxtime,eps,fstop(krun),stopcr);


        % Print results:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'0.5*xQX - cx = %10.3e\n',fxfw);
        fprintf(1,'Number of non-zero components of x = %d\n',...
            sum((abs(xfw)>=0.0001)));
        fprintf(1,'Number of iterations = %d\n',...
            iterfw);
        fprintf(1,'CPU time = %10.3e\n', tottimefw);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %----------------------------------------------------------

        count = count + 1;

        disp('*****************');
        disp('*      PG       *');
        disp('*****************');


        [x_pg,iter_pg,fx_pg,tottime_pg,frun_cell{krun,count},timeVectot_cell{krun,count}]=...
            PG_Q(Q,c,x0,0,maxit,maxtime,eps,fstop(krun),stopcr); 

        % Print results:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'0.5*xQX - cx = %10.3e\n',fx_pg);
        fprintf(1,'Number of non-zero components of x = %d\n',...
            sum((abs(x_pg)>=0.0001)));
        fprintf(1,'Number of iterations = %d\n',...
            iter_pg);
        fprintf(1,'CPU time = %10.3e\n', tottime_pg);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



    end

    fstop=min([fxfwaw fxfw fx_pg]);
    
    maxit = 0;
    for i = 1:nsolvers
        for krun = 1:nrun
            maxit = max(maxit,size(frun_cell{krun,i},2));
        end
    end
    frun = zeros(nrun,nsolvers,maxit);
    timerun = zeros(nrun,nsolvers,maxit);
    for i = 1:nsolvers
        for krun = 1:nrun
            niter = size(frun_cell{krun,i},2);
            fx = frun_cell{krun,i}(niter);
            t = timeVectot_cell{krun,i}(niter);
            frun(krun,i,1:niter) = frun_cell{krun,i};
            frun(krun,i,niter+1:maxit) = fx;
            timerun(krun,i,1:niter) = timeVectot_cell{krun,i};
            timerun(krun,i,niter+1:maxit) = t;
        end
    end

    clear frun_cell timeVectot_cell Q c
    save(strcat('Results_m_',num2str(m),'_n_',num2str(n)));

    %==========================================================================

    frun_plot = zeros(nsolvers,maxit);
    timerun_plot = zeros(nsolvers,maxit);
    for i = 1:nsolvers
        for krun = 1:nrun
            frun_plot(i,:) = frun_plot(i,:) + max(0e0,(squeeze(frun(krun,i,:))-fstop(krun)))';
            timerun_plot(i,:) = timerun_plot(i,:) + squeeze(timerun(krun,i,:))';
        end
    end
    frun_plot = frun_plot/nrun;
    timerun_plot = timerun_plot/nrun;

    %==========================================================================

    plot_threshold = 1e-20;

    % AFW
    plot_fwaw=max(frun_plot(1,1:maxit),plot_threshold);
    plot_tfwaw=timerun_plot(1,1:maxit);
   

    %==========================================================================

    % FW 
    plotfw=max(frun_plot(2,1:maxit),plot_threshold);
    plottfw=timerun_plot(2,1:maxit);
    %==========================================================================

    % PG 
    plot_pg=max(frun_plot(3,1:maxit),plot_threshold);
    plot_tpg=timerun_plot(3,1:maxit);
    %==========================================================================

   
    %==========================================================================

    % All Solvers
    figure
    semilogy(plot_tfwaw(1,:),plot_fwaw(1,:),'b-')
    hold on
    semilogy(plottfw(1,:),plotfw(1,:),'r-')
    semilogy(plot_tpg(1,:),plot_pg(1,:),'g-')
    title('FW variants & PG - objective function error')
    legend('AFW','FW','PG');
    xlim([0,30]);  ylim([10^(-10),10^5]);
    savefig(strcat('Plot_m_',num2str(m),'_n_',num2str(n),'_all_solvers'));
  %  close(gcf)

end