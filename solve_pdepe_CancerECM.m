%% solve_pdepe_CancerECM

% % Ouput are tumor cells density; density of non-crosslink and cross-link
% % ECM; concentration of enzyme MMPs and LOXs 
% % in a metastatic tumor microenvironment.
% % The model describes the remodeling of ECM due to MMP and LOX and the 
% % migration of tumor cells through a remodeling ECM

%  The current model is inspired    
%  from the following primariy previous models of 
%  Anderson (2000), Gerisch (2008), and Andasari (2011)
%  This new model is further improved with extended features 
%  related to LOX and its effect on the whole system.

function sol = solve_pdepe_CancerECM(varargin)
% default values first

ploton = 1;



%% Author: Ye Nguyen 

%% Description of Input and output parameters or variables:
    %% Input
    % x             : Ind. spatial variable 
    % L             : reference length 
    
    % t             : Ind. simulation time variable
    % tau           : residence time 
    
    % c             : Dep. tumor cells density variable
    % co            : ref. value for c       
    
    % D             : reference chemical diffusion coefficient 
    % Dc            : diffusion coefficient of tumor cells 
    % Dc_hat        : dimensionless coefficient of Dc   
   
    % rho           : haptotaxis toward regular ECM 
    % rho_hat       : dimensionless rho
      
    % rho_cl_hat    : dimensionless haptotaxis toward cross-link ECM
    
    % gamma         : dimensionless proliferation of the primary tumor 
        
    % f             : dep. non-cross-link ECM density variable
    % f_cl          : cross-link ECM density 
    % fo            : reference parameter for f and f_cl
    
    % alpha_f_hat   : dimensionless MMPs uptake rate of ECM for degradation
    % beta_f_hat    : dimensionless LOXs uptake rate of ECM for crosslinkings 
    % mu_f_hat      : dimensionless remodeling rate of ECM
    
    % m             : MMPs concentration
    % mo            : ref value for m
    
    % Dm            : diffusion coefficient of enzyme MMPs
    % Dm_hat        : dimensionless Dm
    
    % alpha_m       : decay coefficient of MMPs
    % alpha_m_hat   : dimensionless alpha_m
    
    % beta_m        : secretion rate of MMPs
    % beta_m_hat    : dimensionless production rates of MMPs
       
    % l             : LOXs concentration 
    
    % Dl_hat        : dimensionless diffusion coefficient of enzyme LOX
    
    % alpha_l_hat   : dimensionless decay coefficient of LOXs
    
    % beta_l_hat    : dimensionless production rates of LOXs
    
    % epsilon       : A positive constant used in I.C.s 
    
%% Units of Input and output parameters or varibales:
    % L         : cm
    % x         : cm
    
    % t         : s 
    % t0        : s 
    
    % c         : cells/cm3
    
    % D         : cm2/s       
    % Dc        : cm2/s 
    
    % rho       : cm^2/(s*M)
    
    % gamma     : time^(-1)
    
    % f         : M 
    % f_cl      : M
    % fo        : M
    
    % m, mo     : M 
    
    % Dm        : cm2/s
    
    % l         : M 
    
    % Dl        : cm2/s
    
    % alpha_m   : s^(-1)
    % beta_m    : s^(-1)
        
%% Parameters value (either calculated or tentative)
        co          = 6.7*10^7    ; % cell/cm3  Anderson (2000)  
        fo          = 10^(-11)    ; % 10^-8 to 10^-11  from Anderson, 2000                
        mo          = 0.1*10^(-9) ; % Andasari(2011)
        
        L           = 1           ; % 0.1 to 1 cm 
        D           = 10^(-6)     ; % Bray 1992     
        tau         = 8*3600      ; % Anderson(2000)- 8 to 24 hours, take 8 hrs 
                              
        Dc          = 10^(-9)     ; %Bray 1992 
        Dc_hat      = Dc/D        ; %10^(-3) to 10^(-5) Chaplain (2006)
                
        rho         = 2600            ; % Anderson      
        rho_hat     = rho*fo*tau/(L^2);     
        rho_cl_hat  = 0.05            ; % Estimated dominates over non-cross-link ECM
        
        epsilon     = 0.001       ; % positive constant used in gamma 
        
        alpha_f_hat = 10          ; % Anderson (2000)
        mu_f_hat    = 0.15        ; % Andasari (2000)
        beta_f_hat  = 18          ; % Estimated

        Dm          = 10^(-9) ; % 10^-8 to 10^-10, Anderson(2000), Kumar(2015) 
        Dm_hat      = Dm/D        ;
  
        alpha_m     = 0.002       ; % sec^(-1) Kumar (2018)
        alpha_m_hat = 0.001       ; % estimated based on magnitude ratio between secretion and degration rate in Kumar (2018), i.e. 0.1/0.002~100      
        beta_m      = 0.1         ; % sec^(-1) Kumar (2018) varied rates from 0.005, 0.1, and 0.5 s^-1
        beta_m_hat  = 0.1         ; % estimated (as alpha) by Anderson (2000)
        
        Dl_hat      = 0.002       ; 
        
        alpha_l_hat = alpha_m_hat ; % assumed equals alpha_m_hat         
        beta_l_hat  = 0.1         ; 
        

            
%% PDEs solver "pdepe" set up and get called
n = 0; %rectangular coordinate for pdepe 

xmesh = linspace(0,1,500);
tspan = linspace(0,20,100);

% Read in chang in values if varargin is supplied to 
% solve_pdepe_CancerECM(parameters,tspan,ploton)
if nargin > 0
    parameters = varargin{1};
    
%% Parameters considered for localsensitivity_CancerECM
    % unpack parameters
    Dc_hat      = parameters(1);
    epsilon     = parameters(2);
    rho_hat     = parameters(3);
    rho_cl_hat  = parameters(4);
	alpha_f_hat = parameters(5); 
    mu_f_hat    = parameters(6);
    beta_f_hat  = parameters(7);
    Dm_hat      = parameters(8);
    alpha_m_hat = parameters(9);
    beta_m_hat  = parameters(10);
	Dl_hat      = parameters(11);
    alpha_l_hat = parameters(12);
    beta_l_hat  = parameters(13);
    tspan = varargin{2};
    ploton = varargin{3};
   
end

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
sol = pdepe(n,@mbpde,@mbic,@mbbc,xmesh,tspan,options);


%% Retrive values for each variables from output soln
c           = sol(:,:,1);
f           = sol(:,:,2);
f_cl        = sol(:,:,3);
m           = sol(:,:,4);
l           = sol(:,:,5);
spaceCheck  = 1-c-f-f_cl;

%desired output for global sensitivity may be 
% cMaxwrtXatFinalT = max(sol(tspan(end),:,1)); % max with respect to x, not
% max change in output wrt x
% fMaxwrtXatFinalT = ...
% ... 
% output = [cMaxwrtXatFinalT; fMaxwrtXatFinalT; ...]; % might need to be
% columns instead of rows

%% Plotting 
if ploton == 1
    %% Plot 1D numerical solution for the system at t = 0 
        %subplot_tight(2.5,2,1, [0.1 0.1])
        subplot(2, 2, 1)
        hold on 
        p1 = plot(xmesh,sol(1,:,1),'k',... 
             xmesh,sol(1,:,2),'m--',... 
             xmesh,sol(1,:,3),'b-.',... 
             xmesh,sol(1,:,4),'r.',...
             xmesh,sol(1,:,5),'gh') 
         
%         plot(xmesh, 1-sol(1,:,1)-sol(1,:,2)-sol(1,:,3),'c')
        
        p1(1).LineWidth = 2;
        p1(2).LineWidth = 2;
        p1(3).LineWidth = 2;
        p1(4).LineWidth = 2;
        ylim([0 1.1])
        xlabel('$x$','Interpreter','latex','FontSize',15)
        ylabel('$y$','Interpreter','latex','FontSize',15,'Rotation',0)
        title('$t = 0$','Interpreter','latex','FontSize',15)
        get(gca);set(gca,'FontSize',15,'FontName','Arial');
        
     %% Legend set up:
        ledg = legend('cancer cells',...
                'regular ECM fibers',...
                'cross-linked ECM fibers',...
                'MMP',...
                'LOX',... 
                'space check',...
                'position',[50 50 450 0]) %[left right across height]
%         hL    = subplot(3, 2, 5.5);
%         poshL = get(hL,'position'); % Getting its position
        
%         ledg = legend({'$c$',...
%                       '$f$',...
%                       '$f_{cl}$',...
%                       '$m$',...
%                       '$l$'},'Interpreter','latex')
        set(ledg,'location','east');
        ledg.FontSize = 12;   
        ledg.FontName ='Arial';
        
            
    %% Plot 1D numerical solution for the system at t = 1 
        %subplot_tight(2.5,2,2, [0.1 0.1])
        subplot(2, 2, 2)
        hold on
        p2 = plot(xmesh,sol(5,:,1),'k',... 
             xmesh,sol(5,:,2),'m--',... 
             xmesh,sol(5,:,3),'b-.',... 
             xmesh,sol(5,:,4),'r.',...
             xmesh,sol(5,:,5),'gd') 
        
%         plot(xmesh, 1-sol(5,:,1)-sol(5,:,2)-sol(5,:,3),'c')
         
        p2(1).LineWidth = 2;
        p2(2).LineWidth = 2;
        p2(3).LineWidth = 2;
        p2(4).LineWidth = 2;
        ylim([0 1.1])
        xlabel('$x$','Interpreter','latex','FontSize',15)
        ylabel('$y$','Interpreter','latex','FontSize',15,'Rotation',0)
        title('$t = 1$','Interpreter','latex','FontSize',15)
        get(gca);set(gca,'FontSize',15,'FontName','Arial');
                
    %% Plot 1D numerical solution for the system at t = 10 
        %subplot_tight(2.5,2,3, [0.1 0.1])
        subplot(2, 2, 3)
        hold on 
        p3 = plot(xmesh,sol(50,:,1),'k',... 
             xmesh,sol(50,:,2),'m--',... 
             xmesh,sol(50,:,3),'b-.',... 
             xmesh,sol(50,:,4),'r.',...
             xmesh,sol(50,:,5),'gd') 
        
%         plot(xmesh, 1-sol(50,:,1)-sol(50,:,2)-sol(50,:,3),'c')
        
        p3(1).LineWidth = 2;
        p3(2).LineWidth = 2;
        p3(3).LineWidth = 2;
        p3(4).LineWidth = 2;
        ylim([0 1.1])
        xlabel('$x$','Interpreter','latex','FontSize',15)
        ylabel('$y$','Interpreter','latex','FontSize',15,'Rotation',0)
        title('$t = 10$','Interpreter','latex','FontSize',15)
        get(gca);set(gca,'FontSize',15,'FontName','Arial');
      
    %% Plot 1D numerical solution for the system at t = 15
        %subplot_tight(2.5,2,4, [0.1 0.1])
        subplot(2, 2, 4)
        hold on
        p4 = plot(xmesh,sol(100,:,1),'k',... 
             xmesh,sol(100,:,2),'m--',...
             xmesh,sol(100,:,3),'b-.',... 
             xmesh,sol(100,:,4),'r.',...
             xmesh,sol(100,:,5),'gd') 
         
%         plot(xmesh, 1-sol(100,:,1)-sol(100,:,2)-sol(100,:,3),'c') 
         
        p4(1).LineWidth = 2;
        p4(2).LineWidth = 2;
        p4(3).LineWidth = 2;
        p4(4).LineWidth = 2;
        ylim([0 1.1])
        xlabel('$x$','Interpreter','latex','FontSize',15)
        ylabel('$y$','Interpreter','latex','FontSize',15,'Rotation',0)
        title('$t = 15$','Interpreter','latex','FontSize',15)
        get(gca);set(gca,'FontSize',15,'FontName','Arial');
        
        %% Export_fig - NOTE; Rename cases 
        set(gcf, 'color','w','Units','inches','Position', [0 0 10 7]);% [0 0 1000 600]);
        saveas(gcf, 'Plots/Case.png');
        export_fig ('Plots/Case','-m10','-painters','-png')
end

    %% Define system of PDEs
    function [a1, a2, a3] = mbpde(x_hat,t_hat,u,DuDx)
        
%         %% Unpack variable y into descriptive variables
        c   = u(1);
        f   = u(2);
        f_cl= u(3); 
        m   = u(4);
        l   = u(5);
        
        dcdx    = DuDx(1);
        dfdx    = DuDx(2);
        df_cldx = DuDx(3);
        dmdx    = DuDx(4);
        dldx    = DuDx(5);
                           
        %% Define PDEs system into the framework of "pdepe"
        
        %% Logistic growth/ physical space filling TERM
        gamma      = exp(-x_hat^2/epsilon);
        spaceCheck = (1-c-f-f_cl);
        growth     = gamma*c*spaceCheck;
        
        %% Cases:
        
            %% Case 1: Withoutout LOX effect
%             g = 0;
%             h = 0;
%             D_l = 0;
%             beta_l_hat = 0;
              

            %% Case 2: With LOX, NO haptotaxis effect toward cross-link ECM 
%             g = 0;
%             h = beta_f_hat*f*l;

            %% Case 3: With LOX AND haptotaxis effect toward cross-link ECM 
            g = -rho_cl_hat*spaceCheck*c*df_cldx;
            h = beta_f_hat*f*l;
                
        %% Model set up in "pdepe" form:            
            a1 = [1; 1; 1; 1; 1];
            
            a2 = [Dc_hat*dcdx-rho_hat*spaceCheck*c*dfdx+g;...                   
                  0;...
                  0;...              
                  Dm_hat*dmdx;... 
                  Dl_hat*dldx];
              
            a3 = [growth;...
                  -alpha_f_hat*m*f+mu_f_hat*spaceCheck-h;...
                  -alpha_f_hat*m*f_cl+h;...
                  beta_m_hat*c-alpha_m_hat*m;
                  beta_l_hat*c-alpha_l_hat*l];
           
    end
        
    %% Define the initial conditions at t = t0 
    function u0 = mbic(x)
        sigma = 0.01; % a positive constant, Anderson (2000)
        
        u0 = [exp(-x^2/sigma);...
              1-exp(-x^2/sigma);...
              0;...
              0;...
              0];
  
        %% Anderson (2000) conditions for case 1
%         u0 = [exp(-x^2/sigma);...
%               1-0.5*exp(-x^2/sigma);...
%               0;...
%               0.5*exp(-x^2/sigma);...
%               0];
    end 
             
    %% Define the boundary condions at x = a = 0 and x = b = 1
    function [pa, qa, pb, qb] = mbbc(xa,ua,xb,ub,t)
        %% Zero-flux in the left edge
        pa = [0; 0; 0; 0; 0];
        qa = [1; 1; 1; 1; 1];
        pb = [ub(1); ub(2)-1; ub(3); ub(4); ub(5)];
        qb = [0; 0; 0; 0; 0];
        
          %% Anderson(2000) conditions: Zero-flux in both edges
%         pa = [0; 0; 0; 0; 0];
%         qa = [1; 1; 1; 1; 1];
%         pb = [0; 0; 0; 0; 0];
%         qb = [1; 1; 1; 1; 1];
        
    end
end

        

