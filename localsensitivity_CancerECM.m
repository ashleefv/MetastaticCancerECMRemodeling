%% Create Ycalc matrix for model output at Xdata values for original set of parameters

%% Project: CancerECM mathematical models
clear all

%% Parameters considered and their base case values:
    Dc_hat      = 0.001  ;
    epsilon     = 0.001  ;
    rho_hat     = 0.0075 ;
    rho_cl_hat  = 0.05   ;
    alpha_f_hat = 10     ;
    mu_f_hat    = 0.15   ;
    beta_f_hat  = 18     ;
    Dm_hat      = 0.001  ;
    alpha_m_hat = 0.001  ;
    beta_m_hat  = 0.1    ;
    Dl_hat      = 0.002  ;
    alpha_l_hat = 0.001  ; 
    beta_l_hat  = 0.1    ;
    
vector_of_default_parameters = [Dc_hat, epsilon, rho_hat, rho_cl_hat,...
                                alpha_f_hat, mu_f_hat, beta_f_hat,...
                                Dm_hat, alpha_m_hat, beta_m_hat,...
                                Dl_hat, alpha_l_hat, beta_l_hat];
%% Simulation set up:
xmesh = linspace(0,1,500);
tspan = linspace(0,20,100);

num_cases = 2;
number_of_scenarios = num_cases*length(vector_of_default_parameters); % each of the 13 default params has two cases, up & down by sensitive_range*100% 
sensitivity_range = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:number_of_scenarios+1
    parameters(i,:) = vector_of_default_parameters;
end

%% Create updated parameters
for j = 1:length(vector_of_default_parameters)
    parameters(j*num_cases,j)   = vector_of_default_parameters(j)*(1+sensitivity_range);
    parameters(j*num_cases+1,j) = vector_of_default_parameters(j)*(1-sensitivity_range); %delete this row to remove moving down by sentivity_range
end

Ycalc_baseline = solve_pdepe_CancerECM(parameters(1,:),tspan,0);
[a,b,~] = size(Ycalc_baseline);
% baseline has no change
normalized_change_Y1 = zeros(a,b,2*length(vector_of_default_parameters)+1);
normalized_change_Y2 = zeros(a,b,2*length(vector_of_default_parameters)+1);
normalized_change_Y3 = zeros(a,b,2*length(vector_of_default_parameters)+1);
normalized_change_Y4 = zeros(a,b,2*length(vector_of_default_parameters)+1);
normalized_change_Y5 = zeros(a,b,2*length(vector_of_default_parameters)+1);

for j = 2: number_of_scenarios+1
%% Create Ycalc matrix for model output at time values
    Ycalc = solve_pdepe_CancerECM(parameters(j,:),tspan,0);
    % for this example Ycalc has five output vs. distance: Y1 vs. x, Y2,
    % Y3, Y4 and Y5 corresponding to c, f, f_cl, m and l in CancerECM model 
    Y1(:,:,j) = Ycalc(:,:,1);
    Y2(:,:,j) = Ycalc(:,:,2);
    Y3(:,:,j) = Ycalc(:,:,3);
    Y4(:,:,j) = Ycalc(:,:,4);
    Y5(:,:,j) = Ycalc(:,:,5);
    change_Y1(:,:,j) = abs(Ycalc(:,:,1)-Ycalc_baseline(:,:,1));
    change_Y2(:,:,j) = abs(Ycalc(:,:,2)-Ycalc_baseline(:,:,2));
    change_Y3(:,:,j) = abs(Ycalc(:,:,3)-Ycalc_baseline(:,:,3));
    change_Y4(:,:,j) = abs(Ycalc(:,:,4)-Ycalc_baseline(:,:,4));
    change_Y5(:,:,j) = abs(Ycalc(:,:,5)-Ycalc_baseline(:,:,5));
    
end

for j = 1:length(vector_of_default_parameters)
    normalized_change_Y1(:,:,2*j)   = change_Y1(:,:,2*j)./Ycalc_baseline(:,:,1)/sensitivity_range;
    normalized_change_Y1(:,:,2*j+1) = change_Y1(:,:,2*j+1)./Ycalc_baseline(:,:,1)/sensitivity_range;
    normalized_change_Y2(:,:,2*j)   = change_Y2(:,:,2*j)./Ycalc_baseline(:,:,2)/sensitivity_range;
    normalized_change_Y2(:,:,2*j+1) = change_Y2(:,:,2*j+1)./Ycalc_baseline(:,:,2)/sensitivity_range;
    normalized_change_Y3(:,:,2*j)   = change_Y3(:,:,2*j)./Ycalc_baseline(:,:,3)/sensitivity_range;
    normalized_change_Y3(:,:,2*j+1) = change_Y3(:,:,2*j+1)./Ycalc_baseline(:,:,3)/sensitivity_range;
    normalized_change_Y4(:,:,2*j)   = change_Y4(:,:,2*j)./Ycalc_baseline(:,:,4)/sensitivity_range;
    normalized_change_Y4(:,:,2*j+1) = change_Y4(:,:,2*j+1)./Ycalc_baseline(:,:,4)/sensitivity_range;
    normalized_change_Y5(:,:,2*j)   = change_Y5(:,:,2*j)./Ycalc_baseline(:,:,5)/sensitivity_range;
    normalized_change_Y5(:,:,2*j+1) = change_Y5(:,:,2*j+1)./Ycalc_baseline(:,:,5)/sensitivity_range;
end

%% specify the component Y1-Y5  and time of interest 
    %     comp_no = 1;
    %     comp    = Y1;    
    %     normalized_change_comp = normalized_change_Y1;
    %     change_comp = change_Y1;
    %     time = 100;

%% Figure 1: normalized change vs. space at a time of interest.
%% Ycalc(s) and Ycalc_baseline vs. space at a time of interest 
%     figure(1)
%             subplot(1, 2, 1)       
%             hold on
%             for j = 2:2:number_of_scenarios+1
%                 p1 = plot(xmesh,normalized_change_comp(time,:,j),'linewidth',2)
%             end
%             for j = 3:2:number_of_scenarios+1
%                 p1 = plot(xmesh,normalized_change_comp(time,:,j),'--','linewidth',2)
%             end       
%             hold off
% 
%             subplot(1, 2, 2)
%             hold on
%             for j = 2:2:number_of_scenarios+1
%                 p1 = plot(xmesh,comp(time,:,j),'linewidth',2)
%             end
%             for j = 3:2:number_of_scenarios+1
%                 p1 = plot(xmesh,comp(time,:,j),'--','linewidth',2)
%             end
%             plot(xmesh, Ycalc_baseline(time,:,comp_no),':k','linewidth',3)        
%             hold off
%              
%         leg1 = legend('$\hat{D}_c$','$\epsilon$','$\hat{\rho}$','$\hat{\rho}_{cl}$',...
%                       '$\hat{\alpha}_f$','$\hat{\mu}_f$','$\hat{\beta}_f$',...
%                       '$\hat{D}_m$','$\hat{\alpha}_m$','$\hat{\beta}_m$',...
%                       '$\hat{D}_l$','$\hat{\alpha}_l$','$\hat{\beta}_l$',...
%                       ...
%                       '$\hat{D}_c$','$\epsilon$','$\hat{\rho}$','$\hat{\rho}_{cl}$',...
%                       '$\hat{\alpha}_f$','$\hat{\mu}_f$','$\hat{\beta}_f$',...
%                       '$\hat{D}_m$','$\hat{\alpha}_m$','$\hat{\beta}_m$',...
%                       '$\hat{D}_l$','$\hat{\alpha}_l$','$\hat{\beta}_l$'); % customize this
%         set(leg1,'Interpreter','latex'); 
        
%         %% Export_fig - NOTE; Rename cases 
%         set(gcf, 'Position', [0 0 1000 600]);
%         saveas(gcf, 'Local Sensitivity/Y3.png');
%         export_fig Y3.png  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 2: Plot max/min of normalized_change_Y increasing/decreasing change        
% increasingchange=normalized_change_comp(time,:,2:2:number_of_scenarios+1);
% decreasingchange=normalized_change_comp(time,:,3:2:number_of_scenarios+1);
% [gsort,gindex] = sort(max(increasingchange,[],2))
% [hsort,hindex] = sort(max(decreasingchange,[],2))

%% Figure 2: Plot max/min of normalized_change_Y increasing/decreasing change
%     figure(2)
%         subplot(221)
%         b=bar(1:13,reshape(max(increasingchange,[],2),13,1))
%         title('increasing')
%         subplot(222)
%         b=bar(1:13,reshape(min(increasingchange,[],2),13,1))
%         title('increasing')
%         subplot(223)
%         b=bar(1:13,reshape(max(decreasingchange,[],2),13,1))
%         title('decreasing')
%         subplot(224)
%         b=bar(1:13,reshape(min(decreasingchange,[],2),13,1))
%         title('decreasing')    
        
    %     %% Export_fig - NOTE; Rename cases 
    %     set(gcf, 'Position', [0 0 1000 600]);
    %     saveas(gcf, 'Local Sensitivity/Yb_3.png');
    %     export_fig Yb_3.png  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 3: max normalized change w.r.t position vs. time    

%     for i = 1:length(tspan)  
%     increasingchange=normalized_change_comp(i,:,2:2:number_of_scenarios+1);
%     decreasingchange=normalized_change_comp(i,:,3:2:number_of_scenarios+1);
%     increasingVStime(:,i)=reshape(max(increasingchange,[],2),11,1);
%     decreasingVStime(:,i)=reshape(max(decreasingchange,[],2),11,1);
%     end     
    
%     figure(3)
%     
%     subplot(211)
%     hold on
%     for i = 1:7
%         semilogy(tspan,increasingVStime(i,:),'linewidth',1)
%     end
%     for i = 8:11
%         semilogy(tspan,increasingVStime(i,:),':','linewidth',2)
%     end 
%     grid on 
%     hold off
%      
%     xlabel('$t$','Interpreter','latex','FontSize',20)
%     ylabel('Absolute Value of Sensitivity','Interpreter','latex','FontSize',20)
%     title('Increasing Parameters','Interpreter','latex','FontSize',20)
%     
%     %% set up legend
%     leg1 = legend('$\hat{D}_c$','$\hat{\epsilon}$','$\hat{\rho}$','$\hat{\rho}_{cl}$',...
%               '$\hat{\alpha}_f$','$\hat{\mu}_f$','$\hat{\beta}_f$',...
%               '$\hat{D}_m$','$\hat{\beta}_m$',...
%               '$\hat{D}_l$','$\hat{\beta}_l$',...
%                'Orientation','horizontal'); % customize this
%            
%     set(leg1,'Interpreter','latex'); %'location','northeastoutside');
%     leg1.FontSize = 20;   
%     
%     subplot(212)
%     hold on
%     for i = 1:7
%         semilogy(tspan,decreasingVStime(i,:),'linewidth',1)
%     end 
%     for i = 8:11
%         semilogy(tspan,decreasingVStime(i,:),':','linewidth',2)
%     end
%     grid on
%     hold off    
%     
%     xlabel('$t$','Interpreter','latex','FontSize',20)
%     ylabel('Absolute Value of Sensitivity','Interpreter','latex','FontSize',20)
%     title('Decreasing Parameters','Interpreter','latex','FontSize',20)
%                
%     %% Export_fig - NOTE; Rename cases 
%         set(gcf, 'Position', [0 0 1000 600]);
%         saveas(gcf, 'Local Sensitivity/Y5.png');
%         export_fig Y5.png 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 4: combined Figure 3 and Figure 1 comp. Y5

% for i = 1:length(tspan)  
%     increasingchange=change_comp(i,:,2:2:number_of_scenarios+1);
%     decreasingchange=change_comp(i,:,3:2:number_of_scenarios+1);
%     [maxinc,index_maxinc]=max(increasingchange,[],2);
%     index_maxinc = reshape(index_maxinc,11,1);
%     increasingVStime(:,i)=reshape(maxinc,11,1)./sensitivity_range;%./reshape(Ycalc_baseline(i,index_maxinc,comp_no),11,1);
%     [maxdec,index_maxdec] = max(decreasingchange,[],2);
%     index_maxdec = reshape(index_maxdec,11,1);
%     decreasingVStime(:,i)=reshape(maxdec,11,1)./sensitivity_range;%./reshape(Ycalc_baseline(i,index_maxdec,comp_no),11,1);
% end

    figure(4)
    
    %% Y1-cancer cells
    subplot(4, 2, 1)
    
    for i = 1:length(tspan)  
    increasingchange=change_Y1(i,:,2:2:number_of_scenarios+1);
    decreasingchange=change_Y1(i,:,3:2:number_of_scenarios+1);
    
    [maxinc,index_maxinc]=max(increasingchange,[],2);
    index_maxinc = reshape(index_maxinc,13,1);
    increasingVStime(:,i)=reshape(maxinc,13,1)./sensitivity_range;
    
    [maxdec,index_maxdec] = max(decreasingchange,[],2);
    index_maxdec = reshape(index_maxdec,13,1);
    decreasingVStime(:,i)=reshape(maxdec,13,1)./sensitivity_range;
    end
    
    hold on
    for i = 1:7
        plot(tspan,increasingVStime(i,:),'linewidth',1)
    end
    for i = 8:13
        plot(tspan,increasingVStime(i,:),':','linewidth',2)
    end 
    plot(tspan, ones(size(tspan)),'k--','linewidth',1)
    ylim([0 3.5])
     
    xlabel('$t$','Interpreter','latex','FontSize',15)
    ylabel('Sensitivity','Interpreter','latex','FontSize',15)
    title('(A) cancer cells','Interpreter','latex','FontSize',15)
      
    %% Y2-regular ECM
    subplot(4, 2, 2)
    
    for i = 1:length(tspan)  
    increasingchange=change_Y2(i,:,2:2:number_of_scenarios+1);
    decreasingchange=change_Y2(i,:,3:2:number_of_scenarios+1);
    
    [maxinc,index_maxinc]=max(increasingchange,[],2);
    index_maxinc = reshape(index_maxinc,13,1);
    increasingVStime(:,i)=reshape(maxinc,13,1)./sensitivity_range;
    
    [maxdec,index_maxdec] = max(decreasingchange,[],2);
    index_maxdec = reshape(index_maxdec,13,1);
    decreasingVStime(:,i)=reshape(maxdec,13,1)./sensitivity_range;
    end
    
    hold on
    for i = 1:7
        plot(tspan,increasingVStime(i,:),'linewidth',1)
    end
    for i = 8:13
        plot(tspan,increasingVStime(i,:),':','linewidth',2)
    end 
    plot(tspan, ones(size(tspan)),'k--','linewidth',1)
    ylim([0 1.5])
    
     
    xlabel('$t$','Interpreter','latex','FontSize',15)
    ylabel('Sensitivity','Interpreter','latex','FontSize',15)
    title('(B) regular ECM fibers','Interpreter','latex','FontSize',15)
    
    %% Y3-cross-linked ECM
    subplot(4, 2, 3)
    
    for i = 1:length(tspan)  
    increasingchange=change_Y3(i,:,2:2:number_of_scenarios+1);
    decreasingchange=change_Y3(i,:,3:2:number_of_scenarios+1);
    
    [maxinc,index_maxinc]=max(increasingchange,[],2);
    index_maxinc = reshape(index_maxinc,13,1);
    increasingVStime(:,i)=reshape(maxinc,13,1)./sensitivity_range;
    
    [maxdec,index_maxdec] = max(decreasingchange,[],2);
    index_maxdec = reshape(index_maxdec,13,1);
    decreasingVStime(:,i)=reshape(maxdec,13,1)./sensitivity_range;
    end
    
    hold on
    for i = 1:7
        plot(tspan,increasingVStime(i,:),'linewidth',1)
    end
    for i = 8:13
        plot(tspan,increasingVStime(i,:),':','linewidth',2)
    end 
    plot(tspan, ones(size(tspan)),'k--','linewidth',1)
    ylim([0 1.5])
     
    xlabel('$t$','Interpreter','latex','FontSize',15)
    ylabel('Sensitivity','Interpreter','latex','FontSize',15)
    title('(C) cross-linked ECM fibers','Interpreter','latex','FontSize',15)
    
%     set(gca,'Position',[0.1 .1 0.5 0.1]) %[left bottom width height]
    
    %% Y4-MMP
    subplot(4, 2, 4)
    
    for i = 1:length(tspan)  
    increasingchange=change_Y4(i,:,2:2:number_of_scenarios+1);
    decreasingchange=change_Y4(i,:,3:2:number_of_scenarios+1);
    
    [maxinc,index_maxinc]=max(increasingchange,[],2);
    index_maxinc = reshape(index_maxinc,13,1);
    increasingVStime(:,i)=reshape(maxinc,13,1)./sensitivity_range;
    
    [maxdec,index_maxdec] = max(decreasingchange,[],2);
    index_maxdec = reshape(index_maxdec,13,1);
    decreasingVStime(:,i)=reshape(maxdec,13,1)./sensitivity_range;
    end
    
    hold on
    for i = 1:7
        plot(tspan,increasingVStime(i,:),'linewidth',1)
    end
    for i = 8:13
        plot(tspan,increasingVStime(i,:),':','linewidth',2)
    end 
    plot(tspan, ones(size(tspan)),'k--','linewidth',1)
    ylim([0 1.5])
     
    xlabel('$t$','Interpreter','latex','FontSize',15)
    ylabel('Sensitivity','Interpreter','latex','FontSize',15)
    title('(D) MMP','Interpreter','latex','FontSize',15)
    
%     set(gca,'Position',[0.1 .1 0.5 0.1]) %[left bottom width height]
    
    %% Y5-LOX
    subplot(4, 2, 5)
    
    for i = 1:length(tspan)  
    increasingchange=change_Y5(i,:,2:2:number_of_scenarios+1);
    decreasingchange=change_Y5(i,:,3:2:number_of_scenarios+1);
    
    [maxinc,index_maxinc]=max(increasingchange,[],2);
    index_maxinc = reshape(index_maxinc,13,1);
    increasingVStime(:,i)=reshape(maxinc,13,1)./sensitivity_range;
    
    [maxdec,index_maxdec] = max(decreasingchange,[],2);
    index_maxdec = reshape(index_maxdec,13,1);
    decreasingVStime(:,i)=reshape(maxdec,13,1)./sensitivity_range;
    end
    
    hold on
    for i = 1:7
        plot(tspan,increasingVStime(i,:),'linewidth',1)
    end
    for i = 8:13
        plot(tspan,increasingVStime(i,:),':','linewidth',2)
    end 
    plot(tspan, ones(size(tspan)),'k--','linewidth',1)
    ylim([0 1.5])
     
    xlabel('$t$','Interpreter','latex','FontSize',15)
    ylabel('Sensitivity','Interpreter','latex','FontSize',15)
    title('(E) LOX','Interpreter','latex','FontSize',15)
    
    %% figure1 of Y5-LOX
    p1 = subplot(4, 2, 6)       
    time = 100;
    hold on
    for j = 2:2:14
        plot(xmesh,Y5(time,:,j),'linewidth',1)
    end  
    for j = 16:2:number_of_scenarios+1
        plot(xmesh,Y5(time,:,j),':','linewidth',2)
    end  
    plot(xmesh, Ycalc_baseline(time,:,5),'k--','linewidth',1)
    ylim([0 1.5])
    
    xlabel('$x$','Interpreter','latex','FontSize',15)
    ylabel('Concentration','Interpreter','latex','FontSize',15)
    title('(F) LOX at t = 20','Interpreter','latex','FontSize',15)
    
%     leg1 = legend('$\hat{D}_c$','$\hat{\epsilon}$','$\hat{\rho}$','$\hat{\rho}_{cl}$',...
%           '$\hat{\alpha}_f$','$\hat{\mu}_f$','$\hat{\beta}_f$',...
%           '$\hat{D}_m$','$\hat{\alpha}_m$','$\hat{\beta}_m$',...
%           '$\hat{D}_l$','$\hat{\alpha}_m$','$\hat{\beta}_l$','base line')
%       leg1.FontSize = 18
%       set(leg1,'Interpreter','latex');
%       set(leg1,'location','northeast','Orientation','horizontal','NumColumns',3);
    
    %% set up legend
    hL    = subplot(4, 2, 7.5);
    poshL = get(hL,'position');  
    leg1 = legend([p1],...
            {'$\hat{D}_c$','$\hat{\epsilon}$','$\hat{\rho}$','$\hat{\rho}_{cl}$',...
          '$\hat{\alpha}_f$','$\hat{\mu}_f$','$\hat{\beta}_f$',...
          '$\hat{D}_m$','$\hat{\alpha}_m$','$\hat{\beta}_m$',...
          '$\hat{D}_l$','$\hat{\alpha}_l$','$\hat{\beta}_l$','base line'},...
          'Interpreter','latex',...
          'Orientation','horizontal','NumColumns',7); % customize this        
      
    set(leg1,'Interpreter','latex'); %'location','northeastoutside');
    leg1.FontSize = 12
    set(leg1,'position',poshL);      % Adjusting legend's position
	axis(hL,'off');   

    %% Export_fig - NOTE;  
    set(gcf, 'color','w','Units','inches','Position', [0 0 10 11]);
    saveas(gcf, 'LocalSensitivity/FigureCombinedMATLAB.png')
    export_fig('LocalSensitivity/LocalSensitivity','-m10','-painters','-png') 
    
%     modifiedpwd = pwd;
%     filename = strcat(modifiedpwd,sprintf('/LocalSensitivity/FigureCombined'));
%     if exist(strcat(filename,'.png'),'file')
%         delete(strcat(filename,'.png'))
%     end
%     
% %     saveas(gcf, 'Local Sensitivity/FigureCombinedMATLAB.png');
%     export_fig(filename,'-r1000','-a4','-q101','-painters','-png') 


    
    
    
