%Uses eFAST from Kirschner http://malthus.micro.med.umich.edu/lab/usadata/
%% PARAMETER INITIALIZATION
% set up max and mix matrices
scalevalue = 1e-1;
pmin=[VmaxoverKm*(scalevalue), % VmaxoverKm 
k_cat_Renin*(scalevalue), % k_cat_Renin
k_feedback*(scalevalue), % k_feedback
250, % feedback_capacity
k_cons_AngII*(scalevalue), % k_cons_AngII
1]; % dummy

scalevalue = 1e+1;
pmax=[VmaxoverKm*(scalevalue), % VmaxoverKm 
k_cat_Renin*(scalevalue), % k_cat_Renin
k_feedback*(scalevalue), % k_feedback
feedback_capacity*(scalevalue), % feedback_capacity
k_cons_AngII*(scalevalue), % k_cons_AngII
1]; % dummy
    % Parameter Labels 
efast_var={'V_{max}/K_M','k_R','k_f','k_{AII}','f','dummy'};

%             Si(i,t,u) = mean(Vi)/mean(V);
%             Sti(i,t,u) = 1-mean(Vci)/mean(V);
%             rangeSi(i,t,:,u) = Vi./V;
%             rangeSti(i,t,:,u) = 1-(Vci./V);
Model_efast
x=1:6; %5 input variables + dummy
% for i = x
%     for t = 1:length(Si(1,:,1))
%         for u = 1:length(Si(1,1,:))
%             meanSi(i,t,u) = mean(rangeSi(i,t,:,u));
%             stdSi(i,t,u) = std(rangeSi(i,t,:,u));
%             meanSti(i,t,u) = mean(rangeSti(i,t,:,u));
%             stdSti(i,t,u) = std(rangeSti(i,t,:,u));
%         end
%     end
% end
    
figure(1)
b=bar(x,[Si(1:3,1,1) Sti(1:3,1,1); Si(5,1,1) Sti(5,1,1); Si(4,1,1) Sti(4,1,1); Si(6,1,1) Sti(6,1,1)]);
b(1).FaceColor = [19 106 177]/255; %blue
b(2).FaceColor = [126 162 43]/255; %green
ax = gca;
ax.XTickLabel = {efast_var{1},efast_var{2},efast_var{3},efast_var{4},efast_var{5},efast_var{6},'FontName','Arial','FontSize',10};
% hold on
% errorbar(x,[meanSi(1:3,1,1) ; meanSi(5,1,1) ; meanSi(4,1,1); meanSi(6,1,1) ],...
%     [stdSi(1:3,1,1) ; stdSi(5,1,1) ; stdSi(4,1,1); stdSi(6,1,1) ],'.')
% hold off
legend('S_i','S_{Ti}','Location','NorthEast')
ylabel('eFAST Sensitivity of Ang II','FontName','Arial','FontSize',10)
set(gca,'FontName','Arial','FontSize',10);
axis([0 7 0 1])
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);
    
figure(2)
b=bar(x,[Si(1:3,1,2) Sti(1:3,1,2); Si(5,1,2) Sti(5,1,2); Si(4,1,2) Sti(4,1,2); Si(6,1,2) Sti(6,1,2)]);
b(1).FaceColor = [19 106 177]/255; 
b(2).FaceColor = [126 162 43]/255;
ax = gca;
ax.XTickLabel = {efast_var{1},efast_var{2},efast_var{3},efast_var{4},efast_var{5},efast_var{6},'FontName','Arial','FontSize',10};
legend('S_i','S_{Ti}','Location','NorthEast')
ylabel('eFAST Sensitivity of Ang I','FontName','Arial','FontSize',10)
set(gca,'FontName','Arial','FontSize',10);
axis([0 7 0 1])  
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);
 
figure(3)
b=bar(x,[Si(1:3,1,3) Sti(1:3,1,3); Si(5,1,3) Sti(5,1,3); Si(4,1,3) Sti(4,1,3); Si(6,1,3) Sti(6,1,3)]);
b(1).FaceColor = [19 106 177]/255;
b(2).FaceColor = [126 162 43]/255;
ax = gca;
ax.XTickLabel = {efast_var{1},efast_var{2},efast_var{3},efast_var{4},efast_var{5},efast_var{6},'FontName','Arial','FontSize',10};
legend('S_i','S_{Ti}','Location','NorthEast')
ylabel('eFAST Sensitivity of PRA','FontName','Arial','FontSize',10)
set(gca,'FontName','Arial','FontSize',10);
axis([0 7 0 1])
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 3.5 2.5]);
save(strcat('globalSens_',run_params.drugname,run_params.renalfunction,'.mat'),'-v7.3');    