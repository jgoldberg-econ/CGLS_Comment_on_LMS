% This code generates Figure 3 of our main text and Figure IA.3 of our Online Appendix.

% First, obtain BGP productivity growth and BGP profit share for each
% discount rate in rvec, exactly as in LMS's code.  Lines 9-24 of bgp_figs.m match
% Lines 3-18 of LMS's main script calibration_EMA_submit.m exactly, without modification.
% (These results are used to generate the black solid line in our Figure 3.  The results
% are also used to generate the top row of our Figure IA.3.)

close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
sig=12; lamb = 1.21; rvec = [6:-0.1:0.1,0.05,0.02];
pivec=compute_pi_fast(sig,lamb,n);
pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
c=33.3569^2;
pi=pivec*c;
kap=3.9345;
pishr=pivec(n+1:-1:1)+pivec(n+1:end);

tic;for i=1:length(rvec)
    [xvec, muvec, ~, gvec(i),~,flag] = gen_compute_eqm(lamb,pi,1,kap,rvec(i),xinit);
    if flag>0;xinit=xvec;end
    invcost=[xvec.^2/c,0];
    avgmu(i)=(0:n)*muvec';
    LI(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
end;toc;

black = 'k';

% Second, obtain BGP productivity growth and BGP profit share for each
% discount rate in rvec, solving the model as described in LMS.  That is,
% gen_compute_eqm_correct solves the equations in Result 1 of our paper.
% See our Online Appendix A regarding Result 1.
%
% Result 1 is based on the assumption (on p. 213 of LMS) that the cost of
% achieving an investment success rate eta_s is (c eta_s)^2.
% These results will be used to generate the middle row of
% our Figure IA.3. 

c=33.3569; % per LMS's Table 1, p. 214 of their paper
xinit = zeros(1,2*n); xinit(n+1) = 1/100;

for i=1:length(rvec)
    [xvec_ad, muvec_ad, ~, gvec_ad(i),~,flag_ad] = gen_compute_eqm_correct(lamb,pivec,1,kap/100,rvec(i)/100,xinit,c);
    if flag_ad>0;xinit=xvec_ad;end
    invcost_ad=[c^2 * xvec_ad.^2,0];
    avgmu_ad(i)=(0:n)*muvec_ad';
    LI_ad(i)=(pishr-invcost_ad(n+1:end)-invcost_ad(n+1:-1:1))*muvec_ad';
end

% Third, obtain BGP productivity growth and BGP profit share for each
% discount rate in rvec, assuming that the cost of achieving an 
% investment success rate eta_s is (100/sqrt(2)/33.3569 eta_s)^2. 
% As described in our Online Appendix B, LMS's code solves for BGP investment success rates 
% as if the investment cost were (100/sqrt(2)/33.3569 eta_s)^2. 
% These results are used to obtain the red dashed line in our Figure 3.
% These results are also used to generate the bottom row of our Figure IA.3.
% Note: 
% Under the assumption that the cost of achieving an 
% investment success rate eta_s is (100/sqrt(2)/33.3569 eta_s)^2,
% BGP productivity growth is gvec_mod and BGP investment success rates are
% xvec_mod.  Observe that gvec_mod is equal to gvec / 100, where 
% gvec is the BGP productivity growth vector from LMS's code (Line 18 above).
% Thus, productivity growth under the investment cost assumption
% (100/sqrt(2)/33.3569 eta_s)^2 is indeed equal to productivity growth
% in LMS's replication code, once one accounts for gvec_mod being stored
% not in percent and gvec being stored (as in LMS's code) in percent.  Similarly,
% xvec_mod is equal to xvec / 100, for each interest rate in rvec.  

c = 100/sqrt(2)/33.3569;
xinit = zeros(1,2*n); xinit(n+1) = 1/100;

for i=1:length(rvec)
    [xvec_mod, muvec_mod, ~, gvec_mod(i),~,flag_mod] = gen_compute_eqm_correct(lamb,pivec,1,kap/100,rvec(i)/100,xinit,c);
    if flag_mod>0;xinit=xvec_mod;end
    invcost_mod=[c^2 * xvec_mod.^2,0];
    avgmu_mod(i)=(0:n)*muvec_mod';
    LI_mod(i)=(pishr-invcost_mod(n+1:end)-invcost_mod(n+1:-1:1))*muvec_mod';
end

gvec_outer = [gvec; 100*gvec_ad; 100*gvec_mod];
LI_outer = 1*[LI; LI_ad; LI_mod];

[status,msg,msgID] = mkdir('figures_comment');

% Our Figure 3

figure; 
set(gcf, 'PaperUnits', 'inches');
x_width=5;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

ax=  subplot(1,1,1);
p1 = plot(1*(rvec+gvec_outer(3, :)),LI_outer(3, :),'--','LineWidth',2,'Color','r'); 
hold on 
p2 = plot(1*(rvec+gvec_outer(1, :)), LI_outer(1, :), '-k', 'LineWidth', 1.75);
xtickformat(ax, 'percentage');
ylim([0.135 0.2]); 
ax.YGrid = 'on';
xlim([0, 7])
box off; 
set(gca,'ytick', [0.14 0.15 0.16 0.17, 0.18, 0.19, 0.20]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Profit share', 'Interpreter', 'latex')


l1 = legend([p1, p2], ...
           ["Correcting inconsistency in LMS code", "As reported in LMS Figure 4"], 'Interpreter', 'latex', ... 
           'Location', 'northeast', 'FontSize', 10);
l1.ItemTokenSize = [10; 5];
set(l1, 'box', 'off')
saveas(gcf, "figures_comment/inv_cost_profit_share.eps", 'epsc');


% Our Figure IA.3

figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=7.1;
y_width=2.15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

idx1 = find(abs(rvec - 3.6/1) == min(abs(rvec - 3.6/1)));
idx2 = find(abs(rvec - 0.3/1) == min(abs(rvec - 0.3/1)));

% Our Figure IA.3, Row 1, Left Panel (Productivity growth)
ax = subplot(1,2,1);
plot(1*(rvec+gvec_outer(1, :)),1*gvec_outer(1, :),'-','LineWidth',2,'Color',black); 
xtickformat(ax, 'percentage'); 
ytickformat(ax, 'percentage'); 
ylim([0 1.12]); 
ax.YGrid = 'on';
box off; 
xlim([0, 7]);
set(ax, 'ytick', [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Productivity growth rate', 'Interpreter', 'latex');
xline(1*(rvec(idx1)+gvec_outer(1, idx1)),'linestyle','--','linewidth',1.25);
xline(1*(rvec(idx2)+gvec_outer(1, idx2)),'linestyle','--','linewidth',1.25);
annotation('textarrow',[0.4 0.36],[0.5 0.6],'String','1984-2000');
annotation('textarrow',[0.24 0.19],[0.5 0.4],'String','2001-2016');

% Our Figure IA.3, Row 1, Right Panel (Profit share)
ax=  subplot(1,2,2);
plot(1*(rvec+gvec_outer(1, :)),LI_outer(1, :),'-','LineWidth',2,'Color',black); 
xtickformat(ax, 'percentage');
ylim([0.135 0.2]); 
ax.YGrid = 'on';
xlim([0, 7])
box off; 
set(gca,'ytick', [0.14 0.15 0.16 0.17, 0.18, 0.19, 0.20]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Profit share', 'Interpreter', 'latex')
xline(1*(rvec(idx1)+gvec_outer(1, idx1)),'linestyle','--','linewidth',1.25);
xline(1*(rvec(idx2)+gvec_outer(1, idx2)),'linestyle','--','linewidth',1.25);
sgtitle('\textbf{As reported in LMS Figure 4}', 'Interpreter', 'latex');

saveas(gcf, "figures_comment/IA_profit_share_row1.eps", 'epsc');

figure; 

% Our Figure IA.3, Row 2
set(gcf, 'PaperUnits', 'inches');
x_width=7.1;
y_width=2.15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


ax = subplot(1,2,1);
plot(1*(rvec+gvec_outer(2, :)),1*gvec_outer(2, :),'-','LineWidth',2,'Color',black); 
xtickformat(ax, 'percentage'); 
ytickformat(ax, 'percentage'); 
ylim([0 0.04]); 
ax.YGrid = 'on';
xlim([0, 7])
ax.YAxis.Exponent =0 ;
box off; 
title('Productivity growth rate', 'Interpreter', 'latex');
xlabel('Interest rate', 'Interpreter', 'latex');
xline(1*(rvec(idx1)+gvec_outer(2, idx1)),'linestyle','--','linewidth',1.25);
xline(1*(rvec(idx2)+gvec_outer(2, idx2)),'linestyle','--','linewidth',1.25);


ax=  subplot(1,2,2);
plot(1*(rvec+gvec_outer(2, :)),LI_outer(2, :),'-','LineWidth',2,'Color',black); 
xtickformat(ax, 'percentage');
ylim([0.135 0.2]); 
ax.YGrid = 'on';
box off; 
xlim([0, 7])
set(gca,'ytick', [0.14 0.15 0.16 0.17, 0.18, 0.19, 0.20]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Profit share', 'Interpreter', 'latex')
xline(1*(rvec(idx1)+gvec_outer(2, idx1)),'linestyle','--','linewidth',1.5);
xline(1*(rvec(idx2)+gvec_outer(2, idx2)),'linestyle','--','linewidth',1.5);
sgtitle({'\textbf{When correctly solving the model described in LMS}','\textbf{(i.e., investment cost $(c \eta)^2$ and $c=33.4$)}'}, 'Interpreter', 'latex');

saveas(gcf, "figures_comment/IA_profit_share_row2.eps", 'epsc');

% Our Figure IA.3, Row 3
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=7.1;
y_width=2.15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

ax = subplot(1,2,1);
plot(1*(rvec+gvec_outer(3,: )),1*gvec_outer(3, :),'-','LineWidth',2,'Color',black); 
xtickformat(ax, 'percentage'); 
ytickformat(ax, 'percentage'); 
ylim([0 1.12]); 
set(ax, 'ytick', [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1]); 
ax.YGrid = 'on';
ax.YAxis.Exponent =0 ;
xlim([0, 7])
title('Productivity growth rate', 'Interpreter', 'latex')
box off; 
xlabel('Interest rate', 'Interpreter', 'latex');
xline(1*(rvec(idx1)+gvec_outer(3, idx1)),'linestyle','--','linewidth',1.5);
xline(1*(rvec(idx2)+gvec_outer(3, idx2)),'linestyle','--','linewidth',1.5);

ax=  subplot(1,2,2);
p1 = plot(1*(rvec+gvec_outer(3, :)),LI_outer(3, :),'--','LineWidth',2,'Color','r'); 
hold on 
p2 = plot(1*(rvec+gvec_outer(1, :)), LI_outer(1, :), '-k', 'LineWidth', 1.75);
xtickformat(ax, 'percentage');
ylim([0.135 0.2]); 
ax.YGrid = 'on';
xlim([0, 7])
box off; 
set(gca,'ytick', [0.14 0.15 0.16 0.17, 0.18, 0.19, 0.20]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Profit share', 'Interpreter', 'latex')
xline(1*(rvec(idx1)+gvec_outer(3, idx1)),'linestyle','--','linewidth',1.5);
xline(1*(rvec(idx2)+gvec_outer( 3, idx2)),'linestyle','--','linewidth',1.5);

l1 = legend([p1, p2], ...
           {"Investment cost" + newline + "$\left(\frac{100}{c\sqrt{2}} \eta\right)^2$", "As reported in" + newline +  "LMS Figure 4"}, 'Interpreter', 'latex', ...
           'Position', [0.635, 0.55, 0.17, 0.0369], 'FontSize', 6);
l1.ItemTokenSize = [10; 5];
set(l1, 'box', 'off');
sgtitle({'\textbf{When solving the model with investment cost $\left(\frac{100}{c\sqrt{2}} \eta\right)^2$ and $c=33.4$}'}, 'Interpreter', 'latex');
saveas(gcf, "figures_comment/IA_profit_share_row3.eps", 'epsc');

