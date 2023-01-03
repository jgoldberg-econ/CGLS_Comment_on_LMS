close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
sig=12; lamb = 1.21; rvec = [6:-0.1:0.1,0.05,0.02];
[pivec,~,~,~,~,lL,lF,l0] =compute_pi_fast(sig,lamb,n); 
pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
c=33.3569^2;
pi=pivec*c;
kap=3.9345;
pishr=pivec(n+1:-1:1)+pivec(n+1:end);

% The above code is identical to Lines 3-10 of LMS's script calibration_EMA_submit.m,
% except that compute_pi_fast.m here returns three additional scalar outputs:
% the production cost of a leader in a market with a productivity gap of 1 step (lL), 
% the production cost of a follower in such a market (lF), and 
% the production cost of a tied firm (l0).  See compute_pi_fast.m for details.  

%%%%%%%%%%%%%%%%%%%%%%%%
% Excercises with a discount rate decline from 3% to 1%
%%%%%%%%%%%%%%%%%%%%%%%%

% Solving for the BGP
r=3; dr=2;
xinit = zeros(1,2*n); xinit(n+1) = 1;
[xvec, muvec, ~, g,~] = gen_compute_eqm(lamb,pi,1,kap,r,xinit);
[xvec2, muvec2] = gen_compute_eqm(lamb,pi,1,kap,r-dr,xinit);
invcost=[xvec.^2/c,0]; invcost2=[xvec2.^2/c,0]; pishrvec=(pishr-invcost(n+1:end)-invcost(n+1:-1:1)); pishrvec2=(pishr-invcost2(n+1:end)-invcost2(n+1:-1:1));
LIst = pishrvec*muvec';
must = muvec*(0:n)';

% Lines 21-27 are identical to Lines 237-242 of LMS's main script 

LI_terminal = pishrvec2*muvec2';  % Profit share in the new BGP. Compare with Line 39 (profit share in initial BGP) of LMS's script.

eps=0.001; % As in LMS's code, a year is divided into periods of length epsilon years
horizon_in_years = 200+eps;  % Calculate dynamics for this number of years after the decline in the discount rate
T=horizon_in_years / eps;  % Number of periods over which to calculate transition dynamics

%%%%%%%%%%%% Obtaining transition dynamics

% When PDA_leader = 0, sc18_transition.m calculates transition dynamics
% under the assumption that the Spillover from Followers PDA holds.  
% See our Online Appendix A.6.  
% Under the Spillover from Followers PDA, a follower's productivity index 
% increases by 1 when the follower obtains an investment success or when
% technology diffusion occurs.  

% Obtain transition dynamics
PDA_leader = 0; % Use the Spillover from Follower PDA 
percent_adjust = 0.01; % When percent_adjust = 0.01, sc18_transition calculates the evolution of the market state distribution fixing the time-scale error in LMS's code
[~, transg_comp_fix100, transLI_fix100, ~,mu_fix100,ginit] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

% transg_comp_fix100 is the transition path for productivity growth, correcting the time-scale error and omission of composition effects in LMS's code
% The time-scale error and omission of composition terms in LMS's code are described in Section 1 of our comment

% Obtain transition dynamics, using LMS transition dynamics for the market state
% distribution (i.e., not fixing time-scale and composition-term errors)
PDA_leader = 0; 
percent_adjust = 1; % When percent_adjust = 1, sc18_transition calculates the evolution of the market state distribution as in LMS's code, without fixing the time-scale error
[transg, ~, transLI, ~,mu] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

% Next, obtain transition dynamics, fixing the time-scale and composition-term errors in LMS's code,
% and using the Spillover from *Leader* PDA.  See our Online Appendix A.5.
% Under the Spillover from Leader PDA, a leader's productivity index 
% increases by 1 when the leader obtains an investment success.  
% Results for the Spillover from Leader PDA are used only in Online Appendix Figure IA.2.
PDA_leader = 1;
percent_adjust = 0.01;
[~, transg_comp_fix100_PDAL, transLI_fix100_PDAL] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec2,percent_adjust,lL,lF,l0,sig,PDA_leader);

[status,msg,msgID] = mkdir('figures_comment'); % Figures will be saved in figures_comment folder; create folder if missing

%%% Figure 1 of our main text

T = T - 1;
horizon_in_years = 200;
pre_horizon_in_years = horizon_in_years / 4;
try close(fig1);catch;end
fig1=figure(1);

set(gcf, 'PaperUnits', 'inches');
    x_width=8.5;
    y_width=2.8;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %   

% Figure 1, Panel A, Productivity Growth

zoom = 100;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (0-Tst):(T/zoom);

subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg(1:T/zoom+1)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom+1)],'--r', 'LineWidth', 2); hold on;

    xlim([-2, xaxis(end) / (1/eps) * 4 ]); ylim([0.75 1.69]);
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title('A. Productivity growth', 'Interpreter', 'latex');

% Figure 1, Panel B, Gap

subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu(1:T/zoom+1)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu_fix100(1:T/zoom+1)],'--r', 'LineWidth', 2); hold on;

    xlim([-2, xaxis(end) / (1/eps) * 4 ]); ylim([1.5 7.5]); 
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Gap', 'Interpreter','latex');
    title('B. Average distance between leaders and followers', 'Interpreter', 'latex'); 
    l1 = legend({"LMS Figure 8"; newline+ "Correcting errors"}, ...
                  'Interpreter', 'latex', 'Location', 'northwest');
    set(l1, 'box', 'off');
    saveas(gcf, "figures_comment/corrected_LMS_Fig8.eps", 'epsc');  
   

%%% Figure 2 of our main text

pre_horizon_in_years = horizon_in_years / 4;
Tst = pre_horizon_in_years / eps;
try close(fig2);catch;end
fig2 = figure(2);

set(gcf, 'PaperUnits', 'inches');
    x_width=8.5;
    y_width=2.8;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %   

zoom = 1;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (0-Tst):(T/zoom);
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg(1:T/zoom+1)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,1); plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom+1)],'--r', 'LineWidth', 2); hold on;

    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlim([-200 800]);
    ylim([0.75 1.69]);
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Percent', 'Interpreter','latex');
    title('A. Productivity growth', 'Interpreter', 'latex');
    l1 = legend({"LMS Figure 8"; newline+ "Correcting errors"}, ...
                  'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off');

subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu(1:T/zoom+1)],'-k', 'LineWidth', 1.75); hold on;
subplot(1,2,2); plot(xaxis / (1/eps) * 4,[must*ones(1,Tst),mu_fix100(1:T/zoom+1)],'--r', 'LineWidth', 2); hold on;

    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlim([-200 800]);
    ylim([1.5 7.5]); 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Gap', 'Interpreter','latex');
    title('B. Average distance between leaders and followers', 'Interpreter', 'latex'); 
    hold off;
    saveas(gcf, "figures_comment/corrected_LMS_Fig8_800qtr.eps", 'epsc');  


%%% Figure IA.2

zoom = 1;
Tst = (pre_horizon_in_years / zoom) / eps;
xaxis = (0-Tst):(T/zoom);


try close(fig3);catch;end
fig3=figure(3); 
set(gcf, 'PaperUnits', 'inches');
    x_width=6.5;
    y_width=3.5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100(1:T/zoom+1)], ':b', 'LineWidth', 2); hold on;
plot(xaxis / (1/eps) * 4,[ginit*ones(1,Tst),transg_comp_fix100_PDAL(1:T/zoom+1)],'LineStyle', '-.', 'color', [0, 0.5, 0], 'LineWidth', 2); hold on;
xlim([xaxis(1) / (1/eps) * 4, xaxis(end) / (1/eps) * 4 ]); ylim([0.75 4]); legend('Spillover from Followers PDA','Spillover from Leaders PDA');
xlabel('Quarters since $r$ declined','FontSize',12,'fontweight','normal','interpreter','latex'); 
    ylim([0.7 3.5]); 
    ax = gca;
    ax.YGrid = 'on';
    box off; 
    xlabel('Quarters since $r$ declined','Interpreter','latex'); 
    ylabel('Productivity growth (Percent)', 'Interpreter','latex');
    l1 = legend(["Spillover from Followers PDA"; "Spillover from Leaders PDA"], ...
                 'Interpreter', 'latex', 'Location', 'northeast');
    set(l1, 'box', 'off')

    saveas(gcf, "figures_comment/PDA_figure.eps", 'epsc');
   
try close(fig3);catch;end  
