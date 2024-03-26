% Simulation to produce Figure 1

%parameters
T1=400; %in ms
T2=400; %in ms
TR=20;  %in ms
D=0;    %no diffusion here
np=1000;%number of pulses

alpha_vec=[2.5:2.5:90 acosd(exp(-TR/T1))]; %vector of flip angles. Last entry is the Ernst-angle
psi_vec=[0:0.1:180]; %vector for phase increments psi

%allocate memory
signal_plus=zeros(length(alpha_vec),length(psi_vec));
signal_minus=zeros(length(alpha_vec),length(psi_vec));
ernst=zeros(1,length(alpha_vec));

%signal computation
for k=1:length(alpha_vec)
   alpha=alpha_vec(k)
   ernst(k)=sind(alpha)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(alpha));
   parfor ps=1:length(psi_vec)
       psi=psi_vec(ps);
       [sigp, sigm]=epg_rfsp(alpha,np,T1,T2,TR,D,psi,1);
       signal_plus(k,ps) = sigp(np);
       signal_minus(k,ps) = sigm(np);
   end
end

% plot the signals S+ nad S- in dependency of psi and for alpha=30°
pos30=find(alpha_vec==30);
figure;plot(psi_vec,squeeze(abs(signal_plus(pos30,:))));
hold on;plot(psi_vec,squeeze(abs(signal_minus(pos30,:))));
hold on;plot(psi_vec,squeeze(ernst(pos30)*ones(1,length(psi_vec))),'k');

ylabel('Signal [a.u.]','FontSize',14);
xlabel('\psi [deg]','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0:30:180],'TickDir','out');
axis([0 180 0 0.25]);
legend('S+','S-','Ernst ampl.','Location','northeast');

%plot the signal S+ for selected values of psi and in dependency of alpha
figure;plot([0 alpha_vec(1:end-1)],[0; abs(squeeze(signal_plus(1:end-1,find(psi_vec==50))))],'b','LineWidth',1);
hold on;plot([0 alpha_vec(1:end-1)],[0; abs(squeeze(signal_plus(1:end-1,find(10*psi_vec==1154))))],'r','LineWidth',1);
hold on;plot([0 alpha_vec(1:end-1)],[0; abs(squeeze(signal_plus(1:end-1,find(psi_vec==117))))],'g','LineWidth',1);
hold on;plot([0 alpha_vec(1:end-1)],[0; abs(squeeze(signal_plus(1:end-1,find(psi_vec==150))))],'cy','LineWidth',1);
hold on;plot([0 alpha_vec(1:end-1)],[0; abs(squeeze(signal_plus(1:end-1,find(psi_vec==169))))],'m','LineWidth',1);
hold on; plot([0 alpha_vec(1:end-1)],[0 ernst(1:end-1)],'k','LineWidth',1);

ylabel('Signal [a.u.]','FontSize',14);
xlabel('flip angle [deg]','FontSize',14);
set(gca,'XTick',[10:10:90],'TickDir','out');
axis([0 90 0 0.2]);
set(gca,'FontSize',14);
legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Ernst-curve','Location','northeast');
lgd = legend;
lgd.FontSize = 14;

% compute the epsilon-values (see manuscript)
flip_ind=[2:2:36];
epsilon=zeros(1,length(psi_vec));
a=abs(signal_plus);
b=repmat(ernst',1,length(psi_vec));
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./b(flip_ind,:); %ernst(end) is the ernst amplitude at the ernst angle
epsilon=100*sum(diffs,1)/18;

fprintf('\nepsilon-values in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));
