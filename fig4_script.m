%Simulation of Fig 4 and Fig 5.
%The only difference between them is the TR.
%For Fig 4, use TR=20, for Fig 5, use TR=50

TR=20 %change this to TR=50 for Fig. 5

%values for the diffusion coefficient
D_vec=[0 0.0055e-3 0.0008 1.93e-3];

%spoil increments
psi_vec=[50 115.4 117 150 169];

%flip angles
flip_angles=[5:5:90];% acosd(exp(-TR/T1))]

%T1 and T2 values
T1TR_vec=[-1:0.1:2.3]; %logarithm of T1/TR
T2TR_vec=[-1:0.1:2.3]; %logarithm of T2/TR
T1_vec=[-1:0.1:2.3]+log10(TR); %logarithm of T1 for the given TR
T2_vec=[-1:0.1:2.3]+log10(TR); %logarithm of T2 for the given TR

%allocate memory (NaN - and not zeros(...) - is used for display purposes
signal_plus=NaN(length(D_vec),length(flip_angles),length(psi_vec),length(T1_vec),length(T2_vec)); %signal S+
ernst=NaN(length(flip_angles),length(T1_vec),length(T2_vec)); %ernst amplitudes

%signal computation
for ps=1:length(psi_vec);
    psi=psi_vec(ps);
    for d=1:length(D_vec)     
        tic
        [psi D_vec(d)]
        D=D_vec(d);
        for t1=1:length(T1_vec)
        T1=power(10,T1_vec(t1));
            for t2=1:length(T2_vec)
            T2=power(10,T2_vec(t2));
                if T1>=T2
                np=floor(T1)+50; %adapt pulse number to T1 to speed up computation
                    parfor k=1:length(flip_angles)
                    alpha=flip_angles(k);
                        if ps==1 &&  d==1 %calculate the ernst amlitude.
                           ernst(k,t1,t2)=sind(alpha)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(alpha));
                        end
                    [temp temp2]=epg_rfsp(alpha,np,T1,T2,TR,D,psi,1);
                    signal_plus(d,k,ps,t1,t2) = abs(temp(np));
                    end %eof flip_angles
                end %eof T1>=T2
            end %eof t2
        end %eof t1
        toc
    end %eof D_vec
end %eof psi_vec


%% compute and display the epsilon-values

%allocate memory
epsilon=NaN(length(D_vec),length(psi_vec),length(T2_vec)+1,length(T1_vec)+1); %the +1 is needed for display purposes

%which flip angles shall be considered for epsilon-value calculation (here we use all of them)
flip_indices=[1:length(flip_angles)];

psi_vec_count=[1 2 3 4 5]; %results for which psis should be displayed? (default: all of them)

for psic=1:length(psi_vec_count)
    psi=psi_vec_count(psic);
    for d=1:length(D_vec)
        for t2=1:length(T2_vec)
            for t1=1:length(T1_vec)
                if (T1_vec(t1)>=T2_vec(t2))
                    epsilon(d,psi,t2,t1)=0;
                    for p=1:length(flip_indices);
                        sig=signal_plus(d,flip_indices(p),psi,t1,t2);
                        diff=sig-ernst(flip_indices(p),t1,t2);
                        epsilon(d,psi,t2,t1)=epsilon(d,psi,t2,t1)+abs(diff)*(1/length(flip_indices)*1./ernst(flip_indices(p),t1,t2));
                    end
                end
            end
        end

    a=log10(squeeze(epsilon(d,psi,:,:))); %log of the epsilon-values - this is going to be displayed now
    figure;
    %display the epsilon-values for T1>=TR and T2>=TR (otherwise epsilon is approx. zero anyway)
    zero_index=find(T1TR_vec==0);
    h = pcolor([T1TR_vec(zero_index:end) 2.4],[T2TR_vec(zero_index:end) 2.4],a(zero_index:end,zero_index:end));set(h, 'EdgeColor', 'none');set(gca,'color',[0.65 0.65 0.65])
    colormap(jet); 
    colorbar;
    set(gca,'ColorScale','log');
    caxis([log10(0.005) log10(0.4)]);
    ylabel('log_{10}(T2/TR)','FontSize',14);
    xlabel('log_{10}(T1/TR)','FontSize',14);
  
    cbh = colorbar; 
    cbh.Ticks = [log10(0.005) log10(0.01) log10(0.05) log10(0.1) log10(0.15) log10(0.2) log10(0.25) log10(0.3) log10(0.35) log10(0.4)] ; % 0.1% 1% 2.5% 5% 7.5% 10% 15% 20% 
    set(cbh,'TickLabelInterpreter','none');
    cbh.TickLabels = {'0.5%' '1%' '5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%'};
   
    set(gca,'FontSize',14);
    
    fig_name = ['epsilon for psi=' num2str(psi_vec(psi)) '  D=' num2str(D_vec(d))];
    title(fig_name)
    end %eof psi_vec
end %eof D_vec
            