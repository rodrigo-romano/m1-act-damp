%
%
%

%% Load optical sensitivity matrices
%%

% Radians to mili arc second conversion constant
rad2mas = (180/pi * 3600 * 1000);

% Root of the mean squared value function
rms = @(x,dir)squeeze(sqrt(mean(x.^2,dir)));

% WFE conversion cte
mas2nm = 10.2;  % [nm/mas]

opticalT_dtfolder = '/home/rromano/Workspace/gmt-data/optical_T';

load(fullfile(opticalT_dtfolder,'lom_tt_dt.mat'),'D_seg_tt');
load(fullfile(opticalT_dtfolder,'D_seg_piston_dt.mat'),'D_seg_piston');
Dttp = [rad2mas*D_seg_tt; D_seg_piston];


%%

for i1 = 1:3
    % Load simulation data
    load(sprintf('m1act_damping_%d_cfdwl.mat', i1-1),...
        'm1rbm_dt','m2rbm_dt','m1HPlc_dt');
    N = size(m1rbm_dt,1);
    t = ((1:N)-1)* 1e-3;

    ptt_data = Dttp * [m1rbm_dt, m2rbm_dt]';
    tt_data = rad2mas*ptt_data(1:14,:);

    wfeTT = sqrt((1/7)* sum(tt_data(1:7,:).^2 + tt_data(8:14,:).^2 ,1));
    wfeP = 1e9*sqrt((1/7)* sum((ptt_data(15:21,:) -...
        mean(ptt_data(15:21,:))).^2,1));
    ttp_wfe = sqrt(wfeTT.^2 + wfeP.^2);

    figure(1)
    h1_ = plot(t, ttp_wfe); hold on;
    xlabel('Time (s)'); ylabel('WFE (nm)'); grid on; axis tight;
    fprintf("TTP induced WFE: %.6g nm\n", rms(ttp_wfe,2));

    figure(2)
    h2_ = plot(t, rms(m1HPlc_dt,2)); hold on;
    xlabel('Time (s)'); ylabel('M1 HP-LC rms (N)'); grid on; axis tight;
    fprintf("M1-HP rms: %.6g N\n", rms(rms(m1HPlc_dt,2),1));
end

legend('No M1 act damping','linear damping','quad damping')
hold off;
figure(1); hold off;


