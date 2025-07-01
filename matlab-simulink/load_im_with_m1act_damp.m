%
% Script to load the Simulink model parameters. The model includes the
% structural dynamics, the mount, the M1 actuator force and positon loops.
% This model provides an implementation of the M1 actuator damping.
% 
%#ok<*UNRCH>

simulink_fname = "im_with_m1act_damp";


%% General IM settings
%%
sysDamp = 0.02;     % Set telescope structural dynamics damping
FEM_Ts = 1/1e3;     % Set structural dynamics model sampling period

% Zenith angle string (affects the control model for now)
sZa = "30";

% - - - - - Simulation setting flags (0:disables) - - - - -
clear osim
osim.reduce_model = 0;  % DO NOT ENABLE UNTIL TODO IS DONE![bool] model reduction feature
osim.dc_mm_comp = 1;    % [bool] DC mismatch compensation
osim.bpless_wl = 1;     % [bool] Smoothing wind load
osim.wload_en = 1;      % [bool] Wind load input
% MNT
osim.mntC_en = 1;       % [bool] Mount feedback control loop switch
osim.en_servo = 0;      % [bool] Enable/disable mount trajectory
osim.mnt_FFc_en = 1;    % [bool] Azimuth feedforward action switch      
% M1
osim.m1olf_en = 1;      % [bool] M1 outer force loop switch
osim.m1cpp_en = 0;      % [bool] M1 Command pre-processor activation flag
osim.m1act_damp = 0;    % [numeric] 0: no damping, 1: linear, and 2: quadratic
% % M2
% osim.m2PZT_en = 0;      % [bool] M2 PZT control loop switch
% osim.m2_pos_en = 0;     % [bool] M2 Positioner control loop


%% Telescope structural dynamics model
%%

% mfolder = '/home/rromano/Workspace/gmt-data';
mFolder = '/home/rromano/mnt';
% ModelID = "20250516_1420_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_M1_actDamping";
ModelID = "20250506_1715_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_M1_actDamping";
fName = "modal_state_space_model_2ndOrder.mat";
    
if(~exist('inputTable','var') || 0)
    try
        fprintf('Loading structural model %s from\n%s\n', fName, ModelID);
        load(fullfile(mFolder, ModelID, fName),...
            'inputs2ModalF','modalDisp2Outputs',...
            'eigenfrequencies','proportionalDampingVec',...
            'inputTable','outputTable');
    catch
        error("Unable to load structural model files from \n%s"+...
            "Download files from AWS-S3 (or DROBO).",ModelID);
    end
    % Static solution gain matrix
    staticSolFile = fullfile(mFolder, ModelID,"static_reduction_model.mat");
    try
        load(staticSolFile,'gainMatrixMountControlled');
        gainMatrix = gainMatrixMountControlled;
    catch
        load(staticSolFile,'gainMatrix');
    end
    fprintf('Static gain matrix loaded from\n%s\n', staticSolFile);
end


%% Pick model IOs according to inputTable and outputTables
%%
% INPUTS
desiredInputLabels = [...
    "OSS_ElDrive_Torque"; "OSS_AzDrive_Torque"; "OSS_RotDrive_Torque";...
    "OSS_Harpoint_delta_F"; "M1_actuators_segment_1";...
    "M1_actuators_segment_2"; "M1_actuators_segment_3";...
    "M1_actuators_segment_4"; "M1_actuators_segment_5";...
    "M1_actuators_segment_6"; "M1_actuators_segment_7";...
    "OSS_M1_lcl_6F";"MC_M2_lcl_6F";"CFD_202504_6F"; "OSS_Hardpoint_extension";];
isDesired = zeros(size(inputs2ModalF,2),1);
modelMuxDims = zeros(numel(desiredInputLabels),1);
first_idx_vec = zeros(numel(desiredInputLabels),1);

for i1 = 1:numel(desiredInputLabels)
    aux = inputTable{desiredInputLabels{i1},"indices"}{1}(:);    
    isDesired(aux) = 1;
    modelMuxDims(i1) = length(aux);
    first_idx_vec(i1) = aux(1); 
end
if any(diff(first_idx_vec) < 0)
    error("Labels in desiredInputLabels must follow the inputTable ordering!");
end
indDesInputs = find(isDesired ~= 0);
modelMuxDims(modelMuxDims == 0) = [];

% OUTPUTS
desiredOutputLabels = [...
    "OSS_ElEncoder_Angle"; "OSS_AzEncoder_Angle"; "OSS_RotEncoder_Angle";...
    "OSS_Hardpoint_D"; "M1_actuators_segment_1_D"; "M1_actuators_segment_2_D";...
    "M1_actuators_segment_3_D";"M1_actuators_segment_4_D";...
    "M1_actuators_segment_5_D";"M1_actuators_segment_6_D";...
    "M1_actuators_segment_7_D";"OSS_M1_lcl";"MC_M2_lcl_6D";"OSS_Hardpoint_force"];
isDesired = zeros(size(modalDisp2Outputs,1),1);
modelDemuxDims = zeros(numel(desiredOutputLabels),1);
first_idx_vec = zeros(numel(desiredOutputLabels),1);

for i1 = 1:numel(desiredOutputLabels)
    aux = outputTable{desiredOutputLabels(i1),"indices"}{1}(:);
    if(contains(desiredOutputLabels(i1),"M1_actuators_segment_"))
        isDesired(aux) = 2;
    else, isDesired(aux) = 1;
    end
    modelDemuxDims(i1) = length(aux);
    first_idx_vec(i1) = aux(1);
end
if any(diff(first_idx_vec) < 0)
    error("Labels in desiredOutputLabels must follow the outputTable ordering!");
end
indDesOutputs = find(isDesired ~= 0);
ind_M1act_vel = find(isDesired == 2);
modelDemuxDims(modelDemuxDims == 0) = [];

%
% REMARK: Use
% utils.createStructDynGMTslx(desiredInputLabels,modelMuxDims,...
% desiredOutputLabels,modelDemuxDims,'GMT')
% to create the structural model block with the selected IOs.
%


%% Structural model discretization
%%
if osim.reduce_model
    % Compute the approximate Hankel singular values
    % <-- TO DO: Always keep the first 3 modes
    [gamma,~] = utils.approxhsv(sqrt(om2), sysDamp*ones(size(om2)), phiB, phiC);
    [~,si] = sort(gamma,'descend');
    th1 = 1e-7;
    gammamax = max(gamma(:,1+3));
    
    nr = length(find((gamma./gammamax) >= th1));
    warning('\n-> The number of modes for TH=%.2g is %d\n',th1,nr);
    mode_ind_vec = si(1:nr);
else
    nr = length(eigenfrequencies(:));
    mode_ind_vec = 1:nr;
%     mode_ind_vec = 4:nr;
end

% om^2 and 2*zeta*om vectors 
om2 = (2*pi*eigenfrequencies(mode_ind_vec)).^2;
twice_zom = 2*sysDamp .* sqrt(om2);

fprintf('Plant model damping ratio set to:%.2g\n',(twice_zom(4)/2/sqrt(om2(4))))
% Perform discretization and provide the 2nd order form DT simulation parameters
PG = zeros(length(om2),6);
for i = 1:length(om2)
    PhiGamma = expm([0 1 0; -om2(i) -twice_zom(i) 1; 0 0 0]*FEM_Ts);
    PG(i,:) = [PhiGamma(1,1:3), PhiGamma(2,1:3)];
end

phiB = inputs2ModalF(mode_ind_vec,indDesInputs);
phiC = modalDisp2Outputs(indDesOutputs,mode_ind_vec);
phiCvel = modalDisp2Outputs(ind_M1act_vel,mode_ind_vec);

s_rate_msg = 'Telescope structural model sampling rate set to %gkHz\n';
if((FEM_Ts ~= 1/8e3) && (FEM_Ts ~= 1e-3)), warning(s_rate_msg,1/FEM_Ts/1e3);
else, fprintf(s_rate_msg,1/FEM_Ts/1e3);
end

%% Load Mount parameters and controller models
%%
try
    % Load ODC mount controller and driver parameters from ODC files
    o = get_odc_mnt_dt(sZa);
catch
    error("Unable to load ODC data! "+...
        "If functions are not available try loading from a standalone file.\n");
end

% Mount (FB & FF) controller discretization
c2d_opt = c2dOptions('method','tustin');
mnt_TF_Ts = FEM_Ts;
% AZ FB & FF
mount.az.SSdtHfb = balreal(c2d(ss(o.az.c.Hp), mnt_TF_Ts, c2d_opt));
mount.az.SSdtHff = balreal(c2d(ss(o.az.c.Hff), mnt_TF_Ts, c2d_opt));
% EL FB & FF
mount.el.SSdtHfb = balreal(c2d(ss(o.el.c.Hp), mnt_TF_Ts, c2d_opt));
mount.el.SSdtHff = balreal(c2d(ss(o.el.c.Hff), mnt_TF_Ts, c2d_opt));
% GIR FB & FF
mount.gir.SSdtHfb = balreal(c2d(ss(o.gir.c.Hp), mnt_TF_Ts, c2d_opt));
mount.gir.SSdtHff = balreal(c2d(ss(o.gir.c.Hff), mnt_TF_Ts, c2d_opt));

% Mount driver model parameters
mount.delay = 4.0e-3;     % [s] DRV delay: 4ms (GMT25-ANA-40000-0007-2.0 - Pg26)
drv_delay = ceil(mount.delay/FEM_Ts);
fprintf('Driver delay set to %d sampling periods (Ts=%gms).\n',drv_delay,1e3*FEM_Ts);
if(rem(mount.delay,FEM_Ts))
    warning('Driver delay is not a multiple of the sampling period!')
end

% Current loop dynamics tranfer function
Hdrive_d = c2d(o.Hdrive(1,1), FEM_Ts, 'tustin');


%% Load M1 subsystems data (controllers and parameters)
%%

% File with controller and interface parameters
ctrl_fname = sprintf('controls_5pt1g_z%s_llTT_oad',sZa);
load(fullfile(ctrl_fname),'m1sys');
% Overwrite controller with Lead-lag+PI (as in test cell) and update SA
% dynamics.
m1sys = ovr_ofl_crtl(m1sys);

m1_act_lin_d = osim.m1act_damp*1800;    %[Ns/m] Linear M1 actuator damping
m1_act_quad_d = 9000;   %[Ns^2/m^2] Quadratic M1 actuator damping

% Upsilon maps the actuator forces into forces and moments about the mirror
% CG.
m1_dt_folder = [];
load(fullfile(m1_dt_folder,'OA_SupportActuatorArrayConfig'),'OA_Upsilon');
load(fullfile(m1_dt_folder,'CS_SupportActuatorArrayConfig'),'CS_Upsilon');
Upsilon = blkdiag(OA_Upsilon,OA_Upsilon,OA_Upsilon,OA_Upsilon,...
    OA_Upsilon,OA_Upsilon, CS_Upsilon);

[M_oa_xyz2abc,M_cs_xyz2abc,M_oa_abc2xyz,M_cs_abc2xyz] = calc_ort2cylTs(m1_dt_folder);


%% M1 Command pre-processor (CCP)
%%
% Velocity and acceleration limits -- Symmetric bounds, i.e.: vmin = -vmax
v_max = 50e-6;      % [m/s]
a_max = 250e-6;     % [m/s^2]

% Command pre-processor (CPP) settings

delta_sim = 1e-3;%m1sys{1}.ofl.Ts;	% Command demand sampling period
Tcmd = 12;                      % Command update period

% Bessel filter to limit acc derivative (jerk)
bessel_num = [9.50565136e-07, 3.80226054e-06, 5.70339082e-06, 3.80226054e-06, 9.50565138e-07];
bessel_den = [ 1., -3.80206663, 5.42365878, -3.4403268, 0.81874986];
bessel_f_disc = d2d(tf(bessel_num, bessel_den, 0.01), delta_sim);

% Bessel filter state-space realization
[A_f, B_f, C_f, D_f] = ssdata(balreal(bessel_f_disc));
[xf, ~] = step(ss(A_f, B_f, eye(4), zeros(4,1), delta_sim),10);
x0_f = xf(end,:);

t_gd = 0*0.5;
t_gd_v = 0*t_gd;

% M1 target trajectory
t  = linspace(0, 30, 1001)';
% Switching times
t_s1 = Tcmd; t_s2 = 2*Tcmd; t_s3 = 3*Tcmd;
i_s1_ = find(t < t_s1, 1, 'last');
i_s2_ = find(t < t_s2, 1, 'last');
i_s3_ = find(t < t_s3, 1, 'last');

tj_mode = 'static';%'piston';%'normal';% uniform
rng(1);     % Set RNG seed
switch tj_mode
    case 'static'
        r = zeros(4,6);
    case 'normal'
        mean_r = 0;
        std_r = 2e-5;
        r = mean_r + std_r.*randn(4,6);
    case 'piston'
        r = zeros(4,6);
        r(:,3) = 1e-6;
    otherwise
        r_low = -2e-5;
        r_high = 2e-5;
        r = r_low + (r_high-r_low).*rand(4,6);
end

x_Si = [kron(r(1,:),ones(size(t(1:i_s1_-1))));...
    kron(r(2,:),ones(size(t(i_s1_:i_s2_-1))));...
    kron(r(3,:),ones(size(t(i_s2_:i_s3_-1))));...
    kron(r(4,:),ones(size(t(i_s3_:end))))];
m1_pos_cmd = [t, x_Si, zeros(size(t,1),36)];

% v = zeros(size(x));

n_in = 6;
% Initial position and velocity state
x0 = 0; % For now, one is unable to set a different IC for each HP length.
v0 = 0;


%% Open Simulink file
%%
open(simulink_fname);


%% Static gain mismatch compensation
%%

struct_dyn_label = "/Telescope model/Structural Dynamics GMT";
memory_label = simulink_fname+struct_dyn_label+"/Psi_ss_memory";
matmult_label = simulink_fname+struct_dyn_label+"/Psi_ss";
zoh_label = simulink_fname+struct_dyn_label+"/Psi_ssZOH";
rt_label = simulink_fname+struct_dyn_label+"/Psi_ssRT";

if osim.dc_mm_comp
    try
        K_ss = phiC(:,4:end)* diag(1./((2*pi*eigenfrequencies(4:end)).^2))* phiB(4:end,:);        
        Psi_ss = gainMatrix(indDesOutputs,indDesInputs) - K_ss;
        
        noDCMCOutputLabels = [...
            "OSS_ElEncoder_Angle"; "OSS_AzEncoder_Angle"; "OSS_RotEncoder_Angle";...
            "M1_actuators_segment_1_D"; "M1_actuators_segment_2_D";...
            "M1_actuators_segment_3_D";"M1_actuators_segment_4_D";...
            "M1_actuators_segment_5_D";"M1_actuators_segment_6_D";...
            "M1_actuators_segment_7_D"];
        y_idx = [[1;cumsum(modelDemuxDims(1:end-1))+1], cumsum(modelDemuxDims)];
        for iy = 1:numel(noDCMCOutputLabels)
            idx = find(contains(desiredOutputLabels,noDCMCOutputLabels(iy)));    
            Psi_ss(y_idx(idx,1):y_idx(idx,2), :) = 0;
            if false
                fprintf("Zeroing %d DCMC rows corresponding to output %s\n",...
                    length(y_idx(idx,1):y_idx(idx,2)), noDCMCOutputLabels(iy)); 
            end
        end
                
        Psi_ssTs = FEM_Ts;        
        set_param(memory_label,'Commented','off');
        set_param(matmult_label,'Commented','off');
        set_param(zoh_label,'Commented','off');
        set_param(rt_label,'Commented','off');
    catch
        warning('Unable to compute static compensation matrix.');
        warning('Disabling static compensation matrix.');
        set_param(memory_label,'Commented','on');
        set_param(matmult_label,'Commented','on');
        set_param(zoh_label,'Commented','on');
        set_param(rt_label,'Commented','on');
        StaticModelFolder = [];
    end
else
    set_param(memory_label,'Commented','on');
    set_param(matmult_label,'Commented','on');
    set_param(zoh_label,'Commented','on');
    set_param(rt_label,'Commented','on');
    StaticModelFolder = [];
end



%% Wind load time series
%%

wl_demux = [42,42, 468]; %#ok<*NASGU>

if(osim.wload_en)
    ModelFolder = fullfile(mFolder, ModelID);
    
    [IMLoads, wl_demux] = load_WLdt(ModelFolder,60);
    set_param(simulink_fname+"/Wind Loads",'Commented','off');
else
    set_param(simulink_fname+"/Wind Loads",'Commented','on');
end




%% Auxiliary functions
%%

%% Function to overwrite M1 controller parameters
function m1sys = ovr_ofl_crtl(m1sys)  

% M1 outer force loop (OFL) controller matrix
ofl.Ts = 1/100;  % OFL sampling time

% OFL controller --- PI+Lead-lag (SPIE2024-)
fPI = tf([0.1 5],[1 0]);
flead = tf([1 0.1*2*pi],[1 0.3*2*pi]);
flag = tf([1 15*2*pi],[1 5*2*pi]);
fbH = fPI * flead * flag;

oflC_ss = cell(6,1);

% Lag term to approximate nonlinear pneumatic circuit behavior
HlagDT = c2d(0.2/1 *tf([1 1*2*pi],[1 0.2*2*pi]), m1sys{1}.SAdyn.Ts, 'tustin');

for seg = 1:7
    % Controller channel loop
    for ich = 1:numel(oflC_ss)
        oflC_ss{seg} = balreal(ss(fbH));
        m1sys{seg}.ofl.SSdtC{ich} = c2d(oflC_ss{seg},ofl.Ts,'foh');

        % Display BODE plot to assess the effect of the discretization
        if(false && seg == 1)
            plot_labels = {'F_x','F_y','F_z','M_x','M_y','M_z'};
            hbode = bodeoptions;
            hbode.FreqUnits = 'Hz';
            hbode.XLabel.FontSize = 11;
            hbode.YLabel.FontSize = 11;
            hbode.Ylabel.String{1} = 'Mag';
            hbode.Title.FontSize = 12;
            figure(7000+seg)
            subplot(2,3,ich)
            bode(oflC_ss{seg},hbode); hold on;
            bode(m1sys{seg}.ofl.SSdtC{ich},'r--',hbode); hold off;
            grid on; title(sprintf('OFL %s Controller',plot_labels{ich}));
            xlim([0.05, 0.5/ofl.Ts])
            if ich == numel(oflC_ss)
                set(gcf,'Position',[600   567   1.6*304*3/2   420*0.95]);
                legend('CT','DT','Location','northwest'); legend box off;
            end
        end
    end
    % Update support actuator dynamics
    m1sys{seg}.SAdyn = m1sys{seg}.SAdyn*HlagDT;
end

end

%% Function to get ODC mount controller parameters
function o = get_odc_mnt_dt(sZa)

% oTest.sZa: elevation zenith angle (ZA) as string e.g. '00','30','60'
oTest.sZa = sZa; %'30', %'00' or '60'
% oTest.sVer: FEM version as string e.g. '19'
oTest.sVer = '20';
% oTest.sSubVer: FEM subversion as string e.g. '1'
oTest.sSubVer = '11'; %'2'; %'9';
% oTest.sDamping: now '02' means 2% structural dumping
oTest.sDamping ='02';
% oTest.bUseReducedModel: [true|false] if true: a reduced model is used
%  which was created by the balred method
oTest.bUseReducedModel = true;

odc_file_folder = '/home/rromano/Workspace/gmt-mnt-odc';
odc_main_folder = "fdr2023/MatlabFilesE2E_2023-05-10";
odc_base_util_folder = fullfile(odc_file_folder,odc_main_folder,'base/util');
odc_base_conf_folder = fullfile(odc_file_folder,odc_main_folder,'base/conf');
addpath(odc_base_util_folder, odc_base_conf_folder);
fprintf('+Including folder\n%s\ninto MatLab path.\n',...
    fullfile(odc_file_folder,odc_main_folder));
fprintf('Getting ODC mount model parameters ...\n');

% ODC Simulink model used (located in ../base)
oTest.sRoot = 'root';
% oTest.sHbsConf: name of HBS configuration: e.g. 'HbTp19'
oTest.sHbsConf = 'HaTp19'; %'HbTp19'
% oTest.sViscFrCase: one of 3 cases w.r.t. viscosity: ['ViscFrLow', 'ViscFrMedium', 'ViscFrHigh']  see fun_mueByTempAzEl
oTest.sViscFrCase = 'ViscFrLow'; %lowest viscosity=>lowest damping
% oTest.sModelDirIn: directory relative to ../ss_model/ where the state space models are located
% returns structure [o] with all configuration parameters
oTest.sModelDirIn = 'v20.11/n100HzR800';

o = fun_confBase(oTest);
% Remove folders from Matlab path
fprintf('-Removing folders\n%s\n%s\nfrom MatLab path.\n',...
    odc_base_util_folder, odc_base_conf_folder);
rmpath(odc_base_util_folder, odc_base_conf_folder);

end

%% Function to calculate M1 actuator transformations
function [M_oa_xyz2abc, M_cs_xyz2abc, M_oa_abc2xyz, M_cs_abc2xyz] =...
    calc_ort2cylTs(m1_dt_folder)

% OA
load(fullfile(m1_dt_folder,'OA_SupportActuatorArrayConfig'),...
    'OA_ActData','xyz2abc40','xyz2abc41','xyz2abc5',...
    'abc2xyz40','abc2xyz41','abc2xyz5');

nr = numel(find(OA_ActData(:,5) >= 40))*3 +...
    numel(find(OA_ActData(:,5) > 5 & (OA_ActData(:,5) < 40))) +...
    numel(find(OA_ActData(:,5) == 5))*6 +...
    numel(find(OA_ActData(:,5) < 5));
nc = nr - numel(find(OA_ActData(:,5) == 5))*3;
fprintf("Dimension of Txyz2abc for OA segment: %dx%d\n",nr,nc);
M_oa_xyz2abc = zeros(nr,nc);
M_oa_abc2xyzT = zeros(nr,nc);
i_cols = 1;
i_rows = 1;
for i_act = 1:size(OA_ActData,1)

    if(OA_ActData(i_act,5) == 40)
        M_oa_xyz2abc((0:2)+i_rows, (0:2)+i_cols) = 1/3*xyz2abc40;
        M_oa_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = abc2xyz40';
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(OA_ActData(i_act,5) == 41)
        M_oa_xyz2abc((0:2)+i_rows, (0:2)+i_cols) = 1/3*xyz2abc41;
        M_oa_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = abc2xyz41';
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(OA_ActData(i_act,5) == 5)
        M_oa_xyz2abc((0:5)+i_rows, (0:2)+i_cols) = 1/6*xyz2abc5;
        M_oa_abc2xyzT((0:5)+i_rows, (0:2)+i_cols) = abc2xyz5';
        i_cols = i_cols+3;
        i_rows = i_rows+6;
    else
        M_oa_xyz2abc(i_rows,i_cols) = 1;
        M_oa_abc2xyzT(i_rows,i_cols) = 1;
        i_cols = i_cols+1;
        i_rows = i_rows+1;
    end
end

M_oa_abc2xyz = M_oa_abc2xyzT';

% CS
load(fullfile(m1_dt_folder,'CS_SupportActuatorArrayConfig'),...
    'CS_ActData','xyz2abc40','xyz2abc41','xyz2abc5',...
    'abc2xyz40','abc2xyz41','abc2xyz5');

nr = numel(find(CS_ActData(:,5) >= 40))*3 +...
    numel(find(CS_ActData(:,5) > 5 & (CS_ActData(:,5) < 40))) +...
    numel(find(CS_ActData(:,5) == 5))*6 +...
    numel(find(CS_ActData(:,5) < 5));
nc = nr - numel(find(CS_ActData(:,5) == 5))*3;
fprintf("Dimension of Txyz2abc for CS segment: %dx%d\n",nr,nc);
M_cs_xyz2abc = zeros(nr,nc);
M_cs_abc2xyzT = zeros(nr,nc);
i_cols = 1;
i_rows = 1;
for i_act = 1:size(CS_ActData,1)

    if(CS_ActData(i_act,5) == 40)
        M_cs_xyz2abc((0:2)+i_rows, (0:2)+i_cols) = 1/3*xyz2abc40;
        M_cs_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = abc2xyz40';
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(CS_ActData(i_act,5) == 41)
        M_cs_xyz2abc((0:2)+i_rows, (0:2)+i_cols) = 1/3*xyz2abc41;
        M_cs_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = abc2xyz41';
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(CS_ActData(i_act,5) == 5)
        M_cs_xyz2abc((0:5)+i_rows, (0:2)+i_cols) = 1/6*xyz2abc5;
        M_cs_abc2xyzT((0:5)+i_rows, (0:2)+i_cols) = abc2xyz5';
        i_cols = i_cols+3;
        i_rows = i_rows+6;
    else
        M_cs_xyz2abc(i_rows,i_cols) = 1;
        M_cs_abc2xyzT(i_rows,i_cols) = 1;
        i_cols = i_cols+1;
        i_rows = i_rows+1;
    end
end

M_cs_abc2xyz = M_cs_abc2xyzT';

end

%% Function to load wind-load time series
function [windload_dt, wl_demux] = load_WLdt(ModelFolder, dur)

dtin_path = '/home/rromano/Workspace/dos-actors/clients/windloads';
dtin_file = fullfile(dtin_path,'model_data_1.parquet');
try
    parquetINFO = parquetinfo(dtin_file);
    cfd_in = parquetread(dtin_file, "SampleRate",1e3,...
        "SelectedVariableNames", parquetINFO.VariableNames);
    % "CFDMountWindLoads"    "CFDM1WindLoads"    "CFDM2WindLoads"

    n_mntdist = size(cfd_in.CFDMountWindLoads{1},1);
    mountwl = reshape(cell2mat(cfd_in.CFDMountWindLoads),n_mntdist,[]);
    m1wl = reshape(cell2mat(cfd_in.CFDM1WindLoads),42,[]);
    m2wl = reshape(cell2mat(cfd_in.CFDM2WindLoads),42,[]);
    time = seconds(cfd_in.Time);
catch
    error('Unable to load wind load time series from file \n%s\n',dtin_file);
end

fHz_wl = floor(1/diff(time(1:2)));
fprintf("Wind load time series (%dHz) loaded from\n%s\n",fHz_wl,dtin_file);

% Handle model extra inputs
% Load inputTable
load(fullfile(ModelFolder,'modal_state_space_model_2ndOrder.mat'),'inputTable');

% List of identifiers of the inputs to be disregarded
wl_descr = ["M2 cell";...
    "Ap instrument";...
    "Cable Trays";...
    "Cabs and instavol";...
    "Azimuth disk";...
    "IP platform with IP";...
    "Lab station volume";...
    "instrument volume attached";...
    "Laser Guide Star cabinet"];

wl_sub_inds = zeros(length(wl_descr),2);

% Take the starting index (first column) and the size of those inputs
for i1 = 1:length(wl_descr)
    aux = find(startsWith(inputTable.descriptions{'CFD_202504_6F'}(:),wl_descr(i1)));
    wl_sub_inds(i1,:) = [aux(1), length(aux)];
end

[~,i__] = sort(wl_sub_inds);
% Insert 0 at columns of the time series related to disregarded model inputs
for i1 = 1:length(wl_descr)
    mountwl = [mountwl(1:wl_sub_inds(i__(i1),1)-1,:);...
        zeros(wl_sub_inds(i__(i1),2),size(mountwl,2));...
        mountwl(wl_sub_inds(i__(i1),1):end,:)];
end

% Take data corresponding to the last Xsec (X=120).
N = size(time,1);
tdist_idxs = N-(dur*fHz_wl)+1:N;
windload_dt.time = time(tdist_idxs);
windload_dt.signals.values = ...
    [m1wl(:,tdist_idxs); m2wl(:,tdist_idxs); mountwl(:,tdist_idxs)]';
% Enforce first time instant to zero
windload_dt.time = windload_dt.time - windload_dt.time(1);

wl_demux = [42, 42, size(mountwl,1)];

end
