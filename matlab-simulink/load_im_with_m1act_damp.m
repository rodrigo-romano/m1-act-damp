%
% 
% 

simulink_fname = "im_with_m1act_damp";
% simulink_fname = "im_with_m1act_damp_CT";

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
osim.m1olf_en = 0;      % [bool] M1 outer force loop switch
osim.m1cpp_en = 0;      % [bool] M1 Command pre-processor activation flag
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
    "M1_actuators_segment_6"; "M1_actuators_segment_7"; "OSS_M1_lcl_6F";...
    "CFD_202504_6F"; "OSS_Hardpoint_extension";];
isDesired = zeros(size(inputs2ModalF,2),1);
modelMuxDims = zeros(numel(desiredInputLabels),1);

for i1 = 1:numel(desiredInputLabels)
    aux = inputTable{desiredInputLabels{i1},"indices"}{1}(:);
    isDesired(aux) = 1;
    modelMuxDims(i1) = length(aux);
end
indDesInputs = find(isDesired ~= 0);
modelMuxDims(modelMuxDims == 0) = [];

% OUTPUTS
desiredOutputLabels = [...
    "OSS_ElEncoder_Angle"; "OSS_AzEncoder_Angle"; "OSS_RotEncoder_Angle";...
    "OSS_Hardpoint_D"; "M1_actuators_segment_1_D"; "M1_actuators_segment_2_D";...
    "M1_actuators_segment_3_D";"M1_actuators_segment_4_D";...
    "M1_actuators_segment_5_D";"M1_actuators_segment_6_D"; "M1_actuators_segment_7_D";...
    "OSS_M1_lcl";"OSS_Hardpoint_force"];
isDesired = zeros(size(modalDisp2Outputs,1),1);
modelDemuxDims = zeros(numel(desiredOutputLabels),1);

for i1 = 1:numel(desiredOutputLabels)
    aux = outputTable{desiredOutputLabels(i1),"indices"}{1}(:);
    if(contains(desiredOutputLabels(i1),"M1_actuators_segment_"))
        isDesired(aux) = 2;
    else, isDesired(aux) = 1;
    end
    modelDemuxDims(i1) = length(aux);
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

%%

m1_dt_folder = [];
load(fullfile(m1_dt_folder,'OA_SupportActuatorArrayConfig'),'OA_Upsilon');
load(fullfile(m1_dt_folder,'CS_SupportActuatorArrayConfig'),'CS_Upsilon');
    Upsilon = blkdiag(OA_Upsilon,OA_Upsilon,OA_Upsilon,OA_Upsilon,...
        OA_Upsilon,OA_Upsilon, CS_Upsilon);


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
m1_act_damp = 1800; %[Ns/m] Actuator damping


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
        %
        num_mnt_ax_in = contains(desiredInputLabels,...
            ["OSS_ElDrive_Torque"; "OSS_AzDrive_Torque"; "OSS_RotDrive_Torque"]);
        v_in = [1; 1; 1; zeros(numel(desiredInputLabels)-3,1)];
        num_mnt_ax_out = contains(desiredOutputLabels,...
            ["OSS_ElEncoder_Angle"; "OSS_AzEncoder_Angle"; "OSS_RotEncoder_Angle"]);
        v_out = [1; 1; 1; zeros(numel(desiredOutputLabels)-3,1)];
        
        if(all(num_mnt_ax_in == v_in) && all(num_mnt_ax_out == v_out))
            Psi_ss(1:sum(outputTable.size(1:3)),:) = 0;
            Psi_ss(:,1:sum(inputTable.size(1:3))) = 0;
        else
            warning("No entries of Psi_ss set to zero!"+...
                "\nCheck for the chosen MNT IOs.\n");
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





%%

dtin_path = '/home/rromano/Workspace/dos-actors/clients/windloads';
dtin_file = fullfile(dtin_path,'model_data_1.parquet');
parquetINFO = parquetinfo(dtin_file);
cfd_in = parquetread(dtin_file, "SampleRate",1e3,...
        "SelectedVariableNames", parquetINFO.VariableNames);
    % "CFDMountWindLoads"    "CFDM1WindLoads"    "CFDM2WindLoads"


%% Auxiliary functions
%%

%% Function to get 
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