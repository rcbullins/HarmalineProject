% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 1162.622827677088253 ; 1199.637894432873281 ];

%-- Principal point:
cc = [ 135.946155143267873 ; -75.647232904671114 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.544883027766914 ; 1.466821096073789 ; -0.020644233223987 ; 0.029111863755824 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 69.091285765432602 ; 85.875126627539771 ];

%-- Principal point uncertainty:
cc_error = [ 0.000000000000000 ; 0.000000000000000 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.374996544569132 ; 1.592723920530094 ; 0.051256591527764 ; 0.021170693070805 ; 0.000000000000000 ];

%-- Image size:
nx = 384;
ny = 260;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 7;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 0;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ -1.585955e+00 ; -1.833300e+00 ; 1.044319e-01 ];
Tc_1  = [ -2.365699e+00 ; 1.656761e+01 ; 1.457155e+02 ];
omc_error_1 = [ 9.453010e-03 ; 1.307160e-02 ; 2.660847e-02 ];
Tc_error_1  = [ 7.669856e-02 ; 3.053244e-01 ; 8.705895e+00 ];

%-- Image #2:
omc_2 = [ -1.595220e+00 ; -1.832954e+00 ; 7.945330e-02 ];
Tc_2  = [ -1.644503e-01 ; 1.442030e+01 ; 1.504977e+02 ];
omc_error_2 = [ 9.188647e-03 ; 1.227935e-02 ; 2.760197e-02 ];
Tc_error_2  = [ 5.425499e-02 ; 2.882133e-01 ; 9.063160e+00 ];

%-- Image #3:
omc_3 = [ -2.042977e+00 ; -2.199075e+00 ; 4.994475e-01 ];
Tc_3  = [ -7.456486e+00 ; 2.155566e+01 ; 1.439068e+02 ];
omc_error_3 = [ 4.425635e-02 ; 6.367382e-02 ; 5.093386e-02 ];
Tc_error_3  = [ 1.892144e-01 ; 3.267341e-01 ; 8.140377e+00 ];

%-- Image #4:
omc_4 = [ -1.302035e+00 ; -1.571419e+00 ; 4.465148e-01 ];
Tc_4  = [ 3.080653e+00 ; 1.578615e+01 ; 1.429894e+02 ];
omc_error_4 = [ 8.737508e-03 ; 9.761329e-03 ; 1.268176e-02 ];
Tc_error_4  = [ 7.198283e-02 ; 3.011033e-01 ; 8.404742e+00 ];

%-- Image #5:
omc_5 = [ -1.625043e+00 ; -1.561341e+00 ; 3.718582e-01 ];
Tc_5  = [ -5.797208e+00 ; 1.692927e+01 ; 1.420207e+02 ];
omc_error_5 = [ 9.802118e-03 ; 1.036772e-02 ; 1.837174e-02 ];
Tc_error_5  = [ 1.407845e-01 ; 3.143689e-01 ; 8.185903e+00 ];

%-- Image #6:
omc_6 = [ -1.889128e+00 ; -1.672670e+00 ; 5.893868e-02 ];
Tc_6  = [ -1.379177e+01 ; 1.900148e+01 ; 1.307869e+02 ];
omc_error_6 = [ 1.096028e-02 ; 1.235847e-02 ; 3.779372e-02 ];
Tc_error_6  = [ 2.531829e-01 ; 2.960637e-01 ; 7.446226e+00 ];

%-- Image #7:
omc_7 = [ -2.098691e+00 ; -2.054699e+00 ; 2.705200e-01 ];
Tc_7  = [ -1.229014e+01 ; 1.765015e+01 ; 1.398971e+02 ];
omc_error_7 = [ 3.536366e-02 ; 4.741158e-02 ; 7.175811e-02 ];
Tc_error_7  = [ 2.616653e-01 ; 2.945440e-01 ; 8.033336e+00 ];

