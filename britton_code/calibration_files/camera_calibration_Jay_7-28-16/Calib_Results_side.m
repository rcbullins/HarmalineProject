% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 1299.190060612994785 ; 1279.855705160334310 ];

%-- Principal point:
cc = [ -72.414531199228549 ; -37.673462686674526 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -1.282092454762115 ; 2.944347670543990 ; 0.080215785247430 ; 0.042416124438176 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 104.576978474993084 ; 77.921247419338997 ];

%-- Principal point uncertainty:
cc_error = [ 0.000000000000000 ; 0.000000000000000 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.309400695698750 ; 0.976075773189418 ; 0.035581508184888 ; 0.049772748650250 ; 0.000000000000000 ];

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
omc_1 = [ 1.672330e+00 ; 1.655379e+00 ; 5.473304e-01 ];
Tc_1  = [ 2.068291e+01 ; 1.614209e+01 ; 1.456984e+02 ];
omc_error_1 = [ 2.637997e-02 ; 2.358328e-02 ; 4.138792e-02 ];
Tc_error_1  = [ 9.542373e-01 ; 1.996486e-01 ; 7.260123e+00 ];

%-- Image #2:
omc_2 = [ 1.649430e+00 ; 1.637574e+00 ; 5.722038e-01 ];
Tc_2  = [ 1.554591e+01 ; 1.732126e+01 ; 1.478418e+02 ];
omc_error_2 = [ 2.512243e-02 ; 2.308138e-02 ; 3.809466e-02 ];
Tc_error_2  = [ 7.306767e-01 ; 2.523391e-01 ; 7.590121e+00 ];

%-- Image #3:
omc_3 = [ 1.655158e+00 ; 1.018452e+00 ; 2.797202e-01 ];
Tc_3  = [ 2.466141e+01 ; 1.813283e+01 ; 1.398410e+02 ];
omc_error_3 = [ 1.063351e-02 ; 1.678328e-02 ; 3.973181e-02 ];
Tc_error_3  = [ 8.599475e-01 ; 6.003564e-01 ; 6.344097e+00 ];

%-- Image #4:
omc_4 = [ 2.062502e+00 ; 1.937406e+00 ; 2.823522e-01 ];
Tc_4  = [ 2.470507e+01 ; 1.364972e+01 ; 1.503945e+02 ];
omc_error_4 = [ 6.641193e-02 ; 4.675591e-02 ; 9.694788e-02 ];
Tc_error_4  = [ 1.144723e+00 ; 1.513895e-01 ; 7.719038e+00 ];

%-- Image #5:
omc_5 = [ 2.056807e+00 ; 1.635194e+00 ; 5.189116e-01 ];
Tc_5  = [ 2.463575e+01 ; 1.305752e+01 ; 1.421979e+02 ];
omc_error_5 = [ 4.446065e-02 ; 2.852549e-02 ; 6.387457e-02 ];
Tc_error_5  = [ 1.131653e+00 ; 1.185396e-01 ; 7.038459e+00 ];

%-- Image #6:
omc_6 = [ 1.799103e+00 ; 1.444691e+00 ; 7.956464e-01 ];
Tc_6  = [ 3.379032e+01 ; 8.114720e+00 ; 1.334583e+02 ];
omc_error_6 = [ 2.640292e-02 ; 2.213836e-02 ; 3.381529e-02 ];
Tc_error_6  = [ 1.467171e+00 ; 2.649528e-01 ; 6.384753e+00 ];

%-- Image #7:
omc_7 = [ 1.666789e+00 ; 1.063471e+00 ; 4.729190e-01 ];
Tc_7  = [ 2.525218e+01 ; 1.228816e+01 ; 1.363399e+02 ];
omc_error_7 = [ 1.211807e-02 ; 1.703167e-02 ; 3.636474e-02 ];
Tc_error_7  = [ 1.000795e+00 ; 2.972873e-01 ; 6.050522e+00 ];

