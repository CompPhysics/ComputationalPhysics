%module ising2dim_backend
%{
#include "ising2dim_backend.h"
%}

%include "typemaps.i"
%apply double *OUTPUT {double& E_av}
%apply double *OUTPUT {double& E_variance}
%apply double *OUTPUT {double& M_av}
%apply double *OUTPUT {double& M_variance}
%apply double *OUTPUT {double& Mabs_av}

%include "ising2dim_backend.h"