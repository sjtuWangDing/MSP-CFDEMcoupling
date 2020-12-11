#################################################################
## SETTINGS FOR 5.x                                            ##
#################################################################

# Specify additional include and library paths, as well as libraries for the compilation
CFDEM_ADD_INCLUDE = \
CFDEM_ADD_LIB = \
CFDEM_ADD_LIB_PATH = \

# additional static libraries to be linked to lagrangian library
CFDEM_ADD_STATIC_LIB = \

CFDEM_TRI_SURF = \
  -ltriSurface

CFDEM_SPRAY_LIBS = \
  -lthermophysicalProperties \

#----------------------------------------------------------------
# incompressible turbulence model settings
#----------------------------------------------------------------

# paths for incompressible turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_PATHS = \
  -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
  -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
  -I$(LIB_SRC)/fvOptions/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_LIBS = \
  -lturbulenceModels \
  -lincompressibleTurbulenceModels \
  -lfvOptions \

#----------------------------------------------------------------
# compressible turbulence model settings
#----------------------------------------------------------------

# paths for compressible turbulence models to use
CFDEM_ADD_COMPTURBMOD_PATHS = \
  -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
  -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
  -I$(LIB_SRC)/transportModels/compressible/lnInclude \
  -I$(LIB_SRC)/thermophysicalModels/radiation/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_COMPTURBMOD_LIBS = \
  -lturbulenceModels \
  -lcompressibleTurbulenceModels \