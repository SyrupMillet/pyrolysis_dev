# NGA location if not yet defined
NGA_HOME = /home/ubuntu/NGA_Simu/NGA2

# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_HYPRE = TRUE
USE_LAPACK= TRUE
USE_IRL   = FALSE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# External libraries are defined in .profile/.bashrc/.zshrc, but could be defined here as well
CVODE_INCLUDE_PATH = /usr/local/fortran
CVODE_LIB_PATH = /usr/local/lib

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := variable_density particles core data solver config grid libraries
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

# Inform user of Make.packages used
ifdef Ulocs
   $(info Taking user code from: $(Ulocs))
endif
$(info Taking base code from: $(Bdirs))

# F90FLAGS += -fcheck=array-temps
F90FLAGS += -cpp -g

fincludes += -I$(CVODE_INCLUDE_PATH)
LDFLAGS += -L$(CVODE_LIB_PATH)
libraries += -l:libsundials_fcvode_mod.a -l:libsundials_cvode.a 

# Target definition
all: $(executable)
	@echo COMPILATION SUCCESSFUL

# NGA compilation rules
include $(NGA_HOME)/tools/GNUMake/Make.rules

run:
	mpiexec -n 8 $(executable) -i input -v 2