programs = DeformableSolids2D DeformableSolids3D FrankeTest2D FrankeTest3D GeodesicsInHeat GradientDomainProcessing

COMPILER ?= gcc
#COMPILER ?= clang

# Allow "make -j" to operate in parallel over the programs.
all: $(programs)
$(programs):
	$(MAKE) -C $@ COMPILER=$(COMPILER)

programs_clean = $(foreach n,$(programs),clean_$(n))  # pseudo-dependency to allow "make -j" parallelism
clean: $(programs_clean)
$(programs_clean):
	$(MAKE) -C $(@:clean_%=%) clean

.PHONY: $(programs) $(programs_clean)
