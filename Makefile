LIBRARY := SUHH2TopTagging
#DICT := include/AnalysisModuleRunnerTopTagging.h include/SUHH2TopTagging_LinkDef.h
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -lTMVA -lSUHH2HOTVR
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.local 
include ../Makefile.common
