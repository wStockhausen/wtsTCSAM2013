#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=clang
CCC=clang++
CXX=clang++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=CLang-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug-Mac
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/FisheryData.o \
	${OBJECTDIR}/ModelConfiguration.o \
	${OBJECTDIR}/ModelConstants.o \
	${OBJECTDIR}/ModelData.o \
	${OBJECTDIR}/TCSAM_WTSv01.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-O3 -DSAFE_ALL -Wno-deprecated
CXXFLAGS=-O3 -DSAFE_ALL -Wno-deprecated

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=-static

# Link Libraries and Options
LDLIBSOPTIONS=-L/Users/WilliamStockhausen/Programs/admb/build/dist/lib -L/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/lib -L/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/dist/Debug-MacOSX/CLang-MacOSX -lwtsadmb -ladmb -lcontrib -lwtsadmb -ladmb -lcontrib

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2013

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2013: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2013 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/FisheryData.o: FisheryData.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LAPLACE -D__GNUDOS__ -Dlinux -I. -I/Users/WilliamStockhausen/Programs/admb/build/dist/include -I/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/include -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FisheryData.o FisheryData.cpp

${OBJECTDIR}/ModelConfiguration.o: ModelConfiguration.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LAPLACE -D__GNUDOS__ -Dlinux -I. -I/Users/WilliamStockhausen/Programs/admb/build/dist/include -I/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/include -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ModelConfiguration.o ModelConfiguration.cpp

${OBJECTDIR}/ModelConstants.o: ModelConstants.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LAPLACE -D__GNUDOS__ -Dlinux -I. -I/Users/WilliamStockhausen/Programs/admb/build/dist/include -I/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/include -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ModelConstants.o ModelConstants.cpp

${OBJECTDIR}/ModelData.o: ModelData.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LAPLACE -D__GNUDOS__ -Dlinux -I. -I/Users/WilliamStockhausen/Programs/admb/build/dist/include -I/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/include -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ModelData.o ModelData.cpp

${OBJECTDIR}/TCSAM_WTSv01.o: TCSAM_WTSv01.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LAPLACE -D__GNUDOS__ -Dlinux -I. -I/Users/WilliamStockhausen/Programs/admb/build/dist/include -I/Users/WilliamStockhausen/Programs/admb/build/dist/contrib/include -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TCSAM_WTSv01.o TCSAM_WTSv01.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam2013

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
