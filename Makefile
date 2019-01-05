
PHDLIBBASE=libraries

INCLUDES = \
	$(shell mpicc -showme:compile) \
	-I$(PHDLIBBASE)/log \
	-I$(PHDLIBBASE)/tracking

EXTRA_LIBS = \
	-L$(PHDLIBBASE)/tracking -ltracking \
	-L$(PHDLIBBASE)/log -llog


CXX = g++
CXXFLAGS = -c -g -Wall --std=c++11 $(INCLUDES)

#CXXFLAGS += -O3

INSTALL = install
INSTALLFLAGS = -D

LIBS = $(EXTRA_LIBS) \
	-lm \
	$(shell gsl-config --libs) \
	$(shell mpicxx -showme:link)

OBJS = rng.o \
	prior.o \
	sphericalprior.o \
	hierarchical_model.o \
	pathutil.o \
	attenuationexception.o 

MPI_OBJS = 


SRCS = Makefile \
	attenuationdataS2.hpp \
	attenuationexception.hpp \
	attenuationutil.hpp \
	birthgenericS2Voronoi.hpp \
	chainhistoryVoronoi.hpp \
	coordinate.hpp \
	deathgenericS2Voronoi.hpp \
	globalS2Voronoi.hpp \
	hierarchicalS2Voronoi.hpp \
	hierarchical_model.hpp \
	moveS2Voronoi.hpp \
	pathutil.hpp \
	perturbationS2Voronoi.hpp \
	perturbationcollectionS2Voronoi.hpp \
	prior.hpp \
	rng.hpp \
	sphericalprior.hpp \
	sphericalvoronoimodel.hpp \
	util.hpp \
	valueS2Voronoi.hpp \
	velocitymodel.hpp \
	attenuationexception.cpp \
	attenuationtomoS2Voronoi.cpp \
	attenuationtomoS2VoronoiPT.cpp \
	hierarchical_model.cpp \
	mksynthetic.cpp \
	pathutil.cpp \
	postS2Voronoi_likelihood.cpp \
	postS2Voronoi_mean.cpp \
	postS2Voronoi_mean_mpi.cpp \
	postS2Voronoi_text.cpp \
	prior.cpp \
	rng.cpp \
	sphericalprior.cpp

EXTRADIST = \
	libraries/Makefile \
	libraries/log.tar.gz \
	libraries/tracking.tar.gz \
	data/coreattenuation.txt \
	data/constant_0.1.txt \
	data/cubedsphere_0.1.txt \
	data/eastwest_0.1.txt \
	data/northsouth_0.1.txt \
	dataterrawulf/hierarchical_prior.txt \
	dataterrawulf/synthetic_hierarchical_prior.txt \
	dataterrawulf/pbs_singlechain_coredata_init.sh \
	dataterrawulf/pbs_singlechain_constant_init.sh \
	dataterrawulf/pbs_singlechain_cubedsphere_init.sh \
	dataterrawulf/pbs_singlechain_eastwest_init.sh \
	dataterrawulf/pbs_singlechain_northsouth_init.sh \
	dataterrawulf/pbs_singlechain_coredata_cont.sh \
	dataterrawulf/pbs_singlechain_constant_cont.sh \
	dataterrawulf/pbs_singlechain_cubedsphere_cont.sh \
	dataterrawulf/pbs_singlechain_eastwest_cont.sh \
	dataterrawulf/pbs_singlechain_northsouth_cont.sh \
	dataterrawulf/position_prior.txt \
	dataterrawulf/prior.txt \
	doc/instructions.tex \
	scripts/plot_likelihood_converge.py \
	scripts/plot_image_ortho.py

TARGETS = attenuationtomoS2Voronoi \
	attenuationtomoS2VoronoiPT \
	postS2Voronoi_mean \
	postS2Voronoi_mean_mpi \
	postS2Voronoi_likelihood \
	postS2Voronoi_text \
	mksynthetic \
	randommodelimage

all : $(TARGETS)

attenuationtomoS2Voronoi : attenuationtomoS2Voronoi.o $(OBJS)
	$(CXX) -o attenuationtomoS2Voronoi attenuationtomoS2Voronoi.o $(OBJS) $(LIBS)

attenuationtomoS2VoronoiPT : attenuationtomoS2VoronoiPT.o $(OBJS)
	$(CXX) -o attenuationtomoS2VoronoiPT attenuationtomoS2VoronoiPT.o $(OBJS) $(LIBS) $(MPI_LIBS)

postS2Voronoi_mean : postS2Voronoi_mean.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_mean.o $(OBJS) $(LIBS)

postS2Voronoi_mean_mpi : postS2Voronoi_mean_mpi.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_mean_mpi.o $(OBJS) $(LIBS)

postS2Voronoi_likelihood : postS2Voronoi_likelihood.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_likelihood.o $(OBJS) $(LIBS)

postS2Voronoi_text : postS2Voronoi_text.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_text.o $(OBJS) $(LIBS)

mksynthetic : mksynthetic.o $(OBJS)
	$(CXX) -o $@ mksynthetic.o $(OBJS) $(LIBS)

randommodelimage : randommodelimage.o $(OBJS)
	$(CXX) -o $@ randommodelimage.o $(OBJS) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = AttenuationVoronoi
TGZ = $(DIR).tar.gz

dist : all
	make -C data
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o




