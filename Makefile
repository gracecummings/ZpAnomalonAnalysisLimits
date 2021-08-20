# here we access the root configuration, include files, and libraries
ROOTCFLAGS=$(shell root-config --cflags)
ROOTINC=$(shell root-config --incdir)
ROOTLIBDIR=$(shell root-config --libdir)
ROOTLIBS=$(shell root-config --libs) -lMinuit -lMathMore
ROOTLDFLAGS=$(shell root-config --ldflags)

ROOTC=$(ROOTCFLAGS)
INC=../UHH2/JetMETObjects/interface
ROOTLINK=-L$(ROOTLIBDIR) $(ROOTLIBS) $(ROOTLDFLAGS)


# replace -W with -Wall to get more verbose warnings
TreeMakerTopiary.so: TreeMakerTopiary.C TreeMakerTopiary.h
	g++ -O -W $(ROOTC) -I$(INC) -fPIC -shared -oTreeMakerTopiary.so TreeMakerTopiary.C # build shared library


clean:
	rm -f TreeMakerTopiary.so
