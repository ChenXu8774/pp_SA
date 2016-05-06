AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce

lib_LTLIBRARIES = libgprFit.la

AM_LDFLAGS = -L$(libdir) -L/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/lib/ -L$(MY_INSTALL)/lib `root-config --evelibs` -L${OFFLINE_MAIN}/lib

AM_CXXFLAGS = `root-config --cflags`

INCLUDES = -I/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/include -I$(OFFLINE_MAIN)/include/eigen3 \
					 -I$(MY_INSTALL)/include -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include 


##############################################
# please add new classes in alphabetical order

include_HEADERS = \
	gprFit.h

##############################################
# please add new classes in alphabetical order

libgprFit_la_SOURCES = \
  gprFitDict.C \
  gprFit.C

##############################################
# libraries used in the package

libgprFit_la_LIBADD = @ROOTGLIBS@ -lRooFit -lRooFitCore -lFoam -lMinuit -lMathMore -lGausProc -lhmanip

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libgprFit.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

gprFitDict.C: \
	gprFit.h \
  gprFitLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
