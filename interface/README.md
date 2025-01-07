    export LD_LIBRARY_PATH=/home/dndus0107/CMSSW_14_0_19_patch2/src/interface/:$LD_LIBRARY_PATH

cd interface
rootcling -f AnalysisClassesDict.cxx -c AnalysisClasses.h AnalysisClassesLinkDef.h
g++ -shared -fPIC -o libAnalysisClasses.so AnalysisClassesDict.cxx \
    -I$(root-config --incdir) \
    -I$CMSSW_RELEASE_BASE/src \
    -I$CMSSW_BASE/src \
    $(root-config --libs)

g++ -o muonCorrection_tmp muonCorrection_tmp.cc \
$(root-config --cflags --libs) \
-I$CMSSW_BASE/src \
-I$CMSSW_RELEASE_BASE/src \
-L$CMSSW_BASE/lib \
-L$CMSSW_RELEASE_BASE/lib \
-lAnalysisClasses \
