 # This file is auto generated by scons
g++ $1 -o $(basename "$1" .cpp) -O2 -I/data/snoplus/software/snocave_CentOS7/root-5.34.36/include -I/data/snoplus/software/snocave_CentOS7/geant4.10.0.p02/include -I/data/snoplus/software/snocave_CentOS7/gsl-1.16/include -I/home/lidgard/oxsx_sl7/include -L/data/snoplus/software/snocave_CentOS7/root-5.34.36/lib -L/data/snoplus/software/snocave_CentOS7/geant4.10.0.p02/lib64 -L/data/snoplus/software/snocave_CentOS7/gsl-1.16/lib -L/home/lidgard/oxsx_sl7/build -loxsx -lCore -lRIO -lHist -lGpad -lTree -lRint -lPostscript -lMatrix -lMathCore -lThread -lpthread -lm -ldl -lMinuit2 -lGraf -lPhysics -lCore -lRIO -lHist -lGpad -lTree -lRint -lPostscript -lMatrix -lMathCore -lThread -lpthread -lm -ldl -lMinuit2 -lGraf -lPhysics -lG4run -lG4clhep -lG4run -lG4clhep -larmadillo -larmadillo -lgsl -lgslcblas -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5  -Wl,-rpath,/data/snoplus/software/snocave_CentOS7/root-5.34.36/lib,-rpath,/data/snoplus/software/snocave_CentOS7/geant4.10.0.p02/lib64,-rpath,/data/snoplus/software/snocave_CentOS7/gsl-1.16/lib
