CXX = pgc++
CXXFLAGS = -fast -acc -Minfo=accel -I/opt/cuda-5.5/include -I$(HOME)/opt/include -L/opt/cuda-5.5/lib64
LIBS = -lcurand -lcudart

all : rand_pi rand_pi_curand_DevAPI rand_pi_curand_HostAPI

rand_pi : rand_pi.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

rand_pi_curand_DevAPI : rand_pi_curand_DevAPI.cxx
	$(CXX) $(CXXFLAGS) -Minline=levels:3 -o $@ $< $(LIBS)

rand_pi_curand_HostAPI : rand_pi_curand_HostAPI.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

rand_pi_omp : rand_pi_omp.cxx
	$(CXX) -openmp -I$(HOME)/opt/include -o $@ $<

