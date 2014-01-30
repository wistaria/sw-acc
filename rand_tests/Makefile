CXX = pgc++
CXXFLAGS = -fast -acc -Minfo=accel -I/opt/cuda-5.5/include -L/opt/cuda-5.5/lib64
LIBS = -lcurand -lcudart

all : rand_pi rand_pi_curand_DevAPI rand_pi_curand_HostAPI

rand_pi : rand_pi.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

rand_pi_curand_DevAPI : rand_pi_curand_DevAPI.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

rand_pi_curand_HostAPI : rand_pi_curand_HostAPI.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)



clean:
	rm -f *.o sw_union_find sw_kalentev
