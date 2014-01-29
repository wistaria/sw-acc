CXX = g++
CXXFLAGS = -O3 -I/opt/alps/boost_1_54_0
LIBS =

COMMON = observable.h options.h square_lattice.h

default: sw_union_find sw_kalentev

sw_union_find: sw_union_find.cxx union_find.h
	$(CXX) $(CXXFLAGS) -o sw_union_find sw_union_find.cxx $(LIBS)

sw_kalentev: sw_kalentev.cxx
	$(CXX) $(CXXFLAGS) -o sw_kalentev sw_kalentev.cxx $(LIBS)

clean:
	rm -f *.o sw_union_find sw_kalentev
