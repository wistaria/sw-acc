CXX = g++
CXXFLAGS = -O3 -I/opt/alps/boost_1_54_0
LIBS =

COMMON = observable.h options.h square_lattice.h

default: sw_union_find

sw_union_find : sw_union_find.C union_find.h
	$(CXX) $(CXXFLAGS) -o sw_union_find sw_union_find.C $(LIBS)
