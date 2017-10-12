CXXFLAGS = -std=c++11 -g -O0 -Wall -Wextra -pedantic

default:
	$(CXX) -I. $(CXXFLAGS) -c *.cpp
	ar rcs libsiggen.a *.o
	rm -f *.o
	make nolib -C example/
	# make nolib -C Detectors/ICPC/
