CXXFLAGS = -std=c++11 -g -O0 -Wall -Wextra -pedantic
LIBS = -lsiggen -lpthread

nolib:
	$(CXX) -I .. $(CXXFLAGS) -c *.cpp
	$(CXX) -L .. -o main *.o $(LIBS)
	rm *.o
