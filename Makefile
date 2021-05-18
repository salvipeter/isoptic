all: isoptic

FLAGS=-Wall -std=c++17 -pedantic -O3
INCLUDES=-I../libgeom -I../dual-contouring
LIBS=-L../libgeom/release -lgeom -L../dual-contouring/build -ldualcontour

isoptic: isoptic.cc
	$(CXX) $(FLAGS) $(INCLUDES) -o $@ $< $(LIBS)
