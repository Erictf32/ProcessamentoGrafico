CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra
TARGET = raytracer
SOURCE = main.cpp

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

clean:
	rm -f $(TARGET) image.ppm image.png

render: $(TARGET)
	./$(TARGET)

convert: image.ppm
	convert image.ppm image.png

all: $(TARGET) render convert

.PHONY: clean render convert all