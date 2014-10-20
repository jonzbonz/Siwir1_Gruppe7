CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG

INCLUDES =
LDFLAGS =
LIBS =


# likwid
#CXXFLAGS += -DUSE_LIKWID -pthread
#INCLUDES += -I/usr/local/likwid-3.1.2/include/
#LDFLAGS += -L/usr/local/likwid-3.1.2/lib/
#LIBS += -llikwid

TARGET = matmult
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp Timer.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)

run: $(TARGET)
	./$(TARGET) A.in B.in C.out

perfAll: perf32 perf64 perf128 perf256 perf512 perf1024 perf2048

perf32: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/32x32-1 Matrices/perfMatrices/32x32-2 C.out

perf64: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/64x64-1 Matrices/perfMatrices/64x64-2 C.out

perf128: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/128x128-1 Matrices/perfMatrices/128x128-2 C.out

perf256: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/256x256-1 Matrices/perfMatrices/256x256-2 C.out

perf512: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/512x512-1 Matrices/perfMatrices/512x512-2 C.out

perf1024: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/1024x1024-1 Matrices/perfMatrices/1024x1024-2 C.out

perf2048: $(TARGET)
	./$(TARGET) Matrices/perfMatrices/2048x2048-1 Matrices/perfMatrices/2048x2048-2 C.out
