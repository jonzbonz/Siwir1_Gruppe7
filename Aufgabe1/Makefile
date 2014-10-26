CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -Ofast -DNDEBUG -mavx#-march=native

INCLUDES =
LDFLAGS =
LIBS =


#likwid
CXXFLAGS += -DUSE_LIKWID -pthread
INCLUDES += -I/usr/local/likwid-3.1.2/include/
LDFLAGS += -L/usr/local/likwid-3.1.2/lib/
LIBS += -llikwid

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

likwid: $(TARGET)
	rm -f likwid.log
	./likwidScript >> likwid.log

perfClean:
	rm -f perf.log

perfAll: perfClean perf32 perf64 perf128 perf256 perf512 perf1024 perf2048

perfAllGit: perfClean perf32 perf64 perf128 perf256 perf512 perf1024 perf2048
	 git commit -a -m "automatic commit of the perflog"
	 git push

perf32: $(TARGET)
	echo "Calculating 32x32" >> perf.log
	./$(TARGET) Matrices/perfMatrices/32x32-1 Matrices/perfMatrices/32x32-2 C.out >> perf.log

perf64: $(TARGET)
	echo "Calculating 64x64" >> perf.log
	./$(TARGET) Matrices/perfMatrices/64x64-1 Matrices/perfMatrices/64x64-2 C.out >> perf.log

perf128: $(TARGET)
	echo "Calculating 128x128" >> perf.log
	./$(TARGET) Matrices/perfMatrices/128x128-1 Matrices/perfMatrices/128x128-2 C.out >> perf.log

perf256: $(TARGET)
	echo "Calculating 256x256" >> perf.log
	./$(TARGET) Matrices/perfMatrices/256x256-1 Matrices/perfMatrices/256x256-2 C.out >> perf.log

perf512: $(TARGET)
	echo "Calculating 512x512" >> perf.log
	./$(TARGET) Matrices/perfMatrices/512x512-1 Matrices/perfMatrices/512x512-2 C.out >> perf.log

perf1024: $(TARGET)
	echo "Calculating 1024x1024" >> perf.log
	./$(TARGET) Matrices/perfMatrices/1024x1024-1 Matrices/perfMatrices/1024x1024-2 C.out >> perf.log

perf2048: $(TARGET)
	echo "Calculating 2048x2048" >> perf.log
	./$(TARGET) Matrices/perfMatrices/2048x2048-1 Matrices/perfMatrices/2048x2048-2 C.out >> perf.log

test1: $(TARGET)
	./$(TARGET) Matrices/testMatrices/A.in Matrices/testMatrices/B.in C.out
	diff -s C.out Matrices/testMatrices/C.out

test2: $(TARGET)
	./$(TARGET) Matrices/testMatrices/A2.in Matrices/testMatrices/B2.in C.out
	diff -s C.out Matrices/testMatrices/C2.out
	
test3: $(TARGET)
	./$(TARGET) Matrices/testMatrices/A3.in Matrices/testMatrices/B3.in C.out
	diff -s C.out Matrices/testMatrices/C3.out

testAll: test1 test2 test3
