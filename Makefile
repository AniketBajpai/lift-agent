CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -w    disables warnings
CFLAGS  = -w -Ofast -march=native -std=c++11

# the build target executable:
TARGET = agent
TESTTARGET = agenttest
CPPFILES = Problem.cpp State.cpp UCTGraph.cpp

all: $(TARGET)

$(TARGET): $(CPPFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(CPPFILES)
	cp agent ../lift-sim/

run: $(TARGET)
	./$(TARGET)

test: $(CPPFILES)
	g++ -std=c++11 -fsanitize=address -o $(TESTTARGET) $(CPPFILES)
	./$(TESTTARGET)


clean:
	$(RM) $(TARGET)
