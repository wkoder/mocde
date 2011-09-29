CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		mocde.o

LIBS =

TARGET =	mocde

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
