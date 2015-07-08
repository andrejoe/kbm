# Don't forget to make clean if switching between build systems. I may add a switch in the makefile
# if it is worth the time.

#ifeq ($(SYSTEM_NAME),Erwin)
#	PATH_PREFIX=/Users/erwin/Development/SyncOut/Magnus
#else
#	PATH_PREFIX=/home/andre
#endif

#ifeq ($(SYSTEM_NAME),Magnus)
#	PATH_PREFIX=/home/andre/
#endif

CC=g++ 
#CFLAGS= -std=c99
#CFLAGS= -std=gnu99
CFLAGS= -std=c++11

LIB_PATH=
INC_PATH=

all: kb
	./kb

kSystem.o: kSystem.cpp kSystem.h
	$(CC)$(CFLAGS) -c kSystem.cpp $(INC_PATH)

kCluster.o: kCluster.cpp kCluster.h
	$(CC)$(CFLAGS) -c kCluster.cpp $(INC_PATH)

kSphere.o: kSphere.cpp kSphere.h
	$(CC)$(CFLAGS) -c kSphere.cpp $(INC_PATH)

main.o: main.cpp
	$(CC)$(CFLAGS) -c main.cpp $(INC_PATH)

kb: main.o kSphere.o kCluster.o kSystem.o
	$(CC)$(CFLAGS) -o kb main.o kSphere.o kCluster.o kSystem.o

clean:
		rm -rf *.o
		rm -rf kb
		