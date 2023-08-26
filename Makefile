CC = g++ -std=c++11
CFLAGS = -c -o

_alglibObjs = $(wildcard mining/alglib/*.cpp)
alglibObjs = $(patsubst %.cpp,%.o,$(_alglibObjs))

_evalObjs = $(wildcard eval/*.cpp)
evalObjs = $(patsubst %.cpp,%.o,$(_evalObjs))

all : main 

main : main.cpp ${evalObjs}
	$(CC) -o $@ $<

mining.a: mining/*.cpp alglib.a 
	ar rvs mining.a $<

alglib.a: $(alglibObjs)
	ar rvs alglib.a $(alglibObjs)

mining/alglib/%.o: mining/alglib/%.cpp     
	$(CC) $(CFLAGS) $@ $<

eval/%.o: eval/%.cpp     
	$(CC) $(CFLAGS) $@ $<


.PHONY: clean

clean:
	rm -f *.o