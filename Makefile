CC = g++
CFLAGS = -c -o

_alglibObjs = $(wildcard mining/alglib/*.cpp)
alglibObjs = $(patsubst %.cpp,%.o,$(_alglibObjs))

mining/alglib/%.o: mining/alglib/%.cpp     
	$(CC) $(CFLAGS) $@ $<

alglib.a: $(alglibObjs)
	ar rvs alglib.a $(alglibObjs)

.PHONY: clean

clean:
	rm -f *.o