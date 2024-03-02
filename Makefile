CC = g++ -std=c++17
CFLAGS = -c -o

_alglibObjs = $(wildcard mining/alglib/*.cpp)
alglibObjs = $(patsubst %.cpp,%.o,$(_alglibObjs))

_miningObjs = $(wildcard mining/*.cpp)
miningObjs = $(patsubst %.cpp,%.o,$(_miningObjs))

_evalObjs = $(wildcard eval/*.cpp)
evalObjs = $(patsubst %.cpp,%.o,$(_evalObjs))

_heuristicsObjs = $(wildcard heuristics/*.cpp)
heuristicsObjs = $(patsubst %.cpp,%.o,$(_heuristicsObjs))

_utilsObjs = $(wildcard utils/*.cpp)
utilsObjs = $(patsubst %.cpp,%.o,$(_utilsObjs))

all : solver 

solver : main.cpp eval.a heuristics.a ${utilsObjs}
	$(CC) -o $@ main.cpp ${evalObjs} ${utilsObjs} heuristics.a

heuristics.a: ${heuristicsObjs} ${miningObjs} ${alglibObjs} 
	ar rvs heuristics.a ${heuristicsObjs} ${miningObjs} ${alglibObjs} 

# heuristics.a: $(heuristicsObjs) alglib.a mining.a
#	ar rvs heuristics.a $<

mining.a: ${miningObjs} ${alglibObjs} 
	ar rvs mining.a ${miningObjs} ${alglibObjs} 

eval.a: ${evalObjs}
	ar rvs eval.a ${evalObjs}

alglib.a: ${alglibObjs}
	ar rvs alglib.a ${alglibObjs}

mining/%.o: mining/%.cpp     
	$(CC) $(CFLAGS) $@ $<

mining/alglib/%.o: mining/alglib/%.cpp     
	$(CC) $(CFLAGS) $@ $<

eval/%.o: eval/%.cpp     
	$(CC) $(CFLAGS) $@ $<

heuristics/%.o: heuristics/%.cpp     
	$(CC) $(CFLAGS) $@ $<

utils/%.o: utils/%.cpp     
	$(CC) $(CFLAGS) $@ $<

#heuristics/%.o: heuristics/%.cpp mining.a stdlib.a 
#	$(CC)  $@ $< 

# heuristics/dm_lshade.cpp: heuristics/dm_lshade.cpp alglib.a
# 	g++ -std=c++11 -o dm_lshade.o heuristics/dm_lshade.cpp alglib.a

clean:
	@find ./ -iname "solver" -exec rm {} \;
	@find ./ -iname "*.o" -exec rm {} \;
	@find ./ -iname "*.a" -exec rm {} \;