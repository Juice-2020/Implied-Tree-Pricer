oop2: 	main.o Market.o Date.o Barrier.o
	g++  -o oop.out main.o Market.o Date.o Barrier.o
# the mainline of our example
main.o: main.cpp TreePricer.h TreeProduct.h EuropeanTrade.h AmericanTrade.h BarrierTrade.h Barrier.h
	g++ -c main.cpp
Market.o: Market.cpp Market.h Date.h
	g++ -c Market.cpp 
Date.o: Date.cpp Date.h
	g++ -c Date.cpp
Barrier.o : Barrier.h Barrier.cpp
	g++ -c Barrier.cpp
# EuropeanTrade.o: EuropeanTrade.cpp EuropeanTrade.h Types.h
#	g++ -c EuropeanTrade.cpp
#AmericanTrade.o: AmericanTrade.cpp AmericanTrade.h Types.h
#	g++ -c AmericanTrade.cpp
# clean all output files
clean:
	rm *.o *.out
