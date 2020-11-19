main: frenet_optimal_trajectory.o quintic_polynomials_planner.o cubic_spline_planner.o
			g++ -o main frenet_optimal_trajectory.o quintic_polynomials_planner.o cubic_spline_planner.o -I/usr/include/python2.7 -lpython2.7 -g -std=c++14
frenet_optimal_trajectory.o: frenet_optimal_trajectory.cpp
			g++ -c frenet_optimal_trajectory.cpp -I/usr/include/python2.7 -lpython2.7 -std=c++14
quintic_polynomials_planner.o: quintic_polynomials_planner.cpp
			g++ -c quintic_polynomials_planner.cpp -std=c++14
cubic_spline_planner.o: cubic_spline_planner.cpp 
			g++ -c cubic_spline_planner.cpp -std=c++14
