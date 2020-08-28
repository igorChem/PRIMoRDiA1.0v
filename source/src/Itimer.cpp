//Itimer.cpp
//==================
#include <iostream>
#include <omp.h>
//==================
#include "../include/common.h"
#include "../include/Itimer.h"
/***********************************************************************************/
Itimer::Itimer(){ tot_time = wall_init = omp_get_wtime(); }
/***********************************************************************************/
double Itimer::return_wall_time(){ return omp_get_wtime() - wall_init; }
/***********************************************************************************/
void Itimer::reset(){ wall_init = omp_get_wtime(); }
/************************************************************************************/
void Itimer::print(){ std::cout << ( omp_get_wtime() - wall_init ) <<  std::endl; }
/************************************************************************************/
Itimer::~Itimer(){
	tot_time = omp_get_wtime() - tot_time;
	std::cout << "Total execution time of PRIMoRDiA program: " << tot_time << " seconds" << std::endl; 
}
/************************************************************************************/