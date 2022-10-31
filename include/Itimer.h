//Itimer.h
//-------------------------------------------------------------------
#ifndef TIMER
#define TIMER
//--------------------------------------------------------------------
#include <omp.h>
//================================================
/**
 * @class Itimer
 * @author Igor Barden Grillo
 * @date 10/09/18
 * @file Itimer.h
 * @brief Class to storage and return the wall time of the program or its parts.
 */
class Itimer{
	public:
		//----------------------------------------------------------------
		/**
		 * @brief Double with the time in seconds when the time is beggining to count.
		 */
		double wall_init;
		
		//----------------------------------------------------------------
		/**
		 * @brief Double with the total program execution time in seconds.
		 */
		double tot_time;
		
		//----------------------------------------------------------------
		/**
		 * @brief Default constructor. 
		 */
		Itimer();
		//----------------------------------------------------------------
		/**
		 * @brief Copy constructor deleted
		 */
		Itimer(const Itimer& rhs_timer) = delete;
		
		//----------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading deleted
		 */
		Itimer& operator=(const Itimer& rhs_timer) = delete;
		
		//----------------------------------------------------------------
		/**
		 * @brief Count the wall time and return.
		 * @return Double with the wall time of exected program from the started time.
		 */
		double return_wall_time();
		
		//----------------------------------------------------------------
		/**
		 * @brief Reset the wall_init to start to count the time from the execution of this member function.
		 * @return None.
		 */
		void reset();
		
		//----------------------------------------------------------------
		/**
		 * @brief Return to the screen the wall time.
		 * @return None.
		 */
		void print();
		
		//----------------------------------------------------------------
		/**
		 * @brief Write the elapsed time on log.
		 */
		void write_in_log();
		
		//----------------------------------------------------------------
		/**
		 * @brief Controls when to throw the final message to the screen.
		 */
		void finish();
		
		//----------------------------------------------------------------
		/**
		 * @brief Destructor.
		 */
		~Itimer();
};	

#endif