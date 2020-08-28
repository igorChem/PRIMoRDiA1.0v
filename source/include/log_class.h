//log.h
//------------------------------

/*********************************************************************/
/* This source code file is part of PRIMoRDiA software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
 
#ifndef LOGCLASS
#define LOGCLASS
//------------------------------
#include <string>
#include <fstream>
//=================================================
/**
 * @class Ilog
 * @author Igor Barden Grillo
 * @date 10/09/18
 * @file log_class.h
 * @brief Class with tools to output log of the program.
 * 
 * Class to store the messages for each library class to output 
 * in the scrren or/and in a file.
 */
class Ilog{
	public:
		std::ofstream log_file;
		bool screen_output;
		Ilog();
		Ilog(const Ilog& rhs_log) = delete;
		Ilog& operator=(const Ilog& rhs_log) = delete;
		//-----------------------------------------------------------------------
		/**
		 * @brief Initializes the object information.
		 * @param String with the name of the file to being writting.
		 * @param If the messages will be printed on the screen.
		 * @return None.
		 */
		void initialize(bool sout);
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(std::string message);
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(int message );
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(double message);
		//------------------------------------------------------------------------
		/**
		 * @brief Get the time count and returns to the console or/and to the log file.
		 * @return None.
		 */
		void timer();
		//------------------------------------------------------------------------
		/**
		 * @brief Exit the program execution printing a message to the console 
		 * @param message
		 * @return None.
		 */
		void abort(std::string message);
		//-----------------------------------------------------------------------
		/**
		 * @brief Default destructor. Here the text is written in the output stream
		 * when the object goes out of scope.
		 */
		~Ilog();
		
};

#endif