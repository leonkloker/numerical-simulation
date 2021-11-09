#include "output_writer/output_writer_text.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

OutputWriterText::OutputWriterText(std::shared_ptr <Discretization> discretization) : 
OutputWriter(discretization_){}

void OutputWriterText::writeFile(double currentTime)
{
  /*
  // create "out" subdirectory if it does not yet exist
  int returnValue = system("mkdir -p out");
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl; 

  // Assemble the filename
  std::stringstream "velocities.txt ";
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo << ".txt"; //->GetDefaultFileExtension();

  std::cout << "Write file \"" << fileName.str() << "\"." << std::endl;
  

  // try to open file to read
  
  /*std::ofstream myfile ("output_velocties.txt");
  if (myfile.is_open())
  {
    myfile << currentTime;
    myfile << discretization_->u_;
    myfile << discretization_->v_;
    myfile.close();
  }
  else {
    std::cout << "Unable to open file";
  }
  */
  freopen( "output_velocties.txt", "a", stdout);
  std::cout<< currentTime;
  std::cout<< "\t u:";
  std::cout<< "\t";
  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++){
    for (int j = discretization_->uJBegin(); i <= discretization_->uJEnd(); j++){
      std::cout << discretization_->u(i,j); // goes to file.txt
      std::cout << "\t"; // goes to file.txt
      freopen("CON", "w", stdout);
    }
  }
    
  std::cout << "\n";
  std::cout<< "\t v:";
  std::cout<< "\t";
  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++){
    for (int j = discretization_->vJBegin(); i <= discretization_->vJEnd(); j++){
      std::cout << discretization_->v(i,j); // goes to file.txt
      std::cout << "tn"; // goes to file.txt
      freopen("CON", "w", stdout);
    }
  }
  std::cout << "\n";
}

void OutputWriterText::writePressureFile() 
{
  freopen( "output_velocties.txt", "a", stdout);
  for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++){
    for (int j = discretization_->pJBegin(); i <= discretization_->pJEnd(); j++){
      std::cout << discretization_->p(i,j); // goes to file.txt
      std::cout << "\t"; // goes to file.txt
      
    }
  }
  std::cout << "\n";
}