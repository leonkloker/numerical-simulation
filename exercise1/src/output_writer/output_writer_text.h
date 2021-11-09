#pragma once

#include "output_writer/output_writer.h"

// Write *.txt files that are useful for debugging. All values are written
// to the file as they are stored in the field variables, no interpolation takes place. 
class OutputWriterText : OutputWriter
{
public:

    //! discretization shared pointer to the discretization object that will contain all the data to be written to the file 
    OutputWriterText(std::shared_ptr <Discretization> discretization);
    
    // Write current velocities to file, filename is output_<count>.txt 
    void writeFile (double currentTime) override;

    //! write only current values of pressure to file, filename is pressure_<count>.txt
    void writePressureFile ();

};
