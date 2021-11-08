#pragma once

#include "discretization.h"
#include <memory>

// Inteface class for writing simulation data output. 
class OutputWriter
{
public:

    //! constructor: discretization shared pointer to the discretization object that will 
    // contain all the data to be written to the file 
    OutputWriter(std::shared_ptr <Discretization> discretization);

    //! write current velocities to file, filename is output_<count>.vti
    virtual void OutputWriter::writeFile (double currentTime) = 0;

protected:
    
    std::shared_ptr<Discretization> OutputWriter::discretization_;
    int OutputWriter::fileNo_;

};
