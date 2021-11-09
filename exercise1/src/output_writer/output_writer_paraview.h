#pragma once

#include "output_writer/output_writer.h"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

class OutputWriterParaview : OutputWriter
{
public:

    //! discretization shared pointer to the discretization object that will contain all the data to be written to the file 
    OutputWriterParaview(std::shared_ptr <Discretization> discretization);
    
    // Write current velocities to file, filename is output_<count>.txt 
    void writeFile (double currentTime) override;

private:
    //! write only current values of pressure to file, filename is pressure_<count>.txt
    vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;


};

