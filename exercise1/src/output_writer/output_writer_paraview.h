#pragma once
#include <output_writer.h>
/*
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
*/
//! Write a file out/output_<fileNo>.vti to be visualized in ParaView.
//! It contains 10x10 nodes with an artifical pressure field.
//! This method is only for demonstration purpose and does nothing useful.
//! However, we will provide similar files, e.g. "output_writer_paraview.h", to be used in the submission code.
void writeParaviewOutput(int fileNo);


class OutputWriterParaview : OutputWriter
{
public:

    //! discretization shared pointer to the discretization object that will contain all the data to be written to the file 
    OutputWriterParaview(std::shared_ptr <Discretization> discretization);
    
    // Write current velocities to file, filename is output_<count>.txt 
    void OutputWriterParaview::writeFile (double currentTime);

private:
    //! write only current values of pressure to file, filename is pressure_<count>.txt
    vtkSmartPointer<vtkXMLImageDataWriter> OutputWriterParaview::vtkWriter_;


};
