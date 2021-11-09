#include "output_writer/output_writer.h"

OutputWriter::OutputWriter(std::shared_ptr <Discretization> discretization){
    discretization_ = discretization;
    fileNo_ = 1;
}

