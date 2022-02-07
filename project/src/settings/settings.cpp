#include "settings.h"

void Settings::loadFromFile(std::string filename)
{
   std::fstream newfile;
   newfile.open(filename,std::ios::in); //open a file to perform read operation using file object
   
   if (newfile.is_open()){   //checking whether the file is open
      std::string parameter;
      
      while(getline(newfile, parameter))
           { //read data from file object and put it into string.
	   // remove whitespace at beginning of line (if there is any)
	   // if first character is a '#', skip line (line[0] == '#')
	   // if line does not contain a '=' sign, skip line
	   // parse parameter name
	   // remove trailing spaces from parameterName
      // parse value
		// remove whitespace at beginning of value
	   // remove comments at end of value
  	   // remove whitespace at end of value
	   if (parameter.find_first_of(" \t") != std::string::npos)
	      {
	      if (!parameter.find_first_not_of('#'))
	         {
	         int parameterNameIndex = parameter.find_first_of(" \t");
		      std::string parameterName = parameter;
		      parameterName = parameterName.erase(parameterNameIndex);
		      std::string parameterValue = parameter.erase(0,parameter.find_first_of("="));
	         parameterValue = parameter.erase(0,parameter.find_first_of(" \t")+1);
	         
            if (parameter.find_first_of(" \t") != std::string::npos)
		      {
  	            parameterValue = parameter.erase(parameter.find_first_of(" "));
	         }
	       
		 
	         // parse actual value and set corresponding parameter
	  	 
	         if (parameterName == "physicalSizeX")
	  	    {
		     physicalSize[0] = atof(parameterValue.c_str());
		    }
                 else if (parameterName == "physicalSizeY")
		    {
		     physicalSize[1] = atof(parameterValue.c_str());
		    }
		 else if (parameterName == "endTime")
                    {
                     endTime = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "re")
                    {
                     re = atoi(parameterValue.c_str());
                    }
		 else if (parameterName == "pr")
		    {
		     pr = atoi(parameterValue.c_str());
		    }
		 else if (parameterName == "gX")
                    {
                     g[0] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "gY")
                    {
                     g[1] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletBottomX")
                    {
                     dirichletBcBottom[0] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletBottomY")
                    {
                     dirichletBcBottom[1] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletTopX")
                    {
                     dirichletBcTop[0] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletTopY")
                    {
                     dirichletBcTop[1] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletLeftX")
                    {
                     dirichletBcLeft[0] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletLeftY")
                    {
                     dirichletBcLeft[1] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletRightX")
                    {
                     dirichletBcRight[0] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "dirichletRightY")
                    {
                     dirichletBcRight[1] = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "nCellsX")
                    {
                     nCells[0] = atoi(parameterValue.c_str());
                    }
		 else if (parameterName == "nCellsY")
                    {
                     nCells[1] = atoi(parameterValue.c_str());
                    }
		 else if (parameterName == "useDonorCell")
                    {
                     if (parameterValue == "True" || parameterValue == "true"){
                        useDonorCell = true;
                     }else{
                        useDonorCell = false;
                     }
                    }
		 else if (parameterName == "alpha")
                    {
                     alpha = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "tau")
                    {
                     tau = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "maximumDt")
                    {
                     maximumDt = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "pressureSolver")
                    {
                     pressureSolver = parameterValue.c_str();
                    }
		 else if (parameterName == "omega")
                    {
                     omega = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "epsilon")
                    {
                     epsilon = atof(parameterValue.c_str());
                    }
		 else if (parameterName == "maximumNumberOfIterations")
                    {
                     maximumNumberOfIterations = atof(parameterValue.c_str());
                    }
	         
           	  }
	         }
           }

      newfile.close(); //close the file object.
   }
}

void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << "s, Reynolds number: " << re << ", external forces: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}
