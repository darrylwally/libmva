/*! \file mvfile.h
  \brief Contains methods for storing and retrieving mvmodels
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2014

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#ifndef MVFILE_H
#define MVFILE_H

#include "mvmodel.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! Writes an array of MVModels to disk
  \arg filename the desired output filename
  \arg models pointer to array of MVModel structs
  \arg num_models number of models.
  \return MVReturnCode (SUCCESS or FILE_WRITE_ERROR)
  */
MVReturnCode mvfile_write(const char * filename, MVModel *models, size_t num_models);

/*! Read an array of MVModels from disk
  \arg models reference to pointer of array of models that will be allocated.
  \arg filename the file to read.
  \return number of models that have been returned or -1 if error.
  */
int mvfile_read(MVModel **models, const char *filename);

#ifdef __cplusplus
}
#endif

#endif // MVFILE_H
