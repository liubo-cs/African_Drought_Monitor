#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id: read_once_ldas_forcing.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

float **read_once_nldas_forcing(filenames_struct     *filenames,
                                global_param_struct   global_param)
/**********************************************************************
  read_once_nldas_forcing()   Ming Pan, mpan@princeton.edu  May 2005

  This subroutine reads all grid forcing data in one shot.

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NR, NF;
  extern int              NCELLS, flag;
  
  FILE         *fin, *fp_soil;
  char         fnin[MAXSTRING], dmychar[MAXSTRING], tmpstr[MAXSTRING];
  int          dummy, ncells, nvars;
  int          nldas_i, nldas_j, i_cell, i_rec;
  float        lat[80000], lon[80000];
  float        *tmp_data;

  char                 errorstr[MAXSTRING];
  int                  i;
  long                 nread;
  float             **forcing_data;

  /* Allocate memory for tmp_grid */
  
  tmp_data = (float *) calloc(options.INPUT_GRID_NX*options.INPUT_GRID_NY, sizeof(float));
  if (tmp_data == NULL)
          vicerror("Memory allocation error in read_once_nldas_forcing(): tmp_data");
  
  /* Read soil param file to find out # of cells and lat-lons */
      
  fp_soil = open_file(filenames->soil, "r");
      /*vicerror("Can't open soil parameter file in write_data_nldas_grid().");*/
  
  fgets(tmpstr, MAXSTRING, fp_soil);
  
  ncells = 0;
  
  while (!feof(fp_soil)) {
      
      sscanf(tmpstr, "%d %d %f %f", &flag, &dummy, &lat[ncells], &lon[ncells]);
      
      if (flag) ncells++;
      
      fgets(tmpstr, MAXSTRING, fp_soil);
      
  }
  
  fclose(fp_soil);
  
  nvars = 0;
  
  /** Allocate data arrays for input forcing data **/
  forcing_data = (float **)calloc(N_FORCING_TYPES,sizeof(float*));
  for(i=0;i<N_FORCING_TYPES;i++) {
      if (param_set.TYPE[i].SUPPLIED) {
          forcing_data[i] = (float *)calloc((global_param.nrecs * ncells), sizeof(float));
          if (forcing_data[i] == NULL)
              vicerror("Memory allocation error in read_once_nldas_forcing(): forcing_data");
          nvars++;
      }
  }
  
  fprintf(stderr, "Forcing: Ncells = %d, Nvars = %d, Nrecs = %d\nTotal # of floats/bytes allocated: %d/%d\n",
              ncells, nvars, global_param.nrecs, ncells*nvars*global_param.nrecs, ncells*nvars*global_param.nrecs*sizeof(float));

  /* open input forcing file */
  
  sprintf(dmychar, "%04i%02i%02i", global_param.forceyear[0], global_param.forcemonth[0], global_param.forceday[0]);
  strcpy(fnin, filenames->forcing[0]);
  strcat(fnin, dmychar);
  fin = open_file(fnin, "rb");
  
  /* read from file */
  
  fseek(fin, options.INPUT_GRID_NX*options.INPUT_GRID_NY*nvars*global_param.forceskip[0]*sizeof(float), SEEK_SET);
  
  for (i_rec=0; i_rec<global_param.nrecs; i_rec++) {
      
      for (i=0;i<N_FORCING_TYPES;i++) {
          
          if (param_set.TYPE[i].SUPPLIED) {
          
              nread = fread(tmp_data, sizeof(float), options.INPUT_GRID_NX*options.INPUT_GRID_NY, fin);
              
              if (nread != options.INPUT_GRID_NX*options.INPUT_GRID_NY) {
                  sprintf(errorstr, "Forcing file reading error: %ld floats read in rec %d. read_once_nldas_forcing()", nread, i_rec);
                  vicerror(errorstr);
              }
              
              for (i_cell=0; i_cell<ncells; i_cell++) {
                  
                  nldas_i = (lon[i_cell]-options.INPUT_GRID_XO)/options.INPUT_GRID_DX + 0.5;
                  nldas_j = (lat[i_cell]-options.INPUT_GRID_YO)/options.INPUT_GRID_DY + 0.5;
                  
                  forcing_data[i][global_param.nrecs*i_cell+i_rec] = tmp_data[options.INPUT_GRID_NX*nldas_j+nldas_i];
                  
              }
              
          }
          
      } /* end of var loop */
      
  } /* end of rec loop */
  
  fclose(fin);
  
  NCELLS = ncells;     

  return(forcing_data);

}
