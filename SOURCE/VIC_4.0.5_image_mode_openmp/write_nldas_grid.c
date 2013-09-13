#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_data_nldas_grid.c,v 4.1.2.4 2004/05/06 20:26:41 tbohn Exp $";

void write_data_nldas_grid(out_data_struct *out_data,
        filenames_struct *filenames,
        global_param_struct *global,
        int              cell_cnt,
        int              nrecs,
		int              rec)
/**********************************************************************
  write_data_nldas_grid(), Ming Pan, May 2005, modified from write_data()

  This subroutine writes all energy and moisture balance parameters to
  a single gridded output file.

  OUTPUT:
	evaporation and vapor fluxes in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in cm
	snow depth in cm
	snow water equivlence in mm
	all energy fluxes are in W/m^2

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC
  3/98          Routine modified to output fluxes in PILPS2c 
                ASCII column format                             Dag
  4/30/98       Routine modified to add binary output options for
                improved file speed, and less disk usage for large
		model basins                                    KAC
  7/19/99       modified to output a single binary file containing
                the data selected for the LDAS project         KAC
  8/3/99        modified again to reduce the storage space needed
                for the LDAS output files.  
  1/4/2000      modified to allow both standard and LDAS formatted
                output using a compiler flag                    KAC
  3-12-03   added energy fluxes to snow band output files   KAC
  04-23-2003    modified LDAS SWQ output, so that it is multiplied by
                10 instead of 100 before being converted to a short
                integer.  This reduces stored value precision to 0.1,
                but increases the maximum storable SWQ, which was
                exceeded in previous LDAS simulations.          KAC

**********************************************************************/
{
  extern option_struct options;
  extern int              NCELLS;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  static int          nvars, ncells;
  static float       *tmp_grid;
  static float        lat[80000], lon[80000];
  
  FILE               *fp_soil;
  int                 i, j, k;
  int                 nldas_i, nldas_j, i_var, i_cell, i_rec;
  char                tmpstr[MAXSTRING];
  int                 dummy;
  float               *tmp_out;
  
  FILE               *fctl, *fout;
  char                fnctl[MAXSTRING], fnout[MAXSTRING];
  char                dmychar[MAXSTRING], monthchar[MAXSTRING], timechar[MAXSTRING];
  
  int                 tmpyear, tmpmonth, tmpday, tmphour, tmpjday=0;
  float               output_bytes;

  //fprintf(stderr, "Start writing cell %d, record %d\n", cell_cnt, rec);
  
  /* initialization */
  
  if (rec == -1) {
      
      /* Read soil param file to find out # of cells and lat-lons */
      
      fp_soil = open_file(filenames->soil, "r");
          /*vicerror("Can't open soil parameter file in write_data_nldas_grid().");*/
      
      fgets(tmpstr, MAXSTRING, fp_soil);
      
      ncells = 0;
      
      while (!feof(fp_soil)) {
          
          sscanf(tmpstr, "%d %d %f %f", &dummy, &dummy, &lat[ncells], &lon[ncells]);
          
          ncells++;
          
          fgets(tmpstr, MAXSTRING, fp_soil);
          
      }
      NCELLS=ncells;
      fclose(fp_soil);
      
      /* Determine # of variables to output */
      
      nvars  = 5;                                                  /* prec, evap, runoff, baseflow, Wdew */
      nvars += options.Nlayer;                                     /* moist */
      if (options.FULL_ENERGY || options.FROZEN_SOIL) nvars += 1;  /* rad_temp */
      nvars += 2;                                                  /* net_short, r_net */
      if (options.FULL_ENERGY || options.FROZEN_SOIL) nvars += 1;  /* latent */
      nvars += 5;                                                  /* evap_canop, evap_veg, evap_bare, sub_canop, sub_snow */
      if (options.FULL_ENERGY || options.FROZEN_SOIL) nvars += 2;  /* sensible, grnd_flux */
      nvars += 3;                                                  /* aero_resist, surf_temp, albedo */
      nvars += 3;                                                  /* swq, snow_depth, snow_canopy */
      if(options.FULL_ENERGY) nvars +=4;                           /* advection, deltaCC, snow_flux, refreeze_energy */
      
      /* Allocate memory */
      
      fprintf(stderr, "Output: Ncells = %d, Nvars = %d, Nrecs = %d\nTotal # of floats/bytes to be allocated: %d/%d\n",
              ncells, nvars, nrecs, ncells*nvars*nrecs, ncells*nvars*nrecs*sizeof(float));
      
      tmp_grid = (float *) calloc(ncells*nvars*nrecs, sizeof(float));
      if (tmp_grid == NULL)
          vicerror("Memory allocation error in write_data_nldas_grid().");
      /*
      for (i_rec=0; i_rec<nrecs; i_rec++)
          for (i_var=0; i_var<nvars; i_var++)
              for (i_cell=0; i_cell<ncells; i_cell++)
                  tmp_grid[nvars*ncells*i_rec+ncells*i_var+i_cell] = -9999;
      */
      return;
  }
  
  /* write to file */
  
  if (rec == -2) {
      
      /** Check output file size **/
      //output_bytes = options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY*nvars*nrecs*sizeof(float)*1.0;
      //fprintf(stderr, "Output file size %.0f bytes.\n", output_bytes);
      if (!options.OUTPUT_PER_STEP) {
          output_bytes = options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY*nvars*nrecs*sizeof(float)*1.0;
          //if (output_bytes >= 2.147483648e+9) {
          //    fprintf(stderr, "Output file size %.0f bytes exceeds 2GB limit for ext filesystem. Reset OUTPUT_PER_STEP to TRUE.\n", output_bytes);
          //    options.OUTPUT_PER_STEP = TRUE;
          //}
      }
      
      /* Allocate memory for tmp_out */
      
      tmp_out = (float *) calloc(options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY, sizeof(float));
      if (tmp_out == NULL)
          vicerror("Memory allocation error in write_data_nldas_grid(): tmp_out.");
      
      /* write output rec by rec, var by var */
      
      fprintf(stderr, "\nStart writing to file.\n");
      
      tmpyear  = global->startyear;
      tmpmonth = global->startmonth;
      tmpday   = global->startday;
      tmphour  = global->starthour;

      /* open output file */
      sprintf(dmychar, "%04i%02i%02i%02i", tmpyear, tmpmonth, tmpday, tmphour);
      strcpy(fnout, filenames->result_dir);
      strcat(fnout, "output_grid_");
      strcat(fnout, dmychar);
      fout = open_file(fnout, "wb");
      
      for (i_rec=0; i_rec<nrecs; i_rec++) {
                    
          for (i_var=0; i_var<nvars; i_var++) {
              
              /* initialize grid */
              for (i=0; i<options.OUTPUT_GRID_NY; i++)
                  for (j=0; j<options.OUTPUT_GRID_NX; j++) tmp_out[options.OUTPUT_GRID_NX*i+j] = UNDEF;
              
              for (i_cell=0; i_cell<ncells; i_cell++) {
                  
                  /* find out where the cell is */
                  nldas_i = (lon[i_cell]-options.OUTPUT_GRID_XO)/options.OUTPUT_GRID_DX + 0.5;
                  nldas_j = (lat[i_cell]-options.OUTPUT_GRID_YO)/options.OUTPUT_GRID_DY + 0.5;
                  
                  /*
                  if (i_rec==0) {
                      fprintf(stderr, "lat, lon: %f, %f       ldas i, j: %d, %d        value: %f\n",
                              lat[i_cell], lon[i_cell], nldas_i, nldas_j, tmp_grid[nvars*ncells*i_rec+ncells*i_var+i_cell]);
                  }
                  */
                      
                  tmp_out[options.OUTPUT_GRID_NX*nldas_j+nldas_i] = tmp_grid[nvars*ncells*i_rec+ncells*i_var+i_cell];
                  
              }
              
              fwrite(tmp_out, options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY, sizeof(float), fout);
              
          }
          
          if (options.OUTPUT_PER_STEP && i_rec<nrecs-1) {
              /* close output */
              fclose(fout);              
              /* march one day */
              get_next_time_step(&tmpyear, &tmpmonth, &tmpday, &tmphour, &tmpjday, global->dt);
              /* open file again */
              sprintf(dmychar, "%04i%02i%02i%02i", tmpyear, tmpmonth, tmpday, tmphour);
              strcpy(fnout, filenames->result_dir);
              strcat(fnout, "output_grid_");
              strcat(fnout, dmychar);
              fout = open_file(fnout, "wb");
          }
          
      }
      
      fclose(fout);              
      free(tmp_grid);
            
      /* try to create a GrADS control file for the output */
      
      sprintf(dmychar, "%04i%02i%02i%02i", global->startyear, global->startmonth, global->startday, global->starthour);
      strcpy(fnctl, filenames->result_dir);
      strcat(fnctl, "output_grid_");
      strcat(fnctl, dmychar);
      strcat(fnctl, ".ctl");

      fprintf(stderr, "The GrADS control file is %s.\n", fnctl);
      
      fctl = open_file(fnctl, "w");
      
      switch (global->startmonth) {
          case 1:  strcpy(monthchar, "JAN"); break;
          case 2:  strcpy(monthchar, "FEB"); break;
          case 3:  strcpy(monthchar, "MAR"); break;
          case 4:  strcpy(monthchar, "APR"); break;
          case 5:  strcpy(monthchar, "MAY"); break;
          case 6:  strcpy(monthchar, "JUN"); break;
          case 7:  strcpy(monthchar, "JUL"); break;
          case 8:  strcpy(monthchar, "AUG"); break;
          case 9:  strcpy(monthchar, "SEP"); break;
          case 10: strcpy(monthchar, "OCT"); break;
          case 11: strcpy(monthchar, "NOV"); break;
          case 12: strcpy(monthchar, "DEC"); break;
      }
      
      if (options.OUTPUT_PER_STEP) fprintf(fctl, "dset ^output_grid_%%y4%%m2%%d2%%h2\n");
      else fprintf(fctl, "dset ^output_grid_%s\n", dmychar);
      fprintf(fctl, "options template\n");
      fprintf(fctl, "title Gridded VIC Output\n");
      fprintf(fctl, "undef %f\n", UNDEF);
      fprintf(fctl, "xdef %d linear %lf %lf\n", options.OUTPUT_GRID_NX, options.OUTPUT_GRID_XO, options.OUTPUT_GRID_DX);
      fprintf(fctl, "ydef %d linear %lf %lf\n", options.OUTPUT_GRID_NY, options.OUTPUT_GRID_YO, options.OUTPUT_GRID_DY);
      fprintf(fctl, "zdef 1 linear 1 1\n");
      fprintf(fctl, "tdef %d linear %02i%s%04i %dhr\n", global->nrecs, global->startday, monthchar, global->startyear, global->dt);
      fprintf(fctl, "vars %d\n", nvars);
      fprintf(fctl, "prec            0  61 prec (mm)\n");
      fprintf(fctl, "evap            0  57 evap (mm)\n");
      fprintf(fctl, "runoff          0 235 runoff (mm)\n");
      fprintf(fctl, "baseflow        0 234 baseflow (mm)\n");
      fprintf(fctl, "Wdew            0  99 Wdew (mm)\n");
      for (i=1; i<=options.Nlayer; i++) fprintf(fctl, "sm%1d             0 151 soil moisture in layer %1d (mm)\n", i, i);
      if (options.FULL_ENERGY || options.FROZEN_SOIL) fprintf(fctl, "rad_temp        0 139 rad_temp (K)\n");
      fprintf(fctl, "net_short       0 111 net_short (W/m^2)\n");
      fprintf(fctl, "r_net           0  99 r_net (W/m^2)\n");
      if (options.FULL_ENERGY || options.FROZEN_SOIL) fprintf(fctl, "latent          0 121 latent (W/m^2)\n");
      fprintf(fctl, "evap_canop      0 200 evap_canop (mm)\n");
      fprintf(fctl, "evap_veg        0 210 evap_veg (mm)\n");
      fprintf(fctl, "evap_bare       0 199 evap_bare (mm)\n");
      fprintf(fctl, "sub_canop       0  99 sub_canop (mm)\n");
      fprintf(fctl, "sub_snow        0 173 sub_snow (mm)\n");
      if (options.FULL_ENERGY || options.FROZEN_SOIL) {
          fprintf(fctl, "sensible        0 122 sensible (W/m^2)\n");
          fprintf(fctl, "grnd_flux       0 155 grnd_flux (W/m^2)\n");
      }
      fprintf(fctl, "aero_resist     0 174 aero_resist (s/m)\n");
      fprintf(fctl, "surf_temp       0 138 surf_temp (C)\n");
      fprintf(fctl, "albedo          0  84 albedo (fraction)\n");
      fprintf(fctl, "swq             0  65 swq (mm)\n");
      fprintf(fctl, "snow_depth      0  66 snow_depth (cm)\n");
      fprintf(fctl, "snow_canop      0  99 snow_canop (mm)\n");
      if (options.FULL_ENERGY) {
          fprintf(fctl, "advection       0  99 advection (W/m^2)\n");
          fprintf(fctl, "deltaCC         0  99 deltaCC (W/m^2)\n");
          fprintf(fctl, "snow_flux       0 229 snow_flux (W/m^2)\n");
          fprintf(fctl, "refreeze_energy 0  99 refreeze_energy (W/m^2)\n");
      }
      fprintf(fctl, "endvars\n");

      fclose(fctl);
      
      return;
  }
  
  /************************************
    Attention: NO Frozen Soil or Snow Band Variables will be in the output
  ************************************/
  
  /************************************
    Output Standard Energy and Moisture Flux Variables
  ************************************/
  
  i_var = 0;
  
  //fprintf(stderr, "started writing at byte: %d ", nvars*ncells*rec+ncells*i_var+cell_cnt-1);
  
  /***** Write Binary Fluxes Variables *****/
  //tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->prec;                i_var++;
  //tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->evap;                i_var++;
  //tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->runoff;              i_var++;
  //tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->baseflow;            i_var++;
  //tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->Wdew;                i_var++;
  //for(j=0;j<options.Nlayer;j++) {    
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->moist[0];            i_var++;
  //}
  /*
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->rad_temp;            i_var++;
  }
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->net_short;           i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->r_net;               i_var++;
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->latent;              i_var++;
  }
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->evap_canop;          i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->evap_veg;            i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->evap_bare;           i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->sub_canop;           i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->sub_snow;            i_var++;
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->sensible;            i_var++;
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->grnd_flux;           i_var++;
  }
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->aero_resist;         i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->surf_temp;           i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->albedo;              i_var++;
  */  
  /***** Write Binary Snow Variables *****/
  /*
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->swq[0];              i_var++;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->snow_depth[0];       i_var++;
  
  //fprintf(stderr, "at byte: %d\n", nvars*ncells*rec+ncells*i_var+cell_cnt-1);
  
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->snow_canopy[0];      i_var++;
  if(options.FULL_ENERGY) {
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->advection[0];        i_var++;
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->deltaCC[0];          i_var++;
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->snow_flux[0];        i_var++;
    tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt-1] = (float)out_data->refreeze_energy[0];  i_var++;
  }
  */

  //fprintf(stderr, "finish writing cell %d, record %d\n", cell_cnt, rec);

}
