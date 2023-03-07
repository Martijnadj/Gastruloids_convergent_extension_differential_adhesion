/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/


#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include "output.h"
#include "parse.h"

Parameter::Parameter() {

  T = 50.;
  target_area = 100;
  target_length = 60;
  lambda = 50;
  lambda2 = 5.0;
  Jtable = strdup("J.dat");
  conn_diss = 2000;
  vecadherinknockout = false;
  extensiononly = false;
  chemotaxis = 1000;
  border_energy = 100;
  neighbours = 2;
  periodic_boundaries = false;
  n_chem = 1;
  diff_coeff = new double[1];
  diff_coeff[0] = 1e-13;
  decay_rate = new double[1];
  decay_rate[0] = 1.8e-4;
  secr_rate = new double[1];
  secr_rate[0] = 1.8e-4;
  saturation = 0;
  dt = 2.0;
  dx = 2.0e-6;
  pde_its = 15;
  n_init_cells = 100;
  size_init_cells = 10;
  sizex = 200;
  sizey = 200;
  divisions = 0;
  mcs = 10000;
  rseed = -1;
  subfield = 1.0;
  relaxation = 0;
  storage_stride = 10;
	conv_ext = false;
  graphics = true;
  store = false;
  datadir = strdup("data_film");

	max_links = 3;
	theta_max = 45;
	r_max = 2.5 * target_length;
	lambda_force = 10;
	t_interval = 100;
	std_pol = 0;

	memory = 1;
	infl_redonred = 0;
	infl_redonyellow = 0;
	infl_yellowonred = 0;
	infl_yellowonyellow = 0;

	mean_divtime_epi = 100000;
	mean_divtime_meso = 80000;
	begin_diff = 500000;
	end_diff = 750000;
	mean_difftime = 50000;
	decrease_size_div = true;

	frac_yellow = 0;
	pulling_method = 2;


}

Parameter::~Parameter() {
  
  // destruct parameter object

  // free string parameter

  CleanUp();

}

void Parameter::CleanUp(void) {
  if (Jtable) 
     free(Jtable);
  if (diff_coeff) 
     free(diff_coeff);
  if (decay_rate) 
     free(decay_rate);
  if (secr_rate) 
     free(secr_rate);
  if (datadir) 
     free(datadir);

}

void Parameter::Read(const char *filename) {
  
  static bool ReadP=false;

  if (ReadP) {

    //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
    CleanUp();
    
  } else
    ReadP=true;

  FILE *fp=OpenReadFile(filename);


  T = fgetpar(fp, "T", 50., true);
  target_area = igetpar(fp, "target_area", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  Jtable = sgetpar(fp, "Jtable", "J.dat", true);
  conn_diss = igetpar(fp, "conn_diss", 2000, true);
  vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
  extensiononly = bgetpar(fp, "extensiononly", false, true);
  chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
  border_energy = igetpar(fp, "border_energy", 100, true);
  neighbours = igetpar(fp, "neighbours", 2, true);
  periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
  n_chem = igetpar(fp, "n_chem", 1, true);
  diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
  decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
  secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
  saturation = fgetpar(fp, "saturation", 0, true);
  dt = fgetpar(fp, "dt", 2.0, true);
  dx = fgetpar(fp, "dx", 2.0e-6, true);
  pde_its = igetpar(fp, "pde_its", 15, true);
  n_init_cells = igetpar(fp, "n_init_cells", 100, true);
  size_init_cells = igetpar(fp, "size_init_cells", 10, true);
  sizex = igetpar(fp, "sizex", 200, true);
  sizey = igetpar(fp, "sizey", 200, true);
  divisions = igetpar(fp, "divisions", 0, true);
  mcs = igetpar(fp, "mcs", 10000, true);
  rseed = igetpar(fp, "rseed", -1, true);
  subfield = fgetpar(fp, "subfield", 1.0, true);
  relaxation = igetpar(fp, "relaxation", 0, true);
  storage_stride = igetpar(fp, "storage_stride", 10, true);
	conv_ext = bgetpar(fp, "conv_ext", false, true);
  graphics = bgetpar(fp, "graphics", true, true);
  store = bgetpar(fp, "store", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);

	max_links = igetpar(fp, "max_links", 3, true);
	theta_max = igetpar(fp, "theta_max", 45, true);
	r_max = igetpar(fp, "r_max", 2.5, true);
	r_max = r_max * target_length;
	lambda_force = fgetpar(fp, "lambda_force", 10, true);
	t_interval = igetpar(fp, "t_interval", 100, true);
	std_pol = igetpar(fp, "std_pol", 0, true);

	memory = fgetpar(fp, "memory", 1, true);
	infl_redonred = fgetpar(fp, "infl_redonred", 0, true);
	infl_redonyellow = fgetpar(fp, "infl_redonyellow", 0, true);
	infl_yellowonred  = fgetpar(fp, "infl_yellowonred", 0, true);
	infl_yellowonyellow = fgetpar(fp, "infl_yellowonyellow", 0, true);

	mean_divtime_epi = igetpar(fp, "mean_divtime_epi", 100000, true);
	mean_divtime_meso = igetpar(fp, "mean_divtime_meso", 80000, true);
	begin_diff = igetpar(fp, "begin_diff", 500000, true);
	end_diff = igetpar(fp, "end_diff", 750000, true);
	mean_difftime = igetpar(fp, "mean_difftime", 50000, true);
	decrease_size_div = bgetpar(fp, "decrease_size_div", true, true);

	frac_yellow = fgetpar(fp, "frac_yellow", 0, true);

	pulling_method = igetpar(fp, "pulling_method", 2, true);
}

const char *sbool(const bool &p) {

  const char *true_str="true";
  const char *false_str="false";
  if (p)
    return true_str;
  else
    return false_str;
}

void Parameter::Write(ostream &os) const {
  setlocale(LC_NUMERIC, "C");

  os << " T = " << T << endl;
  os << " target_area = " << target_area << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " lambda2 = " << lambda2 << endl;

  if (Jtable) 
    os << " Jtable = " << Jtable << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
  os << " extensiononly = " << sbool(extensiononly) << endl;
  os << " chemotaxis = " << chemotaxis << endl;
  os << " border_energy = " << border_energy << endl;
  os << " neighbours = " << neighbours << endl;
  os << " periodic_boundaries = " << sbool(periodic_boundaries) << endl;
  os << " n_chem = " << n_chem << endl;
  os << " diff_coeff = "<< diff_coeff[0] << endl;
  os << " decay_rate = "<< decay_rate[0] << endl;
  os << " secr_rate = "<< secr_rate[0] << endl;
  os << " saturation = " << saturation << endl;
  os << " dt = " << dt << endl;
  os << " dx = " << dx << endl;
  os << " pde_its = " << pde_its << endl;
  os << " n_init_cells = " << n_init_cells << endl;
  os << " size_init_cells = " << size_init_cells << endl;
  os << " sizex = " << sizex << endl;
  os << " sizey = " << sizey << endl;
  os << " divisions = " << divisions << endl;
  os << " mcs = " << mcs << endl;
  os << " rseed = " << rseed << endl;
  os << " subfield = " << subfield << endl;
  os << " relaxation = " << relaxation << endl;
  os << " storage_stride = " << storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;

	os << " max_links = " << max_links << endl;
  os << " theta_max = " << theta_max << endl;
  os << " r_max = " << r_max << endl;
  os << " lambda_force = " << lambda_force << endl;
  os << " t_interval = " << t_interval << endl;
  os << " std_pol = " << std_pol << endl;

	os << " memory = " << memory << endl;
	os << " infl_redonred = " << infl_redonred << endl;
	os << " infl_redonyellow = " << infl_redonyellow << endl;
	os << " infl_yellowonred = " << infl_yellowonred << endl;
	os << " infl_yellowonyellow = " << infl_yellowonyellow << endl;

  os << " mean_divtime_epi = " << mean_divtime_epi << endl;
  os << " mean_divtime_meso = " << mean_divtime_meso << endl;
  os << " begin_diff = " << begin_diff << endl;
  os << " end_diff = " << end_diff << endl;
  os << " mean_difftime = " << mean_difftime << endl;
	os << " decrease_size_div = " << decrease_size_div << endl;

	os << " frac_yellow = " << frac_yellow << endl;
	os << " pulling_method = " << frac_yellow << endl;
	

  if (datadir) 
    os << " datadir = " << datadir << endl;
}


ostream &operator<<(ostream &os, Parameter &p) {
  p.Write(os);
  return os;
}

Parameter par;
