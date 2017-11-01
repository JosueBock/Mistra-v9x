#
# Copyright 1996-2017 the Authors
#
# Licensed under the EUPL, Version 1.1 only (the "Licence");
#
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#   https://joinup.ec.europa.eu/software/page/eupl
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the Licence is distributed on an
# "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied.
#
# See the Licence for the specific language governing permissions
# and limitations under the Licence.

# Intel Fotran Compiler
#FC=ifort
#CFLAGS= -O2  -fpe0
#CFLAGS= -O2 -r8 -fpe0
#CFLAGS= -O2 -r8 -fpe0 -traceback
#CFLAGS= -check all -traceback -warn all
#CFLAGS= -O2 -r8 -fpe0 -fltconsistency
# CFLAGS= -O2 #-check all -warn declarations -warn uncalled -warn unused#-g #-warn all
#CFLAGS= -O2 -ipo #-check all -warn all
#CFLAGS= -fast # = -O3 -ipo -static (-xP)

# GNU Fortran
FC=gfortran
#  debug flags
CFLAGS= -fdefault-real-8 -fdefault-double-8 -fbacktrace -Wall -ffpe-trap=invalid,zero,overflow -fbounds-check -O0 -g
#CFLAGS= -fdefault-real-8 -fdefault-double-8 -fbacktrace -ffpe-trap=invalid,zero,overflow
#CFLAGS= -fdefault-real-8 -fdefault-double-8 -fbacktrace -fbounds-check -O0 -g
#CFLAGS= -fdefault-real-8 -fdefault-double-8 -fbacktrace -O3
#  run flags
#CFLAGS= -fdefault-real-8 -fdefault-double-8
#CFLAGS= -fdefault-real-8 -fdefault-double-8 -O3

# Makefile memo
# -------------
# $@ = target name
# $< = first dependency
# $^ = all dependencies


# Inclusion of Makefile in dependencies lists:
# --------------------------------------------
# when compilations option are changed in the Makefile, all files must be recompiled.
# To avoid mistakes in this situation, the Makefile has been added in the dependency list of all modules,
# and in the dependency list of all objects which don't already include a module

mistra: precision.o constants.o global_params.o common_modules.o directories.o \
	str.o radinit.o  nrad.o outp.o nuc.o jrate.o activity.o \
	tot.o bud_t.o bud_s_t.o \
	aer.o bud_a.o bud_s_a.o \
	gas.o bud_g.o bud_s_g.o \
	kpp.o utils.o \
	mod_out_netCDF.o out_netCDF.o
	$(FC) -o  $@ $(CFLAGS) $^ -lnetcdf -lnetcdff

# modules
global_params.o global_params.mod: global_params.f90 Makefile
	rm -f global_params.mod
	$(FC) $(CFLAGS) -c $<

common_modules.o gas_common.mod: common_modules.f global_params.mod Makefile
	rm -f gas_common.mod
	$(FC) $(CFLAGS) -c $<

constants.o constants.mod: constants.f90 Makefile
	rm -f constants.mod
	$(FC) $(CFLAGS) -c $<

directories.o directories.mod: directories.f Makefile
	rm -f directories.mod
	$(FC) $(CFLAGS) -c $<

mod_out_netCDF.o cdf_var_gas.mod: mod_out_netCDF.f90 Makefile
	rm -f cdf_var_gas.mod
	$(FC) $(CFLAGS) -c $<

precision.o precision.mod: precision.f90 Makefile
	rm -f precision.mod
	$(FC) $(CFLAGS) -c $<

# files containing modules
tot.o kpp_tot_parameters.mod kpp_tot_global.mod: tot.f tot_Parameters.h tot_Global.h tot_Sparse.h aer_mk.dat aer_km.dat \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_tot_parameters.mod kpp_tot_global.mod
	$(FC) $(CFLAGS) -c $<

aer.o kpp_aer_parameters.mod kpp_aer_global.mod: aer.f aer_Parameters.h aer_Global.h aer_Sparse.h aer_mk.dat aer_km.dat \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_aer_parameters.mod kpp_aer_global.mod
	$(FC) $(CFLAGS) -c $<

gas.o kpp_gas_parameters.mod kpp_gas_global.mod: gas.f gas_Parameters.h gas_Global.h gas_Sparse.h \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_gas_parameters.mod kpp_gas_global.mod
	$(FC) $(CFLAGS) -c $<

# other files
utils.o: utils.f gas_common.mod global_params.mod \
	kpp_aer_global.mod \
	kpp_gas_global.mod \
	kpp_tot_global.mod \
	Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

str.o: str.f constants.mod directories.mod gas_common.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

radinit.o: radinit.f constants.mod directories.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

nrad.o: nrad.f90 constants.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

outp.o: outp.f constants.mod global_params.mod gas_common.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

out_netCDF.o: out_netCDF.f global_params.mod gas_common.mod cdf_var_gas.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $< -I/usr/include

jrate.o: jrate.f directories.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

nuc.o: nuc.f constants.mod global_params.mod gas_common.mod Makefile
	rm -f mod_nuc.mod
	$(FC) $(CFLAGS) -o $(@F) -c $<

kpp.o: kpp.f tot_Parameters.h tot_Global.h aer_Parameters.h aer_Global.h gas_Parameters.h gas_Global.h \
	kpp_tot_global.mod kpp_tot_parameters.mod \
	kpp_aer_global.mod kpp_aer_parameters.mod \
	constants.mod gas_common.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

activity.o: activity.f90 constants.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_g.o: bud_g.f global_params.mod gas_Parameters.h gas_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_a.o: bud_a.f global_params.mod aer_Parameters.h aer_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_t.o: bud_t.f global_params.mod tot_Parameters.h tot_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_g.o: bud_s_g.f global_params.mod gas_Parameters.h gas_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_a.o: bud_s_a.f global_params.mod aer_Parameters.h aer_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_t.o: bud_s_t.f global_params.mod tot_Parameters.h tot_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<
