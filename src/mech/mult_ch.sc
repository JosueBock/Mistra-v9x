##/bin/csh
# called by make_kpp.sc
# change *.f *.h for each input
#
# created by Roland von Glasow
# updated in June 2016 by Josue Bock when upgrading KPP

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


# define variables
# ----------------
set file_f=$1.f
set file_h=("$1_Global.h" "$1_Parameters.h")
set file_s=$1_Sparse.h

set appendix=$2

# ==============
# work on file_f
# ==============

# delete unnecessary subroutines
# ------------------------------
set del_SR=("Update_SUN" "Update_PHOTO" "KppDecompCmplx" "KppSolveIndirect" "KppSolveCmplx" "KppSolveTR" "WCOPY" "WSCAL" "WLAMCH" "SET2ZERO" "WDOT" "WADD" "Jac_SP_Vec" "JacTR_SP_Vec")
set fold=$file_f
set inn=0
foreach i ($del_SR)
  @ inn=$inn + 1
  set fnew=tmp$inn
# [CEN] accomodates either a comment, the end of word "SUBROUTINE" or the end of word "FUNCTION"
# the special quoting ('"$i"') allows use of shell variables in sed
  sed '/[CEN] '"$i"'/,/End of '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fold $file_f

# delete calls to Update_SUN and Update_PHOTO
# -------------------------------------------
set del_CALL=("Update_SUN" "Update_PHOTO")
set fold=$file_f
set inn=0
foreach i ($del_CALL)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed '/CALL '"$i"'/,/)/d' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fold $file_f

# rename remaining subroutines
# ----------------------------
set ch_SR=("Update_RCONST" "INTEGRATE" "RosenbrockIntegrator" "ros_ErrorNorm" "ros_FunTimeDerivative" "ros_PrepareMatrix" "ros_ErrorMsg" "Ros2" "Ros3" "Ros4" "Rodas3" "Rodas4" "DecompTemplate" "SolveTemplate" "FunTemplate" "JacTemplate" "KppDecomp" "WAXPY" "KppSolve" "WLAMCH_ADD" "Jac_SP")
set fold=$file_f
set inn=0
foreach i ($ch_SR)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end

# Special cases:
# 1) for SR Fun:
#   needs to avoid replace every occurrence of Fun*, which is a substring of "FunTemplate" and is obviously also a substring of "Function"
#   the search pattern accomodates "Fun(", "Fun (", "Fun -" and "Fun f"
@ inn=$inn + 1
set fnew=tmp$inn
# " \?" matches zero or one space (backslash is needed in pattern context)
# [f(-] matches one of the 3 character inside the brackets
# \( ... \) saves the pattern in between, which is recalled with \1
sed 's/Fun\( \?[f(-]\)/Fun_'"$appendix"'\1/g' $fold > ! $fnew
set fold=$fnew

# 2) for SR Rosenbrock : add ( to the searched pattern
@ inn=$inn + 1
set fnew=tmp$inn
sed 's/Rosenbrock(/Rosenbrock_'"$appendix"'(/g' $fold > ! $fnew
mv -f $fnew $file_f

# delete unnecessary LU_IROW data
# -------------------------------
# explanations, if future change is needed:
#   \+ means one or more
#   \ {1\} means exactly one
# searched ending string:
#   - starts by "*"
#   - optional pattern: optional space(s) followed by at least 1 number followed by comma
#   - mandatory pattern: optional space(s) followed by at least 1 number followed by 1 space and /
set fold=$file_f
@ inn=$inn + 1
set fnew=tmp$inn
sed '/DATA( LU_IROW(i)/,/\*\( *[0-9]\+,\)*\( *[0-9]\+ \/\)\{1\}/d' $fold > ! $fnew
mv -f $fnew $file_f

# rename LU_ variables, BLOCK DATA JACOBIAN_SPARSE_DATA and MONITOR_DATA
# ----------------------------------------------------------------------
set ch_VAR=("LU_ICOL" "LU_CROW" "LU_DIAG" "JACOBIAN_SPARSE_DATA" "MONITOR_DATA")
set fold=$file_f
set inn=0
foreach i ($ch_VAR)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_f



# ==============
# work on file_h
# ==============

# delete unnecessary variables
#  (these ones are only declared)
set del_VAR=("SUN" "TEMP" "RTOLS" "TSTART" "TEND" "CFACTOR" "LOOKAT" "MONITOR" "SMASS" "EQN_TAGS")
set fold=$file_h[1]
set inn=0
foreach i ($del_VAR)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed '/C '"$i"' -/,/\/ '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_h[1]

# delete unnecessary variables
#  (these one are defined as parameters, the searched pattern is different)
set del_VAR2=("NONZERO" "CNVAR" "NLOOKAT" "NMONITOR" "NMASS")
set fold=$file_h[2]
set inn=0
foreach i ($del_VAR2)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed '/C '"$i"' -/,/PARAMETER ( '"$i"' =/d' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_h[2]

# change names of common blocks *GDATA
set ch_VAR=("\/GDATA" "\/CHARGDATA" "\/INTGDATA")
set fold=$file_h[1]
set inn=0
foreach i ($ch_VAR)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_h[1]



# ==============
# work on file_s
# ==============
# delete unnecessary variables
set del_VAR3=("LU_IROW")
set fold=$file_s
set inn=0
foreach i ($del_VAR3)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed '/C '"$i"' -/,/\/ '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_s

# change names of common block SDATA
set fold=$file_s
set fnew=tmp$inn
sed 's/SDATA/&_'"$appendix"'/g' $fold > ! $fnew
set fold=$fnew
mv -f $fnew $file_s

# change names of LU_ variables
set ch_VAR2=("LU_ICOL" "LU_CROW" "LU_DIAG")
set fold=$file_s
set inn=0
foreach i ($ch_VAR2)
  @ inn=$inn + 1
  set fnew=tmp$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end
mv -f $fnew $file_s



# ========
# clean up
# ========
rm -f tmp*
unset file_f file_h file_s appendix
unset fold fnew inn
unset del_SR ch_SR ch_VAR del_VAR del_VAR2 del_VAR3 ch_VAR2

