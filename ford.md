project: nlesolver-fortran
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/nlesolver-fortran
summary: Nonlinear Equation Solver with Modern Fortran
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         protected
         private
source: true
graph: true
externalize: true
external: fmin = https://jacobwilliams.github.io/fmin
extra_mods: iso_fortran_env: https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            fmin_module: https://jacobwilliams.github.io/fmin
            lusol_ez_module: https://jacobwilliams.github.io/lusol
            lsqr_module: https://jacobwilliams.github.io/LSQR
            lsmrModule: https://jacobwilliams.github.io/LSMR

{!README.md!}
