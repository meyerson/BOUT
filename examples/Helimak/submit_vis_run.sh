qsub -l h_rt=1:00:00 -V -q development -pe 1way 16 -N CLM_32x32_dt10_vis -cwd ../run_idl
