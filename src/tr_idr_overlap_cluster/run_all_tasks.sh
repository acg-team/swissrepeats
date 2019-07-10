#!/bin/bash
#
# Modules laden – kann sein dass das nicht nötig ist, bin mir aber nicht sicher
module load slurm/17.11.8
module load gcc/7.3.0
module load bzip2/1.0.6
module load freetype/2.9.1
module load curl/7.60.0
module load icu4c/60.1
module load libtiff/4.0.9
module load tcl/8.6.8
module load pcre/8.42
module load readline/7.0
module load glib/2.56.3-py3.6-pe5.26
module load cairo/1.16.0-py3.6-pe5.26
module load tk/8.6.8
module load jdk/11.0.2_9
module load ncurses/6.1
module load zlib/1.2.11
module load pango/1.41.0-py3.6-pe5.26
module load libx11/1.6.7
module load bison/3.0.5
#
Rscript one-hot-encode_for_upsetr_1-7000.R &
Rscript one-hot-encode_for_upsetr_7001-14001.R &
Rscript one-hot-encode_for_upsetr_14002-21002.R &
Rscript one-hot-encode_for_upsetr_21003-28003.R &
Rscript one-hot-encode_for_upsetr_28004-35004.R &
Rscript one-hot-encode_for_upsetr_35005-42005.R &
Rscript one-hot-encode_for_upsetr_42006-49006.R &
Rscript one-hot-encode_for_upsetr_49007-56007.R &
Rscript one-hot-encode_for_upsetr_56008-63008.R &
Rscript one-hot-encode_for_upsetr_63009-70009.R &
Rscript one-hot-encode_for_upsetr_70010-77010.R &
Rscript one-hot-encode_for_upsetr_77011-84011.R &
Rscript one-hot-encode_for_upsetr_84012-91012.R &
Rscript one-hot-encode_for_upsetr_91013-98013.R &
Rscript one-hot-encode_for_upsetr_98014-105014.R &
Rscript one-hot-encode_for_upsetr_105015-112015.R &
Rscript one-hot-encode_for_upsetr_112016-119016.R &
Rscript one-hot-encode_for_upsetr_119017-126017.R &
Rscript one-hot-encode_for_upsetr_126018-133018.R &
Rscript one-hot-encode_for_upsetr_133019-140019.R &
Rscript one-hot-encode_for_upsetr_140020-147020.R &
Rscript one-hot-encode_for_upsetr_147021-154021.R &
Rscript one-hot-encode_for_upsetr_154022-161022.R &
Rscript one-hot-encode_for_upsetr_161023-168023.R &
Rscript one-hot-encode_for_upsetr_168024-175024.R &
Rscript one-hot-encode_for_upsetr_175025-182025.R &
Rscript one-hot-encode_for_upsetr_182026-189026.R &
Rscript one-hot-encode_for_upsetr_189027-196027.R &
Rscript one-hot-encode_for_upsetr_196028-203028.R &
Rscript one-hot-encode_for_upsetr_203029-end.R &

wait
