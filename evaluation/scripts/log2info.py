# Usage:
#    cat *.log | python3 scripts/log2info.py

import sys

seq_name = ''
algorithm = ''
elapsed_wall_clock_time = ''
max_resident_set_size = ''

error = False

for line in sys.stdin:
    if 'Command terminated by signal' in line:
        error = True
    elif 'Time.Alignment' in line:
        elapsed_wall_clock_time = line.strip().split('Time.Alignment    ')[-1].split(' (s)')[0].strip()

        max_resident_set_size = ''
    elif 'Command being timed' in line:
        seq_name = line.split('.seq ')[0].split('/')[-1]

        if 'gap-affine-wfa' in line:
            if 'wfa-bidirectional' in line:
                algorithm = 'biwfa'
            else:
                algorithm = 'wfa-' + line.split('--wfa-memory-mode ')[1].split(' ')[0].strip()
        else:
            algorithm = line.split(' -a ')[1].split(' ')[0].strip()

        max_resident_set_size = ''
    elif 'Maximum resident set size' in line:
        max_resident_set_size = line.strip().split('): ')[-1].strip()

        if error:
            elapsed_wall_clock_time = 'nan'
            max_resident_set_size = 'nan'

        # Check if the basic information is available
        if seq_name and algorithm:
            print('\t'.join([seq_name, algorithm, elapsed_wall_clock_time, 'time_s']))
            print('\t'.join([seq_name, algorithm, max_resident_set_size, 'memory_kb']))

        seq_name = ''
        algorithm = ''
        elapsed_wall_clock_time = ''
        max_resident_set_size = ''

        error = False
