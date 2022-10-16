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
        tmp = line.strip().split('Time.Alignment    ')[-1]
        if '(s)' in tmp:
            # The old align_benchmark emits all runtimes in seconds
            elapsed_wall_clock_time = tmp.split(' (s)')[0].strip()
        else:
            alignment_time, unit_of_measure = tmp.strip().split(' ')[:2]
            alignment_time = float(alignment_time)
            if unit_of_measure == 'm':
                elapsed_wall_clock_time = alignment_time * 60.0
            elif unit_of_measure == 's':
                elapsed_wall_clock_time = alignment_time
            elif unit_of_measure == 'ms':
                elapsed_wall_clock_time = alignment_time / 1000.0
            elif unit_of_measure == 'us':
                elapsed_wall_clock_time = 0.0
            else: #if unit_of_measure == 'ns':
                elapsed_wall_clock_time = 0.0
            elapsed_wall_clock_time = str(elapsed_wall_clock_time)

        max_resident_set_size = ''
    elif 'Command being timed' in line:
        seq_name = line.split('.seq ')[0].split('/')[-1]

        if 'gap-affine-wfa' in line:
            # The old align_benchmark had the --wfa-bidirectional flag
            #if 'wfa-bidirectional' in line:
            #    algorithm = 'biwfa'
            #else:
            algorithm = 'wfa-' + line.split('--wfa-memory-mode ')[1].split(' ')[0].strip()
        else:
            algorithm = line.split(' -a ')[1].split(' ')[0].strip()

        max_resident_set_size = ''
    elif 'Maximum resident set size' in line:
        max_resident_set_size = line.strip().split('): ')[-1].strip()
    elif 'Exit status' in line:
        error = error or line.strip().split(': ')[-1] != '0'

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
