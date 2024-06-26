# Usage:
#    cat *.log | python3 scripts/log2info.py

import sys

seq_name = ''
algorithm = ''
elapsed_wall_clock_time = ''
max_resident_set_size = ''
num_replicate = ''

error = False

for line in sys.stdin:
    if 'Command terminated by signal' in line:
        error = True
    elif 'Time.Alignment' in line:
        alignment_time, unit_of_measure = line.strip().split('Time.Alignment    ')[-1].strip().split(' ')[:2]
        alignment_time = float(alignment_time)
        if unit_of_measure == 'm':
            elapsed_wall_clock_time = alignment_time * 1000000000.0 * 60.0
        elif unit_of_measure == 's':
            elapsed_wall_clock_time = alignment_time * 1000000000.0
        elif unit_of_measure == 'ms':
            elapsed_wall_clock_time = alignment_time * 1000000.0
        elif unit_of_measure == 'us':
            elapsed_wall_clock_time = alignment_time * 1000.0
        else: #if unit_of_measure == 'ns':
            elapsed_wall_clock_time = alignment_time
        elapsed_wall_clock_time = str(elapsed_wall_clock_time)

        max_resident_set_size = ''
    elif 'Command being timed' in line:
        seq_name = line.split('.seq ')[0].split('/')[-1]

        if 'gap-affine-wfa' in line:
            if '--wfa-memory-mode' in line:
                algorithm = 'wfa-' + line.split('--wfa-memory-mode ')[1].split(' ')[0].strip()
            else:
                algorithm = 'wfa-ultralow'
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
            if not num_replicate:
                num_replicate = '1'
            print('\t'.join([seq_name, algorithm, num_replicate, elapsed_wall_clock_time, 'time_ns']))
            print('\t'.join([seq_name, algorithm, num_replicate, max_resident_set_size, 'memory_kb']))

        seq_name = ''
        algorithm = ''
        elapsed_wall_clock_time = ''
        max_resident_set_size = ''
        num_replicate = ''

        error = False
    elif 'Replicate ' in line:
        num_replicate = line.strip().split(' ')[-1]
