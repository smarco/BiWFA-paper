# Usage:
#    cat *.log | python3 scripts/log2info.py

import sys

seq_name = ''
algorithm = ''
elapsed_wall_clock_time = ''
max_resident_set_size = ''
num_replicate = ''
error_rate = 'na'
time_unit_of_measure = ''

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

        time_unit_of_measure = 'time_s'
        max_resident_set_size = ''
    elif 'InfoAlignment' in line:
        # InfoAlignment: 162789	1225.355969	-43276	5.21
        # QUERY_LENTH, AVERAGE_TIME_MS[, SCORE, ERROR_RATE]
        info = line.strip().split('InfoAlignment: ')[1].split('\t')
        elapsed_wall_clock_time = info[1]
        if len(info) > 2:
            # WFA-based alignments report also the score and the error rate
            error_rate = info[3]

        time_unit_of_measure = 'time_ms'
        max_resident_set_size = ''
    elif 'Command being timed' in line:
        seq_name = line.split('.seq ')[0].split('/')[-1]

        if 'gap-affine-wfa' in line:
            # The old align_benchmark had the --wfa-bidirectional flag
            if 'wfa-bidirectional' in line or 'biwfa' in line:
                algorithm = 'biwfa'
            else:
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
            print('\t'.join([seq_name, algorithm, num_replicate, elapsed_wall_clock_time, time_unit_of_measure, error_rate]))
            print('\t'.join([seq_name, algorithm, num_replicate, max_resident_set_size, 'memory_kb', error_rate]))

        seq_name = ''
        algorithm = ''
        elapsed_wall_clock_time = ''
        max_resident_set_size = ''
        num_replicate = ''
        error_rate = 'na'
        time_unit_of_measure = ''

        error = False
    elif 'Replicate ' in line:
        num_replicate = line.strip().split(' ')[-1]

