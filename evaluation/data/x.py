seq_2_score_2_mode_dict = {}

with open('evaluation/data/scores_all.tsv') as f:
    for line in f:
        set, seq, mode, score = line.strip().split('\t')
        if mode != 'bitpal-scored':
            if set+seq not in seq_2_score_2_mode_dict:
                seq_2_score_2_mode_dict[set+seq] = {}
            if score not in seq_2_score_2_mode_dict[set+seq]:
                seq_2_score_2_mode_dict[set+seq][score] = []
            seq_2_score_2_mode_dict[set+seq][score].append(mode)

for key, score_2_mode_dict in seq_2_score_2_mode_dict.items():
    if len(score_2_mode_dict) > 0:
        print(key, score_2_mode_dict)
