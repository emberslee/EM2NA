import os
import numpy as np
try:
    from numba import jit
    print("# numba.jit successfully imported")
    use_jit = True
    use_jit = False
    if use_jit:
        print("# Use numba.jit to accelerate")
    else:
        print("# Do not used numba.jit to acclerate")
except:
    use_jit = False
    print("# numba.jit not available")

def jit_if_available(func):
    if use_jit:
        return jit(func, nogil=True, nopython=True)
    else:
        return func

# Alignments considering geometry of points
@jit_if_available
def global_align_score(
        ldp_gap, seq_gap, matrix, ldp_list, seq_list, pointer, scratch,
        seq_gap_count, ldp_gap_count, ldp_distance_array, ldp_prev_connect,
        seq_gap_total_limit=2,
    ):
    """
        ldp_gap : gap penalty between ldps default=10
        seq_gap : gap penalty between seqs default=25
        matrix : match probability for [i, j] 
        ldp_list : list for ldps
        seq_list : list for seqs
        pointer : directions to recover path default=np.zeros((M, N))
        scratch : score matrix default=np.zeros((M, N)) and scratch[..., 0]=-1e6

        seq_gap_count : counts for seq gaps
        ldp_gap_count : counts for ldp gaps
        ldp_distance_array : ldp distances
        ldp_prev_connect : previous connected ldp
        
        seq_gap_total_limit : gap limits allowed in seq
    """

    # Limits
    low_cut_off_distance = 4
    high_cut_off_distance = 8
    mean_distance = 6

    # For each point
    for i in range(1, len(ldp_list) + 1):
        for j in range(1, len(seq_list) + 1):
            reasonable_skip_ldp = 0  # Init whether we can skip points

            # From diagonal
            if i == 1:
                current_best = scratch[i - 1, j - 1] + matrix[i, j]
            else:
                # Get distance from previous point
                prev_neighbor = int(ldp_prev_connect[i - 1, j - 1]) if ldp_prev_connect[i - 1, j - 1] != -1 else i - 1
                cur_distance = ldp_distance_array[prev_neighbor, i]

                # Penalty by distance, matrix and ss_matrix
                if cur_distance <= low_cut_off_distance or cur_distance >= high_cut_off_distance:
                    current_best = scratch[i - 1, j - 1] + matrix[i, j] - ldp_gap * abs(cur_distance - mean_distance) ** 2
                else:
                    current_best = scratch[i - 1, j - 1] + matrix[i, j]

            # From left
            if i == 1:
                left_cur = scratch[i - 1, j]
            else:
                prev_neighbor = int(ldp_prev_connect[i - 1, j]) if ldp_prev_connect[i - 1, j] != -1 else i - 1
                cur_distance = ldp_distance_array[prev_neighbor, i]

                # Some controls whther we can skip ldp
                if i < len(ldp_list) and cur_distance + ldp_distance_array[i, i + 1] < high_cut_off_distance:
                    reasonable_skip_ldp = 1
                if cur_distance <= low_cut_off_distance:
                    reasonable_skip_ldp = 1

                if reasonable_skip_ldp:
                    left_cur = scratch[i - 1, j]
                else:
                    left_cur = scratch[i - 1, j] - ldp_gap * abs(cur_distance - mean_distance) ** 2
            left_best_step = 1

            # From up
            if seq_gap_count[i, j - 1] <= seq_gap_total_limit:
                up_cur = scratch[i, j - 1] - seq_gap
            else:
                up_cur = scratch[i, j - 1] - seq_gap * (seq_gap_count[i, j - 1] + 1)
            up_best_step = 1

            # Get the best score
            final_best = max(current_best, left_cur, up_cur)
            scratch[i, j] = final_best

            # Record the direction
            if final_best == up_cur:
                pointer[i, j] = 1
                seq_gap_count[i, j] = 1 + seq_gap_count[i, j - 1]
                ldp_prev_connect[i, j] = ldp_prev_connect[i, j - 1]
                ldp_gap_count[i, j] = ldp_gap_count[i, j - 1]

            elif final_best == left_cur:
                pointer[i, j] = -1
                seq_gap_count[i, j] = seq_gap_count[i - 1, j]
                if ldp_prev_connect[i - 1, j] != -1:
                    ldp_prev_connect[i, j] = ldp_prev_connect[i - 1, j]
                else:
                    ldp_prev_connect[i, j] = i - 1
                if reasonable_skip_ldp:
                    ldp_gap_count[i, j] = ldp_gap_count[i - 1, j]
                else:
                    ldp_gap_count[i, j] = ldp_gap_count[i - 1, j] + 1

            else:
                pointer[i, j] = 999999  # Move from diagonal
                seq_gap_count[i, j] = seq_gap_count[i - 1, j - 1]
                ldp_gap_count[i, j] = ldp_gap_count[i - 1, j - 1]
    return scratch, pointer, seq_gap_count, ldp_prev_connect, ldp_gap_count


def dynamic_assign_multi(updated_base_list, fragment_distance_array, current_chain,
                         path_ss=None, seq_ss=None,
                         ldp_gap_penalty=10, seq_gap_penalty=25, save_path=None, lower_bound=-1000, sort=True, top=20):
    """
    Sequence alignment

    Parameters:
    - updated_base_list (list): List for base probability (A,G,C,U) for each point
    - current_chain (list): List of current chain
    - fragment_distance_array (numpy.ndarray): Distance matrix of points
    - ldp_gap_penalty (float): Gap penalty default=10
    - seq_gap_penalty (float): Gap penalty default=25
    - save_path (str, optional): Save path

    Returns:
    - max_score_list (list): List of max scores
    - match_seq_list (list): List of matched seq
    - match_interval_list (list): List of matched interval

    Raises:
    Possible exception caused by global_align_score

    Example:
    >>> dynamic_assign_multi(updated_base_list, current_chain, 0.1, 0.2, fragment_distance_array, save_path='results')

    """
    # Mapping from digital to str
    map_dict = {0: "A", 1: "G", 2: "C", 3: "U"}
    
    # Create (M, N) match matrix
    match_matrix = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
   
    score_by_prob = True

    # 2024-01-24
    aa_match_matrix = np.array(
        [[0.20466549, 0.13293310, 0.07639524, 0.08726946],
         [0.12993139, 0.23787625, 0.07880677, 0.06438896],
         [0.07319305, 0.07168299, 0.22714286, 0.11792232],
         [0.08726946, 0.05300044, 0.17080019, 0.18721251]]
    )
    aa_match_matrix = (aa_match_matrix - 0.05) * 10

    # Fill the matrix
    if score_by_prob:
        for i in range(len(updated_base_list)):
            current_prob = np.asarray(updated_base_list[i], dtype=np.float32)
            max_p = np.argmax(np.asarray(current_prob), axis=0)
            current_prob = np.asarray(current_prob) * aa_match_matrix[max_p]

            for j in range(len(current_chain)):
                cur_label = int(current_chain[j])
                match_matrix[i + 1, j + 1] = current_prob[cur_label] * 100 # [0.0, 1.0] * 100.0
    else:
        pred_base_types = np.argmax(np.asarray(updated_base_list, dtype=np.float32), axis=-1) # (L, )
        seq_base_types = np.asarray(current_chain, dtype=np.int32) # (L, )
        for i in range(len(pred_base_types)):
            pred_base_type = pred_base_types[i]
            for j in range(len(seq_base_types)):
                seq_base_type = seq_base_types[j]
                match_matrix[i + 1, j + 1] = aa_match_matrix[seq_base_type, pred_base_type] * 100

    # Fill the ss_matrix
    ss_match_matrix = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
    if  path_ss is not None and \
        seq_ss is not None:
        ss_score = np.asarray(
            [
                [0.0, 0.0],
                [0.0, 1.0],
            ]
        )
        for i in range(len(updated_base_list)):
            cur_path_ss_type = path_ss[i]
            for j in range(len(current_chain)):
                cur_seq_ss_type = seq_ss[j]
                ss_match_matrix[i + 1, j + 1] = ss_score[cur_path_ss_type][cur_seq_ss_type] * 200
        match_matrix = match_matrix + ss_match_matrix

    # Save match matrix
    if save_path is not None:
        score_path = os.path.join(save_path, "match_score.txt")
        np.savetxt(score_path, match_matrix)
    
    match_score = match_matrix
    pointer = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
    scratch = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
    
    # Init scratch matrix
    for k in range(1, len(updated_base_list) + 1):
        scratch[k, 0] = -999999  # Not allow gaps
    
    fragment_final_distance = np.zeros([len(fragment_distance_array) + 1, len(fragment_distance_array) + 1])
    fragment_final_distance[1:, 1:] = fragment_distance_array
    fragment_final_distance[0, :] = 20
    fragment_final_distance[:, 0] = 20
    
    seq_gap_accumulation = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
    ldp_gap_accumulation = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1])
   
    for k in range(1, len(updated_base_list) + 1):
        ldp_gap_accumulation[k, 0] = k
 
    seq_gap_total_limit = int(len(updated_base_list) * 0.05)
    ldp_prev_connect = np.zeros([len(updated_base_list) + 1, len(current_chain) + 1]) - 1
    ldp_prev_connect[0] = 0
    
    scratch, pointer, seq_gap_count, ldp_prev_connect, ldp_gap_count = global_align_score(
        ldp_gap_penalty, seq_gap_penalty, 
        match_matrix, 
        np.arange(len(updated_base_list)), np.arange(len(current_chain)), 
        pointer, scratch, 
        seq_gap_accumulation, ldp_gap_accumulation, fragment_final_distance, ldp_prev_connect, seq_gap_total_limit,
    )

    # Save infos
    if save_path is not None:
        score_path = os.path.join(save_path, "optimal_score.txt")
        np.savetxt(score_path, scratch)
        score_path = os.path.join(save_path, "optimal_direction.txt")
        np.savetxt(score_path, pointer)
        score_path = os.path.join(save_path, "sequence_miss_count.txt")
        np.savetxt(score_path, seq_gap_count)
        score_path = os.path.join(save_path, "ldp_prev_connect.txt")
        np.savetxt(score_path, ldp_prev_connect)
        score_path = os.path.join(save_path, "ldp_distance_array.txt")
        np.savetxt(score_path, fragment_final_distance)
        score_path = os.path.join(save_path, "ldp_gap_count.txt")
        np.savetxt(score_path, ldp_gap_count)
    
    input_seq_line = ""
    
    for i in range(len(updated_base_list)):
        choice = int(np.argmax(updated_base_list[i]))
        input_seq_line += map_dict[choice]
    
    # Backtrace all choices that have >0 scores
    max_score_candidate_index = np.argwhere(scratch[len(updated_base_list)] > lower_bound)
    max_score_list = []
    match_seq_list = []
    match_interval_list = []
    count_match = 0
    
    # Clip too many candidates
    if len(max_score_candidate_index) > 50:
        max_score_candidate_index = max_score_candidate_index[:50]
 
    for candidate_index in max_score_candidate_index:
        match_matrix = np.zeros(len(updated_base_list)) - 1
        candidate_index = int(candidate_index)
        max_index = candidate_index
        max_index = [len(updated_base_list), max_index]
        check_pointer = pointer[max_index[0], max_index[1]]
        cur_x = int(max_index[0])
        cur_y = int(max_index[1])
        
        while check_pointer != 0:
            if check_pointer == 999999:
                match_matrix[cur_x - 1] = cur_y - 1
                current_begin_index = cur_y - 1
                cur_x -= 1
                cur_y -= 1
            elif check_pointer > 0:
                cur_y -= int(check_pointer)
            elif check_pointer < 0:
                cur_x -= int(abs(check_pointer))
            
            check_pointer = pointer[cur_x, cur_y]
        
        if save_path is not None:
            score_path = os.path.join(save_path, "match_result%d.txt" % count_match)
            np.savetxt(score_path, match_matrix)
        
        match_seq_line = ""
        count_replace = 0
        
        for i in range(len(match_matrix)):
            if match_matrix[i] == -1:
                match_seq_line += "-"
            else:
                point_index = int(match_matrix[i])
                match_label = int(current_chain[point_index])
                current_ldp_score = updated_base_list[i]
                current_ldp_prob = current_ldp_score[match_label]
                match_seq_line += map_dict[match_label]
                
                if current_ldp_prob <= 0.25:
                    count_replace += 1
        
        current_score = float(scratch[len(updated_base_list), candidate_index])
        max_score_list.append(current_score)
        match_seq_list.append(match_seq_line)
        match_interval_list.append([current_begin_index, candidate_index])
        count_match += 1
    
    all_seq_line = ""
    
    for i in range(len(current_chain)):
        match_label = int(current_chain[i])
        all_seq_line += map_dict[match_label]

    # Sort the result
    if sort:
        indices = np.argsort(max_score_list)
        indices = indices[::-1]
        max_score_list = [max_score_list[x] for x in indices]
        match_seq_list = [match_seq_list[x] for x in indices]
        match_interval_list = [match_interval_list[x] for x in indices]
    
    # Save the match result
    if save_path is not None:
        score_path = os.path.join(save_path, "match_seq.txt")
        with open(score_path, 'w') as file:
            file.write(input_seq_line + "\n")
            
            for j in range(len(match_seq_list)):
                current_interval = match_interval_list[j]
                file.write(match_seq_list[j] + "\t" + "%.2f" % max_score_list[j] +
                           "\t%d,%d\n" % (current_interval[0], current_interval[1]))
            
            file.write(all_seq_line + "\n")

    # Return top
    if len(max_score_list) > top:
        max_score_list = max_score_list[:top]
        match_seq_list = match_seq_list[:top]
        match_interval_list = match_interval_list[:top]
    
    return max_score_list, match_seq_list, match_interval_list

