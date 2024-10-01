def find_non_conflicting_brackets(pair_list, bracket_types):
    """
    对于每一对配对，找到一个不产生矛盾的括号表示。
    """
    # 初始化配对的括号类型
    pair_bracket_map = {}
    for i, j in pair_list:
        for left_bracket, right_bracket in bracket_types:
            # 检查当前括号是否产生矛盾
            conflict = False
            for k, l in pair_list:
                if i < k < j < l or k < i < l < j:
                    # 存在交叉配对，检查是否已经使用了当前的括号类型
                    if pair_bracket_map.get((k, l), (None, None)) == (left_bracket, right_bracket):
                        conflict = True
                        break
            if not conflict:
                pair_bracket_map[(i, j)] = (left_bracket, right_bracket)
                break
    return pair_bracket_map

def matrix_to_dbn(pair_matrix):
    """
    将二维配对矩阵转换为RNA二级结构的DBN格式，处理交叉配对情况。
    """
    n = len(pair_matrix)
    dbn_list = ['.'] * n
    pair_list = [(i, j) for i in range(n) for j in range(i+1, n) if pair_matrix[i][j] == 1]
    bracket_types = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>'), 
                     ('A', 'a'), ('B', 'b'), ('C', 'c'), ('D', 'd')]  # 可以根据需要添加更多类型
    
    pair_bracket_map = find_non_conflicting_brackets(pair_list, bracket_types)
    
    for (i, j), (left_bracket, right_bracket) in pair_bracket_map.items():
        dbn_list[i] = left_bracket
        dbn_list[j] = right_bracket
    
    return ''.join(dbn_list)

# 示例使用
pair_matrix = [
    [0, 0, 0, 0, 0, 1, 0, 0],  # 5 with 7
    [0, 0, 0, 0, 0, 0, 1, 0],  # 6 with 1
    [0, 0, 0, 0, 0, 0, 0, 1],  # 7 with 2
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],  # 5 with 7
    [0, 1, 0, 0, 0, 0, 0, 0],  # 6 with 1
    [0, 0, 1, 0, 0, 0, 0, 0]   # 7 with 2
]

print(matrix_to_dbn(pair_matrix))

