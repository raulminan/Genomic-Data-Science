"""functions for global aligment"""

def global_alignment(x: str, y: str) -> int:
    """Performs global aligment on two DNA strings, returning the edit distance 
    between them. To caculate the edit distance, a penalty matrix is used, giving
    transitions, transversions and skips a penalty of 2, 4 and 8, respectively

    Parameters
    ----------
    x : str
        first dna string
    y : str
        second dna string

    Returns
    -------
    int
        edit distance between x and y 
    """
    alphabet = ["A", "C", "G", "T"]
    score = [
        [0, 4, 2, 4, 8], # penalties for A
        [4, 0, 4, 2, 8], # penalties for C
        [2, 4, 0, 4, 8], # penalties for G
        [4, 2, 4, 0, 8], # penalties for T
        [8, 8, 8, 8, 8]  # penalties for -
        ] # TODO create build_penalty_matrix() function to assign custom penalties

    D = []

    # initialize empty matrix
    for i in range(len(x) + 1):
        D.append([0] * (len(y)+1))
    
    # fill first row and first column
    for i in range(1, len(x) + 1):
        # alphabet.index gives the index of the current char in x (what row to look at),
        # then return the last one (skip, for a penalty of 8), since the 1st row
        # corresponds to skipping the first i chars in y
        D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]


    for i in range(1, len(y) + 1):
        # same thing as above applies here
        D[0][i] = D[0][i-1] + score[-1][alphabet.index(y[i-1])]
        

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            dist_hor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            dist_ver = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                dist_diag = D[i-1][j-1]
            else:
                dist_diag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            D[i][j] = min(dist_hor, dist_ver, dist_diag)

    return D[-1][-1]