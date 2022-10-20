"""utility functions for DNAToolKit"""

def edit_distance(x: str, y: str) -> int:
    """Computes the edit distance between two strings using dynamic programming.

    The function fills up a matrix D of the edit distances of every prefix in x and y:
    In this matrix, every element D[i, j] represents the edit distance between x[:i] and
    y[:j].

    Matrix D:

     e G C T A T A C
   e 0 1 2 3 4 5 6 7
   G 1 
   C 2
   G 3
   T 4
   A 5
   T 6
   G 7

   The final edit distance is given by the element in the bottom right (D[-1][-1])

    Parameters
    ----------
    x : str
        first string
    y : str
        second string

    Returns
    -------
    int
        edit distance between x and y
    """
    D = []

    # initialize empty matrix
    for i in range(len(x) + 1):
        D.append([0] * (len(y)+1))
    
    # fill first row and first columns
    for i in range(len(x) + 1):
        D[i][0] = i

    for i in range(len(y) + 1):
        D[0][i] = i

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            dist_hor = D[i][j-1] + 1
            dist_ver = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                dist_diag = D[i-1][j-1]
            else:
                dist_diag = D[i-1][j-1] + 1  

            D[i][j] = min(dist_hor, dist_ver, dist_diag)

    return D[-1][-1]

def edit_distance_fewest_edits(pattern: str, text: str) -> int:
    """Computes the edit distance of the match between P and T with 
    the fewest edits

    Parameters
    ----------
    pattern : str
        pattern to match against text
    text : str
        text to find pattern in

    Returns
    -------
    int
        edit distance
    """
    D = []
    # init 0-filled matrix
    for i in range(len(pattern) + 1):
        D.append([0] * (len(text) + 1))
    
    # fill first column
    for i in range(len(pattern)+1):
        D[i][0] = i 
    
    # note: 
    # There's no need to init the first row with ascending values
    # We'll init the first row with all 0 because we don't know ahead
    # of time where pattern will occur within text, so every offset
    # is equally likely.

    for i in range(1, len(pattern) + 1):
        for j in range(1, len(text) + 1):
            dist_hor = D[i][j-1] + 1
            dist_ver = D[i-1][j] + 1
            if pattern[i-1] == text[j-1]:
                dist_diag = D[i-1][j-1]
            else:
                dist_diag = D[i-1][j-1] + 1  

            D[i][j] = min(dist_hor, dist_ver, dist_diag)

    # the lowest edit distance is given my the min in the last row
    min_distance = min(D[-1])
    return min_distance