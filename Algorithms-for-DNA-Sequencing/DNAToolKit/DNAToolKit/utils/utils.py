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
