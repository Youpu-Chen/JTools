'''
Description: The module could be applied to conduct basic Needleman-Wunsch algorithm and Smith-Waterman algorithm,
and print the best alignment to terminal as well as visualize the score/trace matrix in a new win
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


global match
global mismatch
global gap_open
global gap_extension
global MIN

match = 1
mismatch = -1
gap_open = -2
gap_extension = -1
MIN = -float("inf")



def identity_match(base1, base2, mat, msm):
    '''this function is used to compare the bases and return match point or mismatch point'''
    if base1 == base2:
        return mat
    else:
        return msm


class PairwiseSequenceAlignment():    
    def __init__(self, seq1, seq2, mat=match, mism=mismatch, g_o=gap_open, g_e=gap_extension, local=False) -> None:
        # Input sequences
        self.seq1 = seq1.upper()  # convert all the bases into Upper format
        self.seq2 = seq2.upper()

        # Set parameters, if not, it will use preset values directly
        self.mat = mat
        self.mism = mism
        self.g_o = g_o
        self.g_e = g_e
        self.local = local


    def print_sequence_len(self):
        print(f'The seq1 has length of {len(self.seq1)}')
        print(f'The seq2 has length of {len(self.seq2)}')

    
    def createscorematrix(self):
        '''this function is used to generate the original score function'''
        n = len(self.seq1)
        m = len(self.seq2)
        # Create match matrix, x matrix and y matrix
        m_mat = np.zeros((n+1, m+1), dtype=float)
        x_mat = np.zeros((n+1, m+1), dtype=float)
        y_mat = np.zeros((n+1, m+1), dtype=float)

        # Create trace matrix
        dt = np.dtype([('diagonal', np.str_, 1),
                    ('up', np.str_, 1), ('left', np.str_, 1)])
        trace_matrix = np.zeros((n+1, m+1), dtype=dt)

        return m_mat, x_mat, y_mat, trace_matrix
    
    def matrix_init(self, m_mat, x_mat, y_mat, trace_matrix):
        '''this function conduct the score matrix initialization'''
        local = self.local
        
        '''Global Alignment'''
        if local == False:
            '''trace matrix filling'''
            for i in range(0, len(trace_matrix)):
                trace_matrix[i][0]['up'] = 'U'
            for j in range(0, len(trace_matrix[0])):
                trace_matrix[0][j]['left'] = 'L'

            '''match matrix filling'''
            for i in range(0, len(m_mat)):
                for j in range(0, len(m_mat[0])):
                    if i == 0 and j == 0:
                        m_mat[i][j] = 0
                    elif i == 0 and j > 0:
                        m_mat[i][j] = MIN
                    elif i > 0 and j == 0:
                        m_mat[i][j] = MIN

            '''x_matrix filling'''
            for i in range(0, len(x_mat)):
                for j in range(0, len(x_mat[0])):
                    if i == 0 and j == 0:
                        x_mat[i][j]=0
                    if i > 0 and j == 0:
                        x_mat[i][j] = gap_open+gap_extension*(i-1)
            x_first_row = [0]
            x_first_row.extend([MIN]*(len(x_mat[0])-1))  # type: ignore
            x_mat[0] = x_first_row

            '''y_matrix filling'''
            for i in range(0, len(y_mat)):
                for j in range(0, len(y_mat[0])):
                    if i == 0 and j == 0:
                        y_mat[i][j]=0
                    elif i > 0 and j == 0: 
                        y_mat[i][j] = MIN
            y_first_row = [0]
            y_first_row.extend([gap_open+gap_extension*(i-1) for i in range(1, len(y_mat[0]))])
            y_mat[0] = y_first_row

            return m_mat, x_mat, y_mat, trace_matrix
        
        '''Local Alignment: Initialization step for Smith-Watermen is useless'''
        if local == True:
            return m_mat, x_mat, y_mat, trace_matrix
    
    
    def matrix_filling(self, m_mat, x_mat, y_mat, trace_matrix):
        '''this function is used to create the scoring matrix using three dynamic programming,
        and building a tracing matrix to restore the paths for the retrieve of aliignments'''
        '''Global Alignment Activation'''
        mat = self.mat
        mism = self.mism
        local = self.local

        if local == False:
            # Filling score matrix and record the trace
            for i in range(1, len(m_mat)):
                for j in range(1, len(m_mat[0])):
                    m_mat[i][j] = max(
                        m_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism),
                        x_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism),
                        y_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism)
                    )
                    x_mat[i][j] = max(m_mat[i-1][j] + gap_open, x_mat[i-1][j] + gap_extension)
                    y_mat[i][j] = max(m_mat[i][j-1] + gap_open, y_mat[i][j-1] + gap_extension)

            # Take the greatest values in these three matrix,
            # merge into one matrix,
            # and record the path
            new_mat = np.zeros((len(m_mat), len(m_mat[0])), dtype=int)
            for i in range(0, len(m_mat)):
                for j in range(0, len(m_mat[0])):
                    new_mat[i][j] = max(m_mat[i][j], x_mat[i][j], y_mat[i][j])
                    # Fill the trace matrix
                    # Note: from match/mismatch is 0, from x_mat (open a gap in seq2) is 1, from y_mat (open a gap in seq1)
                    if m_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['diagonal'] = 'D'
                    elif x_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['up'] = 'U'
                    elif y_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['left'] = 'L'

            return new_mat, trace_matrix
        
        '''Local Alignment Activation'''
        if local == True:
            # Filling score matrix and record the trace
            for i in range(1, len(m_mat)):
                for j in range(1, len(m_mat[0])):
                    m_mat[i][j] = max(
                        m_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism),
                        x_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism),
                        y_mat[i-1][j-1] + identity_match(self.seq1[i-1], self.seq2[j-1], mat, mism),
                        0
                    )
                    x_mat[i][j] = max(m_mat[i-1][j] + gap_open, x_mat[i-1][j] + gap_extension)
                    y_mat[i][j] = max(m_mat[i][j-1] + gap_open, y_mat[i][j-1] + gap_extension)

            # Take the Greatest values in these three matrix
            new_mat = np.zeros((len(m_mat), len(m_mat[0])), dtype=int)
            for i in range(0, len(m_mat)):
                for j in range(0, len(m_mat[0])):
                    new_mat[i][j] = max(m_mat[i][j], x_mat[i][j], y_mat[i][j])
                    if new_mat[i][j] == 0:   # only non-zero to zero direction will be considered
                        continue
                    if m_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['diagonal'] = 'D'
                    elif x_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['up'] = 'U'
                    elif y_mat[i][j] == max(m_mat[i][j], x_mat[i][j], y_mat[i][j]):
                        trace_matrix[i][j]['left'] = 'L'

            return new_mat, trace_matrix

    def global_backtracking(self, trace_matrix, score_matrix):
        '''this function is used to trace back the input matrix and output the final alignment
        Note: 
        - the input matrix is trace matrix
        - Only print out the best alignment'''
        # Create path list
        path_list = []

        # Retrieving
        ti = len(self.seq1)
        tj = len(self.seq2)
        alignment1 = ''
        alignment2 = ''
        while (ti > 0 or tj > 0):
            # Choose to go left, up or diagonal
            cell = trace_matrix[ti][tj]
            if 'D' in cell:
                alignment1 = self.seq1[ti-1] + alignment1
                alignment2 = self.seq2[tj-1] + alignment2
                path_list.append([tj, ti, tj-1, ti-1])
                ti -= 1
                tj -= 1
                
            elif 'U' in cell:
                alignment1 = self.seq1[ti-1] + alignment1
                alignment2 = '-' + alignment2
                path_list.append([tj, ti, tj, ti-1])
                ti -= 1
                
            elif 'L' in cell:
                alignment1 = '-' + alignment1
                alignment2 = self.seq2[tj-1] + alignment2
                path_list.append([tj, ti, tj-1, ti])
                tj -= 1
                
        info = f"======The Global======\n     {alignment1}\n     {alignment2}\nSCORE: {score_matrix[len(score_matrix)-1][len(score_matrix[0])-1]}"
        print(info)

        return path_list
    
    def local_backtracking(self, trace_matrix, score_matrix):
        '''this function does backtracking like FUNCTION global_backtracking, but in the way of local aligment
        Note: Only print out the best alignment (start from the best score)
        '''
        # Create path list
        path_list = []

        # Convert score matrix into Numpy array to find maximum value
        # new_score_matrix = np.array(score_matrix)
        pos = np.unravel_index(np.argmax(score_matrix, axis=None), score_matrix.shape)  # retrieve the maximum value
        
        # Retrieving
        ti = pos[0]
        tj = pos[1]
        print(f'{ti}\t{tj}')
        alignment1 = ''
        alignment2 = ''

        while (ti > 0 or tj > 0):
            if score_matrix[ti][tj] == 0: # stop local alignment back tracking when 0 values met
                break
            cell = trace_matrix[ti][tj]
            if 'D' in cell:
                alignment1 = self.seq1[ti-1] + alignment1
                alignment2 = self.seq2[tj-1] + alignment2
                path_list.append([tj, ti, tj-1, ti-1])
                ti -= 1
                tj -= 1
                
            elif 'U' in cell:
                alignment1 = self.seq1[ti-1] + alignment1
                alignment2 = '-' + alignment2
                path_list.append([tj, ti, tj, ti-1])
                ti -= 1
                
            elif 'L' in cell:
                alignment1 = '-' + alignment1
                alignment2 = self.seq2[tj-1] + alignment2
                path_list.append([tj, ti, tj-1, ti])
                tj -= 1
                
        info = f"======The Local======\n     {alignment1}\n     {alignment2}\nSCORE: {np.ndarray.max(score_matrix)}"
        print(info)

        return path_list

    def path_plot(self, score_matrix, trace_matrix, path_list):
        '''this function is use to visualize the score_matrix and trace_matrix in combination
        Note: 
        - adapted from kevinah95/plot_needleman_wunsch.py
        - links: https://gist.github.com/kevinah95/f6f9d16ebd17a3a28208471f4e4bb878'''
        # 
        # print(path_list)
        path_list = np.array(path_list)
        
        # Set plot parameters
        plt.rcParams["figure.figsize"] = 4, 5
        param = {"grid.linewidth": 1.6,
                "grid.color": "lightgray",
                "axes.linewidth": 1.6,
                "axes.edgecolor": "lightgray"}
        plt.rcParams.update(param)
        
        # Set matrix plot labels
        headv = self.seq1
        headh = self.seq2

        # Plot
        fig, ax = plt.subplots()
        ax.set_xlim(-1.5, len(score_matrix[0]) - .5)
        ax.set_ylim(-1.5, len(score_matrix) - .5)
        ax.invert_yaxis()
        # for i in range(score_matrix.shape[0]):
        for i in range(len(score_matrix)):
            # for j in range(score_matrix.shape[1]):
            for j in range(len(score_matrix[0])):
                ax.text(j, i, score_matrix[i, j], ha="center", va="center")
        for i, l in enumerate(headh):
            ax.text(i + 1, -1, l, ha="center", va="center", fontweight="semibold")
        for i, l in enumerate(headv):
            ax.text(-1, i + 1, l, ha="center", va="center", fontweight="semibold")

        ax.xaxis.set_minor_locator(ticker.FixedLocator(
            np.arange(-1.5, len(score_matrix[0]) - .5, 1)))
        ax.yaxis.set_minor_locator(ticker.FixedLocator(
            np.arange(-1.5, len(score_matrix[0]) - .5, 1)))
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        left="off", right="off", labelbottom='off', labelleft='off')
        ax.grid(True, which='minor')

        # Set arrorw tags property
        arrowprops = dict(facecolor='blue', alpha=0.5, lw=0,
                        shrink=0.2, width=2, headwidth=7, headlength=7)

        # Plot all paths
        for i in range(1, trace_matrix.shape[0]):
            for j in range(1, trace_matrix.shape[1]):
                # if i == 0 and j == 0:
                #     continue
                if(trace_matrix[i][j]['left'] != ''):
                    ax.annotate("", xy=(j - 1, i),
                                xytext=(j, i), arrowprops=arrowprops)
                if(trace_matrix[i][j]['diagonal'] != ''):
                    ax.annotate("", xy=(j - 1, i - 1),
                                xytext=(j, i), arrowprops=arrowprops)
                if(trace_matrix[i][j]['up'] != ''):
                    ax.annotate("", xy=(j, i - 1),
                                xytext=(j, i), arrowprops=arrowprops)
        # optimal path
        arrowprops.update(facecolor='crimson')
        for i in range(path_list.shape[0]):
            ax.annotate("", xy=path_list[i, 2:],
                        xytext=path_list[i, :2], arrowprops=arrowprops)
        plt.show()
        
        return fig

    @staticmethod
    def Runalignment(seq1, seq2, m, mism, g_o, g_e, local):
        '''main function to run sequence alignment module'''
        if local == False:
            alignment = PairwiseSequenceAlignment(seq1, seq2, m, mism, g_o, g_e, local=False)
            m_mat, x_mat, y_mat, trace_matrix = alignment.createscorematrix()
            m_mat, x_mat, y_mat, trace_matrix = alignment.matrix_init(m_mat, x_mat, y_mat, trace_matrix)
            score_matrix, trace_matrix = alignment.matrix_filling(m_mat, x_mat, y_mat, trace_matrix)
            print(score_matrix)
            path_list = alignment.global_backtracking(trace_matrix, score_matrix)
            alignment.path_plot(score_matrix, trace_matrix, path_list)
        
        if local == True:
            alignment = PairwiseSequenceAlignment(seq1, seq2, m, mism, g_o, g_e, local=True)
            m_mat, x_mat, y_mat, trace_matrix = alignment.createscorematrix()
            m_mat, x_mat, y_mat, trace_matrix = alignment.matrix_init(m_mat, x_mat, y_mat, trace_matrix)
            score_matrix, trace_matrix = alignment.matrix_filling(m_mat, x_mat, y_mat, trace_matrix)
            # print(score_matrix)
            path_list = alignment.local_backtracking(trace_matrix, score_matrix)
            alignment.path_plot(score_matrix, trace_matrix, path_list)


if __name__ == "__main__":
    # ---------------------------------------------------------
    # Pairwise Sequence Alignment
    # ---------------------------------------------------------
    print('------------------------------- Pairwise Sequence Alignment -------------------------------')
    seq1 = "TCGTAGACGA"
    seq2 = "ATAGAATGCGG"
    # Global Alignment
    PairwiseSequenceAlignment.Runalignment(seq1, seq2, 1, -1, -2, -1, local=False)

    # Local Alignment
    PairwiseSequenceAlignment.Runalignment(seq1, seq2, 1, -1, -2, -1, local=True)