# AUTHORS : YOUSSEF EMAD SAAD AND RAMEZ NABIL
# DATE: 13/5/2021
import operator
import pandas as pd
import matplotlib.pyplot as plt

matrix = []
# parameters
match = 1
mismatch = -1
gap = -1
# inputs
sequence_a = "GCATGCU"
sequence_b = "GATTACA"
N = len(sequence_a) + 1
M = len(sequence_b) + 1

# STEP 1: CONSTRUCTING AND FILLING THE GRID
# initialization with zeros
matrix = [[0] * N for _ in range(M)]

# initialize first row and coloumn
for i in range(M):
    matrix[i][0] = i * mismatch
for i in range(N):
    matrix[0][i] = i * mismatch

# initialize the tree dictionary, which will keep the paths from nodes
tree = {}
for i in range(1, M):
    tree.setdefault((i, 0), [(i - 1, 0)])
for i in range(1, N):
    tree.setdefault((0, i), [(0, i - 1)])

# fill the numbers of each element
for i in range(1, M):
    for j in range(1, N):

        diagonal = 0
        if sequence_a[j - 1] == sequence_b[i - 1]:
            diagonal += match
        else:
            diagonal += mismatch

        diagonal += matrix[i - 1][j - 1]

        # Getting values of neighbours
        left = matrix[i][j - 1] + gap
        up = matrix[i - 1][j] + gap
        matrix[i][j] = max(diagonal, up, left)

        # Getting indices of neighbours
        left_id = (i, j - 1)
        up_id = (i - 1, j)
        diagonal_id = (i - 1, j - 1)
        indices = [diagonal_id, up_id, left_id]

        # Getting maximum indices
        max_indices = []
        values = [diagonal, up, left]
        for k in range(len(values)):
            if matrix[i][j] == values[k]:
                max_indices.append(indices[k])

        tree.setdefault((i, j), max_indices)


# Step 2: GET THE POSSIBLE ALIGNMENTS, AND TRACE BACK TO THE FIRST ELEMENT.

# queue variable represent a path that can be taken from bottom right to top left cell
# queues variable represent all optimal paths
# tree is traversed in depth-first recursive manner
def traverse(tree, node, queue, queues):
    if node == (0, 0):
        queues.append(list(queue))
        queue.clear()
        return
    for nodes in tree[node]:
        old_queue = list(queue)
        queue.append(nodes)
        traverse(tree, nodes, queue, queues)
        queue = old_queue


queues = []
traverse(tree, (len(sequence_b), len(sequence_a)), [], queues)
for queue in queues:
    queue.insert(0, (M - 1, N - 1))

# STEP 3: FROM ARROWS TO STRING (DETERMINE GAP POSITIONS)
optimal = []

for queue in queues:
    seq_a = []
    seq_b = []
    ind_a = N - 1
    ind_b = M - 1
    for index in range(len(queue)):
        if index == len(queue) - 1:
            optimal.append([seq_a, seq_b])
            break
        move = tuple(map(operator.sub, queue[index], queue[index + 1]))
        if move == (1, 1):
            seq_a.insert(0, sequence_a[ind_a - 1])
            seq_b.insert(0, sequence_b[ind_b - 1])
            ind_a -= 1
            ind_b -= 1
        if move == (0, 1):
            seq_a.insert(0, sequence_a[ind_a - 1])
            seq_b.insert(0, "-")
            ind_a -= 1
        if move == (1, 0):
            seq_a.insert(0, "-")
            seq_b.insert(0, sequence_b[ind_b - 1])
            ind_b -= 1

# STEP 4: PRINT THE MAP USING MATPLOTLIB


fig, ax = plt.subplots()
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

df = pd.DataFrame(matrix, columns=list(" " + sequence_a))
colors = [['w'] * N for _ in range(M)]

# COLOR THE MAP
x = [j for sub in queues for j in sub]
for i in range(M):
    for j in range(N):
        if (i, j) in x:
            colors[i][j] = 'g'

ax.table(cellText=df.values, colLabels=df.columns,
         loc='center', rowLabels=" " + sequence_b,
         cellColours=colors, colWidths=[0.05 for x in df.columns], )
fig.tight_layout()
plt.show()

# print alignments
print("There are " + str(len(optimal)) + " Optimal Alignments:")
for item in optimal:
    print(' '.join(item[0]))
    print(' '.join(item[1]))
    print("-----------------------------------")
