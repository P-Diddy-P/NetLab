import numpy as np


class NestednessCalculator(object):
    """Calculates the nestedness of the input matrix.
    The algorithms that have been implemented are:
        - NODF (Nestedness based on Overlap and Decreasing Fill)
    """
    def __init__(self, mat):
        """Initialize the Nestedness calculator and check the input matrix.
        :param mat: binary input matrix
        :type mat: numpy.array
        """
        self.check_input_matrix_is_binary(mat)
        self.check_degrees(mat)

    @staticmethod
    def check_input_matrix_is_binary(mat):
        """Check that the input matrix is binary, i.e. entries are 0 or 1.
        :param mat: binary input matrix
        :type mat: numpy.array
        :raise AssertionError: raise an error if the input matrix is not
            binary
        """
        assert np.all(np.logical_or(mat == 0, mat == 1)), \
            "Input matrix is not binary."

    @staticmethod
    def check_degrees(mat):
        """Check that rows and columns are not completely zero.
        :param mat: binary input matrix
        :type mat: numpy.array
        :raise AssertionError: raise an error if the input matrix has
            completely zero rows or columns.
        """
        assert np.all(mat.sum(axis=1) != 0), \
            "Input matrix rows with only zeros, abort."
        assert np.all(mat.sum(axis=0) != 0), \
            "Input matrix columns with only zeros, abort."

################################################################################
# NODF - Nestedness based on Overlap and Decreasing Fill
################################################################################

    @staticmethod
    def get_paired_nestedness(mat, rows=True):
        """Calculate the paired nestedness along the rows or columns of the.
        :param mat: binary input matrix
        :type mat: numpy.array
        :param rows: if True, pairs are calculated along the rows, if False
            along the columns
        :type rows: bool
        :returns: degree of paired nestedness
        :rtype: float
        The method uses the algorithm described in the `BiMat framework for
        MATLAB <https://bimat.github.io/alg/nestedness.html>`_.
        """
        if rows:
            # consider rows
            po_mat = np.dot(mat, mat.T)
            degrees = mat.sum(axis=1)
        else:
            # consider cols
            po_mat = np.dot(mat.T, mat)
            degrees = mat.sum(axis=0)
        assert len(degrees) == len(po_mat)

        neg_delta = (degrees != degrees[:, np.newaxis])
        deg_matrix = degrees * np.ones_like(po_mat)
        deg_minima = np.minimum(deg_matrix, deg_matrix.T)
        n_pairs = po_mat[neg_delta] / (2. * deg_minima[neg_delta])
        return n_pairs.sum()

    @classmethod
    def nodf(cls, mat):
        """Calculate the NODF nestedness of the input matrix [AlmeidaNeto]_.
        :param mat: binary input matrix
        :type mat: numpy.array
        :returns: NODF nestedness of the input matrix
        :rtype: float
        The algorithm has been tested by comparison with the `online tool
        provided at <http://ecosoft.alwaysdata.net/>`_
        """
        n_pairs_rows = cls.get_paired_nestedness(mat, rows=True)
        n_pairs_cols = cls.get_paired_nestedness(mat, rows=False)
        print(f"row n_pairs: {n_pairs_rows}")
        print(f"col n_pairs: {n_pairs_cols}")
        norm = np.sum(np.array(mat.shape) * (np.array(mat.shape) - 1) / 2.)
        nodf = (n_pairs_rows + n_pairs_cols) / norm
        return nodf


def nodf_basic(network_table):
    """
    Basic implementation of NODF algorithm without linear algebra black magic
    as a verification implementation. Based on maxnodf package paper:
    https://www.biorxiv.org/content/10.1101/2020.03.20.000612v1.full.pdf
    :param network_table:
    :return:
    """
    zo_table = network_table
    nrows, ncols = zo_table.shape
    row_degrees, col_degrees = zo_table.sum(axis=1), zo_table.sum(axis=0)
    nodfr, nodfc = [], []

    for r1 in range(nrows):
        row_contribution = 0
        for r2 in range(r1+1, nrows):
            r1_deg, r2_deg = row_degrees[r1], row_degrees[r2]
            if r1_deg > r2_deg:
                r12_overlap = sum(zo_table[r1][j] * zo_table[r2][j] for j in range(ncols))
                row_contribution += r12_overlap / r2_deg

        nodfr.append(row_contribution)

    for c1 in range(ncols):
        col_contribution = 0
        for c2 in range(c1+1, ncols):
            c1_deg, c2_deg = col_degrees[c1], col_degrees[c2]
            if c1_deg > c2_deg:
                c12_overlap = sum(zo_table[i][c1] * zo_table[i][c2] for i in range(nrows))
                col_contribution += c12_overlap / c2_deg

        nodfc.append(col_contribution)

    print(nodfr, nodfc)




def network_nodf_explained(network_table):
    """
    Calculate network NODF index for the entire network, as well as
    per species NODF contributions.

    Base algorithm: http://guimaraes.bio.br/032.pdf
    Reference implementation: https://github.com/tsakim/nestedness/blob/master/nestedness_calculator.py
    :param network_table: a pandas dataframe of the network
    :return: A 2-tuple of (network NODF score, species relative NODF dictionary)
    """
    zo_table = network_table

    # Calculate row N_paired
    # Calculate total overlap (not divided by lower row degree yet)
    row_overlap = np.dot(zo_table, zo_table.T)

    # Calculate a matrix where each row contains the row degrees of zo_table
    row_degrees = zo_table.sum(axis=1)
    row_deg_matrix = row_degrees * np.ones_like(row_overlap)

    # Calculate a binary matrix M where M_ij=1 if row i has higher degree (marginal sum)
    # than row j
    df_paired_row = (row_degrees > row_degrees[:, np.newaxis])

    # hadamard product of df_paired matrix and a lower triangular matrix will give a binary
    # matrix where only entries which should be calculated have value of 1, i.e. entries M_ij
    # where i < j (i.e. for each row consider lower rows) and deg(i) > deg(j)
    valid_matrix = np.multiply(df_paired_row, np.tril(np.ones_like(row_overlap), k=-1))

    # partial overlap matrix is an entry-wise division of each overlap entry with the degree of
    # the lower row. This calculates correct partial overlaps only for the upper triangle, but the
    # rest doesn't matter
    row_partial_overlap = np.divide(row_overlap, row_deg_matrix)

    # finally, multiply partial_overlap with valid matrix, the trace will contain NODF value donated
    # by each row, this can be summed to get a total NODF value and divided to get individual NODF
    # contributions per species
    row_nodf = np.dot(row_partial_overlap, valid_matrix).diagonal()

    print(f"explained row n_pairs contribution: {row_nodf}")

    # Calculate column N_paired
    # calculate total column overlap by transposing first
    col_overlap = np.dot(zo_table.T, zo_table)

    # Calculate a matrix where each row contains the column degrees of zo_tablegfd
    col_degrees = zo_table.sum(axis=0)
    col_degrees_matrix = col_degrees * np.ones_like(col_overlap)

    # Calculate a binary matrix with 1 in M_ij if the marginal total of col i is strictly larger
    # than that of col j. There's no assumption on whether i>j or j<=i yet, but note that entries
    # below the diagonal are where i>j
    df_paired_col = (col_degrees > col_degrees[:, np.newaxis])
    valid_matrix = np.multiply(df_paired_col, np.tril(np.ones_like(col_overlap), k=-1))

    col_partial_overlap = np.divide(col_overlap, col_degrees_matrix)
    col_nodf = np.dot(col_partial_overlap, valid_matrix).diagonal()
    print(f"explained column n_pairs contribution: {col_nodf}")


if __name__ == "__main__":
    base_matrix = np.array([[1, 1, 0, 1], [0, 1, 1, 0], [0, 0, 0, 1]])

    print(f"base matrix:\n{base_matrix}")
    NestednessCalculator.nodf(base_matrix)
    print()
    network_nodf_explained(base_matrix)
    print()
    nodf_basic(base_matrix)
