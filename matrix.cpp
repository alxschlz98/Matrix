#include <vector>
#include <iostream>
#include <iomanip>

class matrix {
public:
    int m, n;
    std::vector<std::vector<double> > matr;

    matrix(int rowL, int columnL, std::vector<double> a) {
        //Constructs a matrix of size [rowL, columnL]. 
        //Elements are defined in a.
        m = rowL;
        n = columnL;
        std::vector<double> blankVec;
        if (m && n) {
            matr.push_back(blankVec);
        }

        int i = 0, im = 0, in = 0;
        while (i < rowL * columnL) {
            //while loop pushes the ith element in a if it exists, and 0 otherwise
            if (i < a.size()) {
                matr[im].push_back(a[i]);
            }
            else {
                matr[im].push_back(0);
            }

            if (in == columnL - 1) {
            //we'll need to push another blank vector to the end of the vector vector
            //before we can add more elements to it
                in = 0; im++;
                matr.push_back(blankVec);
            }

            else {
                in++;
            }
            i++;
        }
    }

    matrix(int rowL, int columnL) {
        //constructs a matrix of size rowL, columnL.
        m = rowL;
        n = columnL;

        std::vector<double> blankVec(n,0);
        std::vector<std::vector<double> > tempVec(m, blankVec);
        matr = tempVec;
    }

    void print() {
        //Prints the matrix to stdout.
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << std::fixed << std::setprecision(3) << matr[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    matrix operator+(matrix& addend) {
        //Adds two matrices together. If they aren't equal size, return the 0 matrix
        if (m != addend.m || n != addend.n)
            return matrix(0,0);

        matrix sum(m,n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                sum.matr[i][j] = matr[i][j] + addend.matr[i][j];
            }
        }
        return sum;
    }

    matrix operator*(double c) {
        //Multiplies a matrix by a scalar.
        matrix newMatr(m,n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                newMatr.matr[i][j] = matr[i][j]*c;
            }
        }
        return newMatr;
    }

    matrix operator/(double c) {
        //Divides a matrix by a scalar
        return (*this) * (1/c);
    }

    matrix operator*(matrix& factor) {
        //Multiplies two matricies together
        
        //Return a matrix with size 0 if the matricies aren't compatible.
        if (n != factor.m)
            return matrix(0,0);
        
        matrix product(m, factor.n);

        for (int row = 0; row < m; row++) {
            for (int column = 0; column < factor.n; column++) {
                for (int i = 0; i < n; i++) {
                    product.matr[row][column] += matr[row][i]*factor.matr[i][column];
                }
            }
        }
        return product;
    }

    matrix cutRow(int cm, int cn) {
        //Removes a row and a column from a matrix.
        int i = 0, j = 0, nI = 0, nJ = 0;
        matrix newMatr(m - 1, n - 1);

        //Go through each element in the matrix. If it's in row cm or column cn, skip it.
        while (i < m) {
            if (i == cm) { i++; continue; }
            j = 0; nJ = 0;
            while (j < n) {
                if (j == cn) { j++; continue; }
                newMatr.matr[nI][nJ] = matr[i][j];
                j++; nJ++;
            }
            i++; nI++;
        }
        return newMatr;
    }

    double minor(int cm, int cn) {
        //Returns the determinant of the matrix excluding rows cm and cn
        matrix newMatr = cutRow(cm, cn);
        return newMatr.determinant();
    }

    double cofactor(int cm, int cn) {
        //Returns the cofactor of [cm, cn]
        //Sign determined by (-1)^(cm+cn).
        //Note that this will be -1 if cm + cn is odd, and 1 otherwise.
        if ((cm + cn) % 2) {
            return minor(cm, cn) * -1;
        }
        return minor(cm, cn);    
    }

    double determinant() {
        //Returns the determinant of the matrix. If it isn't a square matrix, return 0.
        //Note that the determinant of a 1x1 matrix is the element of the matrix.
        if (m != n) { return 0; }
        if (m == 1) { return matr[0][0]; }

        double sum = 0;
        for (int i = 0; i < n; i++) {
            //Calculate the determinant using the first row:
            sum += cofactor(0, i) * matr[0][i];
        }
        return sum;
    }

    void multRow(int row, double scalar) {
        //Multiply a row by a scalar
        for (int i = 0; i < n; i++) {
            matr[row][i] *= scalar;
        }
    }

    void addRowToRow(int row1, int row2, double scalar) {
        //Add a scalar multiple of a row to another row
        double tempVal;
        for (int i = 0; i < n; i++) {
            matr[row2][i] += matr[row1][i] * scalar;
        }
    }

    void swapRow(int row1, int row2) {
        //Swap two rows
        std::vector<double> tempRow = matr[row1];
        matr[row1] = matr[row2];
        matr[row2] = tempRow;
    }

    matrix rref() {
        //Calculates the reduced row echelon form of the matrix
        //Algorithm pseudocode from https://rosettacode.org/wiki/Reduced_row_echelon_form
        matrix newMatrix = *this;
        int lead = 0, i;

        for (int r = 0; r < m; r++) {
            if (n <= lead)
                return newMatrix;

            i = r;
            while (newMatrix.matr[i][lead] == 0) {
                std::cout << "loop" << std::endl;
                i++;
                if (m == i) {
                    i = r;
                    lead++;

                    if (n == lead)
                        return newMatrix;
                }
            }
            
            newMatrix.swapRow(i,r);

            if (newMatrix.matr[r][lead] != 0)
                newMatrix.multRow(r, 1 / newMatrix.matr[r][lead]);

            for (i = 0; i < m; i++) {
                if (i != r) {
                    newMatrix.addRowToRow(r, i, -newMatrix.matr[i][lead]);
                }
            }
            lead++;
        }
        return newMatrix;
    }

    matrix T() {
        //Returns the transpose of the matrix
        matrix newMatr(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = n-1; j >= 0; j--) {
                newMatr.matr[j][i] = matr[i][j];
            }
        }
        return newMatr;
    }

    matrix inverse() {
        //Returns the inverse of the matrix
        //Algorithm used: Get the cofactor matrix, transpose it, then divide it by the determinant
        if (m != n) { return matrix(0,0); }
        matrix newMatr(m,n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                newMatr.matr[i][j] = cofactor(i,j);
            }
        }
        newMatr = newMatr.T();
        newMatr = newMatr / determinant();

        return newMatr;
    }
};
