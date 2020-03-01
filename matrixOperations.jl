#=
# filename: matrixOperations.jl
# Created by Spartan City
#   Nick Guthrie, Nick Filerman, Brendan Walsh, Michael Moran
# 
# A library to provide several operations for matrices, such
# as add, subtract, multiply, normalize, inverse, and more!
=#

#=
# A function to add two same size matrices together.
# \param matrix1 Is the first matrix in the operation
# \param matrix2 Is the second matrix in the operation
# \returns The sum of the two matrices.
=#
function add(matrix1, matrix2)
    if size(matrix1) == size(matrix2)
        matrix3 = Matrix{Float64}(undef, size(matrix1))
        for i = 1:length(matrix1)
            matrix3[i] = matrix1[i] + matrix2[i]
        end
        return matrix3

    else
        throw(DomainError("Dimensions of matrices must be equal"))
    end
end



#=
# Function to subtract two matrices.
# \param mat1 The first matrix to use in the operation
# \param mat2 The second matrix to use in the operation
# \returns A matrix resulting from the subtraction
=#
function subtract(mat1, mat2)
    if size(mat1) == size(mat2)
        mat3 = Matrix{Float64}(undef, size(mat1))
        for i = 1:length(mat1)
            mat3[i] = mat1[i] - mat2[i]
        end
        return mat3
    else
        throw(DomainError("Dimensions of matrices must be equal"))
    
    end
end



#=
# Function to multiply a matrix by a scalar
# \param matr The matrix to use in the operation
# \param factor Number to multiply the matrix by
# \returns A matrix that has been multiplied.
=#
function scalarMultiplication(matr, factor)
    n = length(matr)
    returnMatrix = Matrix{Float64}(undef, size(matr)[1], size(matr)[2])
    for i = 1:n
        if matr[i] != undef
            returnMatrix[i] = matr[i] * factor
        else
            throw(DomainError("Index undefined, cannot multiply"))
    
        end
    end
    return float(returnMatrix)
end



#=
# Returns the minor of a given matrix
# \param matrix Is the matrix to take the minor from
# \param row The row to exclude from the matrix
# \param column The column to exclude from the matrix
# \returns A sub matrix smaller than the given matrix
=#
function minor(matrix, row, column)
    if size(matrix)[1] == 1 || size(matrix)[2] == 1
        throw(DomainError("Cannot create a submatrix with a 1xN or Nx1 matrix"))
    end
    sub_matrix = Matrix{Float64}(undef, size(matrix)[1]-1, size(matrix)[2]-1)
    for j in 1:size(matrix)[2]-1
        if j >= row
            for k in 1:size(matrix)[1]-1
                if k >= column
                    sub_matrix[j, k] = matrix[j+1, k+1]
                else
                    sub_matrix[j, k] = matrix[j+1, k]
                end
            end
        else
            for k in 1:size(matrix)[1]-1
                if k >= column
                    sub_matrix[j, k] = matrix[j, k+1]
                else
                    sub_matrix[j, k] = matrix[j, k]
                end
            end
        end
    end
    return sub_matrix
end



#=
# Function to multiply two matrices together.
# \param matrix1 Is the first matrix of the operation
# \param matrix2 Is the second matrix of the operation
# \returns The product of the two matrices.
=#
function multiply(matrix1, matrix2)
    if size(matrix1)[1] == 1 || size(matrix2)[1] == 1 || size(matrix1)[2] == 1 || size(matrix2)[2] == 1
        throw(DomainError("Neither matrices may be vectors"))

    elseif size(matrix1)[2] == size(matrix2)[1]
        matrix3 = Matrix{Float64}(undef, size(matrix1)[1], size(matrix2)[2])
        vector1 = Matrix{Float64}(undef, 1, size(matrix1)[2])
        vector2 = Matrix{Float64}(undef, 1, size(matrix1)[2])

        n = 1
        sum = 0
        for i = 1:size(matrix1)[1]
            m = 1
            for j = 1:size(matrix1)[2]
                vector1[1, j] = matrix1[i, j]
            end

            for j = 1:size(matrix2)[2]
                for k = 1:size(matrix2)[1]
                    vector2[1, k] = matrix2[k, j]
                    sum = 0

                    for l = 1:size(matrix2)[1]
                        sum +=  vector1[1, l] * vector2[1, l]
                    end
                    
                    matrix3[n, m] = sum
                end
                m += 1
            end
            n += 1
        end

        return matrix3

    else
        throw(DomainError("Inner dimensions of matrices must be equal"))
    
    end

end



#=
# This function inverts a given matrix.
# \param matrix Is the matrix that the program will invert
# \returns An inverted matrix
=#
function inverse(matrix)

    scalar = determinant(matrix)
    if scalar == 0
        throw(DomainError("Matrix does not have an inverse"))
    end

    inverse = Matrix{Float64}(undef, size(matrix)[1], size(matrix)[2])
    scalar = 1.0/scalar

    for i in 1:size(matrix)[2]
        for j in 1:size(matrix)[1]
            inverse[i, j] = determinant(minor(matrix, i, j)) * scalar
        end
    end

    return inverse
end



#=
# Gives the product of a matrix and another matrix's inverse
# 
=#
function divide(matrix1, matrix2)
    return multiply(matrix1, inverse(matrix2))
end



#=
# A function to normalize a given matrix
# \param matrix Is the matrix to normalize
# \returns Nothing. This function edits the given matrix directly
=#
function normalize(original)
    matrix = original
    scalar = determinant(matrix)
    if scalar == 0
        throw(DomainError("Matrix does not have an inverse"))
    end
    scalar = 1.0 / scalar

    for i in 1:size(matrix)[i]
        for j in 1:size(matrix)[i]
            matrix[i, j] *= scalar
        end
    end
    
    return matrix
end



#=
# Function to find the determinant of a matrix.
# \param matrix Is the matrix to get the determinant of
# \returns The value of the determinant
=#
function determinant(matrix)
    if length(size(matrix)) == 1
        return matrix[1, 1]

    elseif size(matrix)[2] == size(matrix)[1] == 2
        returnScalar = 0
        returnScalar += matrix[1, 1] * matrix[2, 2] - matrix[2, 1] * matrix[1, 2]

        return returnScalar

    elseif size(matrix)[2] == size(matrix)[1] == 3
        returnScalar = 0
        returnScalar += matrix[1, 1] * (matrix[2, 2] * matrix[3, 3] - matrix[2, 3] * matrix[3, 2])
        returnScalar += matrix[1, 2] * (matrix[2, 3] * matrix[3, 1] - matrix[2, 1] * matrix[3, 3])
        returnScalar += matrix[1, 3] * (matrix[2, 1] * matrix[3, 2] - matrix[2, 2] * matrix[3, 1])

        return returnScalar

    elseif size(matrix)[2] == size(matrix)[1]
        returnScalar = 0
        for i in 1:size(matrix)[2]
            sub_matrix = minor(matrix, 1, i)

            total = matrix[1, i] * determinant(sub_matrix)
            if i % 2
                returnSchalar += total
            else
                returnSchalar -= total
            end
        end

        return returnScalar
        
    else
        throw(DomainError("Matrix must be square of size NxN"))
    end
    
end