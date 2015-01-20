/*==================================================================================================
 *  pcabigFn.c
 *
 *  Edited by William Halsey and Scott Rodgers
 *  whalsey@g.clemson.edu
 *  srodger@g.clemson.edu
 *  
 *  This file contains
 *      pcabigFn
 *  
 *  Lasted Edited: Jul. 3, 2013
 *
 *  Changes made: major changes to last section 
 *
 */

 #include "levelTwoOps.c"
 #include "matrix_manipulation.h" 
/*==================================================================================================
 *  pcabigFn
 *
 *  Parameters
 *      double pointer, type double = matrix
 *      value, type integer         = rows
 *      value, type integer         = cols
 *
 *  Returns
 *      N/A
 *      Implicitly returns a matrix through variable 'matrix.'
 *
 *  Description: WHAT NEEDS WORK
 *      LINE 94 - figure out matrix division
 *
 *  THIS FUNCTION CALLS
 *
 *  THIS FUNCTION IS CALLED BY
 *      
 */
void pcabigFn(data_t **U, data_t **R, data_t **E, int rows, int cols, data_t **B, int rows, int cols) {
                                                //  %function [U,R,E] = pcabigFn(B);
                                                //  %Compute PCA by calculating smaller covariance matrix and reconstructing
                                                //  %eigenvectors of large cov matrix by linear combinations of original data
                                                //  %given by the eigenvecs of the smaller cov matrix. 
                                                //  %Data in Cols of B. Third version.  
                                                //  %
                                                //  %***** justification
                                                //  %
                                                //  %B = N x P data matrix.  N = dim of data  (Data in Cols of B, zero mean)
                                                //  %                        P = #examples
                                                //  %                        N >> P
                                                //  %
                                                //  %Want eigenvectors ui of BB' (NxN)
                                                //  %Solution:
                                                //  %Find eigenvectors vi of B'B (PxP)
                                                //  %From def of eigenvector: B'Bvi = di vi ---> BB'Bvi = di Bvi
                                                //  %Eigenvecs of BB' are Bvi
                                                //  %-------------------------------
                                                //  %[V,D] = eig (B'B)
                                                //  %Eigenvecs are in cols of V.    (Sorted cols)
                                                //  %
                                                //  %U = BV;  Cols of U are Bvi (eigenvecs of lg cov mat.) (Gave unit length)
                                                //  %R = B'U; Rows of R are pcarep of each observation.
                                                //  %E = eigenvalues        (eigenvals of small and large cov mats are equal)
                                                //  %*****
                                                //
                                                //  function [U,R,E] = pcabigFn(B)
                                                //
												//  %Read data into columns of B;
    data_t *length_ones;  
	data_t *matrix_t;							//  %B = datamat';
	data_t *temp;
	data_t *temp_vec;
	data_t *index;
	data_t *Vsort;
	data_t *outmatrix;
	data_t *U;
	data_t *righteignm;
	data_t *length_matrix;
	data_t *R;
	//data_t *wr_matrix;
	data_t *wr;
	data_t *wi;
	data_t *lefteigenvect;
	data_t *righteigenvect;
	data_t *workmatrix;
	data_t *vector;
	int integer;
	
	/*  Allocation of all arrays    */
	// Rows by 1 array of ones
	allocate_matrix(&length_ones, rows, 1);
	// temp vars are used to move matrix from on function to another
	allocate_matrix(&temp, rows, cols);
	allocate_matrix(&temp_vec, 1, cols);
	
    allocate_matrix(&index, 1, cols);			//  [N,P] = size(B);    -   sizes found in matlab code                 
	allocate_matrix(&Vsort, cols, cols);
	// image data transposed
	allocate_matrix(&matrix_t, cols, rows);
	allocate_matrix(&outmatrix, cols, cols);
	allocate_matrix(&U, rows, cols);
	allocate_matrix(&righteignm, cols, cols);
	allocate_matrix(&length_matrix, rows, cols);
	allocate_matrix(&R, cols, cols);
	//allocate_matrix(&wr_matrix, cols, cols);
	
	
	// there seems to be a problem here
	allocate_matrix(&wr, cols);
	allocate_matrix(&wi, cols);
	allocate_matrix(&lefteigenvect, cols);
	allocate_matrix(&righteigenvect, cols);
	allocate_matrix(&workmatrix, cols);
	allocate_matrix(&vector, cols);
	/*  END OF ARRAY ALLOCATION */
	
    /*  zero_mean() subtracts out the mean  */  //  %********subtract mean
    zero_mean(matrix, matrix, rows, cols);      //  mb=mean(B');
                                                //  B=B-(ones(P,1)*mb)';
                                                //
												//  %********Find eigenvectors vi of B'B (PxP)
                                                //  [V,D] = eig (1/(P-1)*(B'*B));   %scale factor gives eigvals correct
    transpose(matrix_t, matrix, rows, cols);
	multiply_matrices(outmatrix, matrix_t, matrix, cols, cols, rows);
	divide_by_constant(outmatrix, outmatrix, cols, cols, (data_t)cols - 1);
	matrix_to_vector(vector, outmatrix, rows, cols);
	
	matrix_eig(righteignm, wr?, vector, rows, cols);
	//DGEEV('N', 'V', cols, vector, cols, wr, wi, lefteigenvect, 0, righteigenvect, cols + 1, workmatrix, -1, integer);
    /*  from lapacke    */                      //  %magnitude for large cov mat 
												//  %(assuming sample cov)
                                                //  %********Sort eigenvectors
                                                //  eigvalvec = max(D); -   handled by lapack function
	                                        	//  [seigvals, index] = sort(eigvalvec); % sort goes low to high
	eigSort(righteignm, (data_t*)wr, cols, rows); 
	fliplr(index, index, rows, cols);
	
	fliplr(Vsort, righteignm, rows, cols);
												//Vsort = V(:,fliplr(index));
    /*  !will use mergesort to sort eigenvalues... when two eigenvalues are switched    !
        !the corresponding eigenvectors also will need to be switched                   !   */

                                                //  %********Reconstruct
    multiply_matrices(U, matrix, Vsort, rows, cols, cols);        //  U = B*Vsort;  % Cols of U are Bvi. (N-dim Eigvecs)
                                                //
                                                //  %********Give eigvecs unit length.  Improves classification.
												//  length = sqrt (sum (U.^2)); 
												//temp = sum_along_rows(temp); <<worng?
	raise_matrix_to_power(temp, U, rows, cols, 2);
	sum_columns(temp, temp_vec, rows, cols);
	matrix_sqrt(temp_vec, temp_vec, rows, cols);
    
	ones(length_ones, rows, 1);					//  U = U ./ (ones(N,1) * length);
	multiply_matrices(length_matrix, length_ones, temp_vec, rows, cols, 1);
	matrix_dot_division(U, U, length_ones, rows, cols); 
	multiply_matrices(R, matrix_t, U, cols, cols, rows);      //  R = B'*U;  % Rows of R are pcarep of each image.
    
	// IS wr matrix = sorted eigvals but not set convert wr to matrix and here?
	//fliplr(wr, wr, cols, cols);          //  E = fliplr(seigvals); 
	fliplr(E, wr, cols, cols);
	
	free_matrix(&temp);
	free_matrix(&temp_vec);
	free_matrix(&workmatrix);
    return;
}