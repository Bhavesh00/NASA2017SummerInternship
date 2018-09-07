# This program models the model the self-consistent field method which is essential to the Electronic Structure Theory. - Bhavesh Manivannan

import os
import numpy as np
import scipy.linalg as sp
import sys
import time
import math

np.set_printoptions(threshold=np.nan)

class ProjectFour:
    enuc = 0 # 9.779406144413407
    eelec = 0
    etotal_new = 0
    etotal_old = 99999999999
    one_e = 0
    two_e = 0
    doubly_occupied_orbitals = 5
    overlap_array = [[]]
    overlap_raw_input = []
    l_eig_vector = [[]]
    lambda_eig_value = [[]]
    kinetic_array = [[]]
    kinetic_raw_input = []
    potential_array = [[]]
    potential_raw_input = []
    h_array = [[]]
    eri_array =[[[[]]]]
    eri_raw_input = []
    inverse_square_root = [[]]
    transposed_overlap = [[]]
    fock_array = [[]]
    epsilon_eig_value = [[]]
    c_eig_vector = [[]]
    scf_array = [[]]
    density_array = [[]]
    new_density_array = [[]]
    rms = 0
    convergence = False

    def read_enuc(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.enuc = float(line)

    def create_2d_array(self, row, column):
        n = row
        m = column
        array = [[0] * m for i in range(n)]
        return array

    def read_overlap(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.overlap_raw_input.append(line.strip())

    # Placing the Overlap value inside the 2d list
    def extract_values_from_overlap(self):
        self.overlap_array = self.create_2d_array(24, 24)
        for x in range(3, len(self.overlap_raw_input)):
            # Gets AO
            par = self.overlap_raw_input[x].partition(' ')
            temp_AO = int(par[0])

            # Gets MO
            par = par[2].lstrip().partition(' ')
            temp_MO = int(par[0])

            # Gets Overlap
            par = par[2].lstrip().partition(' ')
            temp_overlap = par[0]

            # Sets the value at index (AO, MO) to the correct value in overlap_array
            self.overlap_array[temp_AO][temp_MO] = temp_overlap

    def read_kinetic(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.kinetic_raw_input.append(line.strip())

    # Placing the Kinetic value inside the 2d list
    def extract_values_from_kinetic(self):
        self.kinetic_array = self.create_2d_array(24, 24)
        for x in range(3, len(self.kinetic_raw_input)):
            # Gets AO
            par = self.kinetic_raw_input[x].partition(' ')
            temp_AO = int(par[0])

            # Gets MO
            par = par[2].lstrip().partition(' ')
            temp_MO = int(par[0])

            # Gets Kinetic
            par = par[2].lstrip().partition(' ')
            temp_kinetic = par[0]

            # Sets the value at index (AO, MO) to the correct value in kinetic_array
            self.kinetic_array[temp_AO][temp_MO] = temp_kinetic

    def read_potential(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.potential_raw_input.append(line.strip())

    # Placing the Potential value inside the 2d list
    def extract_values_from_potential(self):
        self.potential_array = self.create_2d_array(24, 24)
        for x in range(3, len(self.potential_raw_input)):
            # Gets AO
            par = self.potential_raw_input[x].partition(' ')
            temp_AO = int(par[0])

            # Gets MO
            par = par[2].lstrip().partition(' ')
            temp_MO = int(par[0])

            # Gets Potential
            par = par[2].lstrip().partition(' ')
            temp_potential = par[0]

            # Sets the value at index (AO, MO) to the correct value in overlap_array
            self.potential_array[temp_AO][temp_MO] = temp_potential

    def form_h(self):
        self.h_array = self.create_2d_array(24, 24)
        for u in range (0, 24):
            for v in range(0, 24):
                self.h_array[u][v] = float(self.kinetic_array[u][v]) + float(self.potential_array[u][v])

    # Creates the 4d array to store the AOs in
    def create_four_d_array(self, value, *dim):
        """
        Create 4D-array
        :param dim: a tuple of dimensions - (w, x, y, z)
        :param value: value with which 4D-array is to be filled
        :return: 4D-array
        """
        return [[[[value for x in range(dim[3])] for x in range(dim[2])] for x in range(dim[1])] for x in range(dim[0])]

    def read_eri(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.eri_raw_input.append(line.strip())

    def extract_values_from_eri(self):
        m = 24
        n = 24
        o = 24
        p = 24
        
        # Creates the array
        self.eri_array = self.create_four_d_array(0, *(m, n, o, p))
        for x in range(0, len(self.eri_raw_input)):
            # Gets w
            par = self.eri_raw_input[x].partition(' ')
            temp_w = int(par[0])

            # Gets x
            par = par[2].lstrip().partition(' ')
            temp_x = int(par[0])

            # Gets y
            par = par[2].lstrip().partition(' ')
            temp_y = int(par[0])

            # Gets z
            par = par[2].lstrip().partition(' ')
            temp_z = int(par[0])

            # Gets mo_value
            par = par[2].lstrip().partition(' ')
            temp_ao_value = par[0]

            # Set the value in the current iteration to index in 4d array
            self.eri_array[temp_w][temp_x][temp_y][temp_z] = temp_ao_value
            self.eri_array[temp_x][temp_w][temp_y][temp_z] = temp_ao_value
            self.eri_array[temp_x][temp_w][temp_z][temp_y] = temp_ao_value
            self.eri_array[temp_w][temp_x][temp_z][temp_y] = temp_ao_value
            self.eri_array[temp_y][temp_z][temp_w][temp_x] = temp_ao_value
            self.eri_array[temp_y][temp_z][temp_x][temp_w] = temp_ao_value
            self.eri_array[temp_z][temp_y][temp_x][temp_w] = temp_ao_value
            self.eri_array[temp_z][temp_y][temp_w][temp_x] = temp_ao_value

    def convert_lists_to_arrays(self):
        self.overlap_array = np.asarray(self.overlap_array, float)
        self.kinetic_array = np.asarray(self.kinetic_array, float)
        self.potential_array = np.asarray(self.potential_array, float)
        self.h_array = np.asarray(self.h_array, float)
        self.eri_array = np.asarray(self.eri_array, float)

    # Finds Λ (Eigenvector value) and L(Eigenvector matrix)
    def diagonalize_overlap(self):
        #L^−1 (inverse eigenvector matrix)  * S(overlap_array) * L(eigenvector matrix) = Λ (eigenvalue matrix)
        self.lambda_eig_value, self.l_eig_vector, = np.linalg.eigh(self.overlap_array)
        print("Overlap Eigenvalues:", '\n', self.lambda_eig_value, '\n')
        temp_eig_value = self.create_2d_array(24,24)
        temp_eig_value = np.asarray(temp_eig_value, float)
        temp_col = 0
        for row in range (0, 24):
            temp_eig_value[row][temp_col] = self.lambda_eig_value[row]
            temp_col+=1
        self.lambda_eig_value = temp_eig_value

    def double_check_overlap(self):
        temp_transposed_eig_vector = np.transpose(self.l_eig_vector)
        temp_overlap = np.matmul(temp_transposed_eig_vector, self.l_eig_vector)
        arr = np.copy(self.lambda_eig_value)
        self.overlap_array = temp_overlap

    # Finds S^(-1/2)
    def inverse_square_root_matrix(self):
        arr = np.copy(self.overlap_array)
        arr = sp.cholesky(arr)
        arr = np.linalg.inv(arr)
        self.overlap_array = arr
        print("S^(-1/2):", '\n',self.overlap_array)

    # Creates initial guess of the Fock Matrix
    def fock_initial_guess(self):
        temp_array = np.matmul(np.transpose(self.overlap_array), self.h_array)
        temp_array = np.matmul(temp_array, self.overlap_array)
        self.fock_array = temp_array

    # Finds epsilon (Eigenvector value) and C(Eigenvector matrix) of the Fock Matrix
    def diagonalize_fock(self):
        #C^−1 (inverse eigenvector matrix)  * F(fock_array) * C(eigenvector matrix) = epsilon (eigenvalue matrix)
        self.epsilon_eig_value, self.c_eig_vector, = np.linalg.eigh(self.fock_array)
        temp_eig_value = self.create_2d_array(24,24)
        temp_eig_value = np.asarray(temp_eig_value, float)
        temp_col = 0
        for row in range (0, 24):
            temp_eig_value[row][temp_col] = self.epsilon_eig_value[row]
            temp_col+=1
        self.epsilon_eig_value = temp_eig_value

    # C Eigenvector Matrix
    def scf(self):
        self.scf_array = np.matmul(self.overlap_array, self.c_eig_vector)

    # Summation for density matrix (D)
    def form_density_matrix(self):
        self.density_array = self.create_2d_array(24,24)
        self.density_array = np.asarray(self.density_array, float)
        for u in range (0, 24):
            for v in range (0, 24):
                temp_value = 0
                for m in range (0, self.doubly_occupied_orbitals):
                    temp_value += (self.scf_array[u][m] * self.scf_array[v][m])
                self.density_array[u][v] = temp_value

    # Form the new Fock Matrix
    def form_new_fock(self):
        for u in range (0, 24):
            for v in range (0, 24):
                temp_value = 0
                for p in range (0, 24):
                    for o in range(0, 24):
                        temp_value += (self.density_array[p][o] * (2 * self.eri_array[u][v][p][o] - self.eri_array[u][p][v][o]))
                self.fock_array[u][v] = temp_value + self.h_array[u][v]


    # Computes the Electronic energy (eelec)
    def compute_electronic_energy(self):
        temp_value = 0
        e_one = 0
        e_two = 0
        for u in range (0, 24):
            for v in range (0, 24):
                temp_value += (self.density_array[u][v] * (self.h_array[u][v] + self.fock_array[u][v]))
                e_one += self.density_array[u][v] * self.h_array[u][v]
                e_two += self.density_array[u][v] * self.fock_array[u][v]
        self.eelec = temp_value
        print("Electronic Energy:", self.eelec)
        print("Total Energy:", self.eelec + self.enuc)
        self.one_e = e_one
        self.two_e = e_two
        print("E One(Density Matrix Value * H Value):", e_one)
        print("E Two(Density Matrix Value * Fock Matrix Value):", e_two)

    # Adds the Nuclear energy and the Electrical energy together
    def compute_total_energy(self):
        self.etotal_new = self.eelec + self.enuc

    def transform_new_fock(self):
        temp_array = np.matmul(np.transpose(self.overlap_array), self.fock_array)
        temp_array = np.matmul(temp_array, self.overlap_array)
        self.fock_array= temp_array

    def diagonalize_transformed_fock(self):
        # C^−1 (inverse eigenvector matrix)  * F(fock_array) * C(eigenvector matrix) = epsilon (eigenvalue matrix)
        self.epsilon_eig_value, self.c_eig_vector, = np.linalg.eigh(self.fock_array)
        temp_eig_value = self.create_2d_array(24,24)
        temp_eig_value = np.asarray(temp_eig_value, float)
        temp_col = 0
        for row in range (0, 24):
            temp_eig_value[row][temp_col] = self.epsilon_eig_value[row]
            temp_col+=1

    def scf_transform(self):
        self.scf_array = np.matmul(self.overlap_array, self.c_eig_vector)

    def form_new_density_matrix(self):
        self.new_density_array = self.create_2d_array(24, 24)
        self.new_density_array = np.asarray(self.new_density_array, float)
        for u in range(0, 24):
            for v in range(0, 24):
                temp_value = 0
                for m in range(0, self.doubly_occupied_orbitals):
                    temp_value += (self.scf_array[u][m] * self.scf_array[v][m])
                self.new_density_array[u][v] = temp_value

    # returns true if passes test and false if it doesnt (doesnt check difference in total energy)
    def test_convergence(self):
        temp_value = 0
        for u in range (0, 24):
            for v in range (0, 24):
                temp_value += ((self.new_density_array[u][v] - self.density_array[u][v])**2)
        self.rms = math.sqrt(temp_value)
        if(self.rms<10**-10):
            return True
        else:
            return False

    def difference_energy(self):
        if((self.etotal_new - self.etotal_old) < 10**-10):
            return True
        else:
            return False

    # Checks difference in total energy + iterates up to 10
    def iteration(self):
        # if test_convergence returns false then set the old electronic energy to the new electronic energy and do the TEST
        for i in range (0, 50):
            print(i, '\n')
            test.form_new_fock()
            test.compute_electronic_energy()
            test.compute_total_energy()
            if self.etotal_old != 99999999999:
                if self.convergence and self.difference_energy():
                    print("Final RMS:", self.rms)
                    print("Final Difference:", self.etotal_new-self.etotal_old)
                    print("Final Energy:", self.two_e + self.one_e + self.enuc)
                    print ("Final Iterations", i)
                    break
            else:
                self.etotal_old = self.etotal_new
            test.transform_new_fock()
            test.diagonalize_transformed_fock()
            test.scf_transform()
            test.form_new_density_matrix()
            self.convergence = test.test_convergence()
            self.density_array = self.new_density_array

start = time.time()
test = ProjectFour()
test.read_enuc("enuc.txt")
test.read_overlap("ao_overlap.txt")
test.extract_values_from_overlap()
test.read_kinetic("ao_kinetic.txt")
test.extract_values_from_kinetic()
test.read_potential("ao_potential.txt")
test.extract_values_from_potential()
test.form_h()
test.read_eri("ao_eri.txt")
test.extract_values_from_eri()
test.convert_lists_to_arrays()
test.diagonalize_overlap()
test.inverse_square_root_matrix()
test.fock_initial_guess()
test.diagonalize_fock()
test.scf()
test.form_density_matrix()
test.iteration()
