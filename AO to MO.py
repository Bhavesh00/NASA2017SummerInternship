# This program transforms two-electron integrals from the atomic orbital basis to the molecular orbital basis.
# - Bhavesh Manivannan

import time

class ProjectTwo:
    c_array = [[]]
    c_raw_input = []
    ao = [[[[]]]]
    ao_raw_input = []
    mo = [[[[]]]]
    mo_raw_input = []
    final_answer = [[[[]]]]
    first_loop_time = 0
    second_loop_time = 0
    third_loop_time = 0
    fourth_loop_time = 0

    def create_2d_array(self):
        n = 58
        m = 58
        global c_array
        c_array = [[0] * m for i in range(n)]

    def read_c(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.c_raw_input.append(line.strip())

    # Placing the C value inside the 2d list
    def extract_values_from_c(self):
        global c_array
        for x in range(3, len(self.c_raw_input)):
            # Gets AO
            par = self.c_raw_input[x].partition(' ')
            temp_AO = int(par[0])

            # Gets MO
            par = par[2].lstrip().partition(' ')
            temp_MO = int(par[0])

            # Gets C
            par = par[2].lstrip().partition(' ')
            temp_C = par[0]

            # Sets the value at index (AO, MO) to the correct value in C
            c_array[temp_AO][temp_MO] = temp_C

    # Testing to see if print works
    def print_C_raw_input(self, indexOne):
        print(self.c_raw_input[indexOne])

    # Creates the 4d array to store the AOs in
    def create_four_d_array(self, value, *dim):
        """
        Create 4D-array
        :param dim: a tuple of dimensions - (w, x, y, z)
        :param value: value with which 4D-array is to be filled
        :return: 4D-array
        """
        return [[[[value for x in range(dim[3])] for x in range(dim[2])] for x in range(dim[1])] for x in range(dim[0])]

    # Reads the AO data into a list
    def read_ao(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.ao_raw_input.append(line.strip())

    def extract_values_from_ao(self):
        global ao
        m = 58
        n = 58
        o = 58
        p = 58

        # Creates the array
        ao = self.create_four_d_array(0, *(m, n, o, p))

        for x in range(0, len(self.ao_raw_input)):
            # Gets a
            par = self.ao_raw_input[x].partition(' ')
            temp_a = int(par[0])

            # Gets b
            par = par[2].lstrip().partition(' ')
            temp_b = int(par[0])

            # Gets c
            par = par[2].lstrip().partition(' ')
            temp_c = int(par[0])

            # Gets d
            par = par[2].lstrip().partition(' ')
            temp_d = int(par[0])

            # Gets ao_value
            par = par[2].lstrip().partition(' ')
            temp_ao_value = par[0]

            # Set the value in the current iteration to index in 4d array
            ao[temp_a][temp_b][temp_c][temp_d] = temp_ao_value
            ao[temp_b][temp_a][temp_c][temp_d] = temp_ao_value
            ao[temp_b][temp_a][temp_d][temp_c] = temp_ao_value
            ao[temp_a][temp_b][temp_d][temp_c] = temp_ao_value
            ao[temp_c][temp_d][temp_a][temp_b] = temp_ao_value
            ao[temp_c][temp_d][temp_b][temp_a] = temp_ao_value
            ao[temp_d][temp_c][temp_b][temp_a] = temp_ao_value
            ao[temp_d][temp_c][temp_a][temp_b] = temp_ao_value

    def print_ao_value_at_index(self, index_a, index_b, index_c, index_d):
        print(ao[index_a][index_b][index_c][index_d])

    # Testing function created to see value in raw input
    def print_ao_raw_input(self, index):
        print(self.ao_raw_input[index])

    # Reads the MO data into a list
    def read_mo(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.mo_raw_input.append(line.strip())

    def extract_values_from_mo(self):
        m = 58
        n = 58
        o = 58
        p = 58

        # Creates the array
        global mo
        mo = self.create_four_d_array(0, *(m, n, o, p))

        for x in range(0, len(self.mo_raw_input)):
            # Gets i
            par = self.mo_raw_input[x].partition(' ')
            temp_i = int(par[0])

            # Gets j
            par = par[2].lstrip().partition(' ')
            temp_j = int(par[0])

            # Gets k
            par = par[2].lstrip().partition(' ')
            temp_k = int(par[0])

            # Gets l
            par = par[2].lstrip().partition(' ')
            temp_l = int(par[0])

            # Gets mo_value
            par = par[2].lstrip().partition(' ')
            temp_mo_value = par[0]

            # Set the value in the current iteration to index in 4d array
            mo[temp_i][temp_j][temp_k][temp_l] = temp_mo_value

    def print_mo_value_at_index(self, index_i, index_j, index_k, index_l):
        print(mo[index_i][index_j][index_k][index_l])

    # Testing function created to see value in raw input
    def print_mo_raw_input(self, index):
        print(self.mo_raw_input[index])

    # Truncates floats
    def truncate_float(self, x, y):
        # Truncates/pads a float f to n decimal places without rounding
        s = '{}'.format(x)
        if 'e' in s or 'E' in s:
            return '{0:.{1}f}'.format(x, y)
        i, p, d = s.partition('.')
        return '.'.join([i, (d + '0' * y)[:y]])

    # Create nested for loops that go through all of the indexes (a, b, c, d, i, j, k, l) from 0 to 57
    def summation(self):
        global final_answer
        output_name = 'output.dat'
        output = open(output_name, 'w')
        output.truncate()
        output.close()
        output = open(output_name, 'w')

        for mo in range(0, 2):
            sum_a = 0
            sum_b = 0
            sum_c = 0
            sum_d = 0
            for a in range(0, 58):
                for b in range(0, 58):
                    for c in range (0, 58):
                        for d in range (0, 58):
                            sum_d += (float(c_array[d][mo]) * float(ao[a][b][c][d]))
                            sum_c += (float(c_array[c][mo]) * sum_d)
                            sum_d = 0
                    sum_b += (float(c_array[b][mo]) * sum_c)
                    sum_c = 0
                sum_a += (float(c_array[a][mo]) * sum_b)
                sum_b = 0

            spacer = "    "
            final_output = self.truncate_float(sum_a, 10)
            wrote = spacer.join([str(int(mo)), str(int(mo)), str(int(mo)), str(int(mo)), str(final_output), '\n'])
            output.write(wrote)
        output.close()

    def summationTwo(self):
        output_name = 'output.dat'
        output = open(output_name, 'w')
        output.truncate()
        output.close()
        output = open(output_name, 'w')
        mo_value_limit = 3

        for i in range(0, mo_value_limit):
            for j in range(0, mo_value_limit):
                for k in range(0, mo_value_limit):
                    for l in range(0, mo_value_limit):
                        num_four = 0
                        for a in range(0, 58):
                            num_three = 0
                            for b in range(0, 58):
                                num_two = 0
                                for c in range(0, 58):
                                    num_one = 0
                                    for d in range(0, 58):
                                        num_one = float(num_one) + (float(c_array[d][l]) * float(ao[a][b][c][d]))
                                    num_two = float(num_two) + (float(c_array[c][k]) * float(num_one))
                                num_three = float(num_three) + (float(c_array[b][j]) * float(num_two))
                            num_four += (float(c_array[a][i]) * float(num_three))

                        x = "    "
                        final_output = self.truncate_float(num_four, 10)
                        y = x.join([str((i)), str(int(j)), str(int(k)), str(int(l)), str(final_output), '\n'])
                        output.write(y)
        output.close()

    def summation_three(self):
        """
                i, a, b, c, d
                j, i, b, c, d
                k, i, j, c, d
                l, i, j, k, d
        """
        global final_answer
        global first_loop_time
        global second_loop_time
        global third_loop_time
        global fourth_loop_time

        intermediate_array_one = self.create_four_d_array(0, *(58, 58, 58, 58))
        intermediate_array_two = self.create_four_d_array(0, *(58, 58, 58, 58))
        intermediate_array_three = self.create_four_d_array(0, *(58, 58, 58, 58))
        final_answer = self.create_four_d_array(0, *(58, 58, 58, 58))
        output_name = 'output.dat'
        output = open(output_name, 'w')
        output.truncate()
        output.close()
        output = open(output_name, 'w')
        number_output = 3
        start= time.time()
        for i in range (0, number_output):
            for a in range (0, 58):
                for b in range (0,58):
                    for c in range(0, 58):
                        for d in range (0,58):
                            intermediate_array_one[i][b][c][d] = float(intermediate_array_one[i][b][c][d]) + (float(c_array[a][i]) * float(ao[a][b][c][d]))
        end = time.time()
        first_loop_time = end-start

        start = time.time()
        for j in range (0,number_output):
            for i in range (0, number_output): # This should be (0, number_output)
                for b in range (0,58):
                    for c in range (0, 58):
                        for d in range (0,58):
                            intermediate_array_two[i][j][c][d] = float(intermediate_array_two[i][j][c][d]) + (float(c_array[b][j]) * float(intermediate_array_one[i][b][c][d]))
        end = time.time()
        second_loop_time = end-start

        start = time.time()
        for k in range (0,number_output):
            for i in range (0, number_output):
                for j in range (0, number_output):
                    for c in range(0, 58):
                        for d in range (0, 58):
                            intermediate_array_three[i][j][k][d] = float(intermediate_array_three[i][j][k][d]) + (float(c_array[c][k]) * float(intermediate_array_two[i][j][c][d]))
        end = time.time()
        third_loop_time = end-start

        start = time.time()
        for l in range (0, number_output):
            for i in range (0, number_output):
                for j in range (0, number_output):
                    for k in range (0, number_output):
                        for d in range(0, 58):
                            final_answer[i][j][k][l] = float(final_answer[i][j][k][l]) + (float(c_array[d][l]) * float(intermediate_array_three [i][j][k][d]))
                        final_output = self.truncate_float(final_answer[i][j][k][l], 10)
                        if((float(final_output) != 0.0000000000 and float(final_output)!= -0.0000000000) and final_answer[j][i][k][l] != 0
                            and final_answer[j][i][l][k] != 0 and final_answer[i][j][l][k] != 0):
                            g = "    "
                            f = g.join([str((i)), str(int(j)), str(int(k)), str(int(l)), str(final_output), '\n'])
                            output.write(f)
        end = time.time()
        fourth_loop_time = end-start
        output.close()

start = time.time()
test = ProjectTwo()
test.read_c('C.dat')
test.create_2d_array()
test.extract_values_from_c()
test.read_ao('ao.dat')
test.extract_values_from_ao()
test.summationThree()
end = time.time()
print("Total Time: " + str(end-start))
