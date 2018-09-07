# This program calculates the second Møller–Plesset perturbation theory energy of water. - Bhavesh Manivannan

import time
class ProjectThree:
    f_array = []
    f_raw_input = []
    answer_array = [[[[]]]]
    answer_raw_input = []

    def create_1d_array(self):
        n = 58
        global f_array
        f_array = [0] * n

    def read_f(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as f_obj:
            lines = f_obj.readlines()
        for line in lines:
            self.f_raw_input.append(line.strip())

    # Placing the F value inside the 1d list
    def extract_values_from_f(self):
        global f_array
        for x in range(0, len(self.f_raw_input)):

            # Gets Index
            par = self.f_raw_input[x].partition(' ')
            temp_index = int(par[0])

            # Gets Energy
            par = par[2].lstrip().partition(' ')
            temp_energy = float(par[0])

            # Sets the value at temp_index to the correct value in F
            f_array[temp_index] = temp_energy

    def read_answer(self, __file_name__):
        file_name = __file_name__
        with open(file_name) as answer_obj:
            lines = answer_obj.readlines()
        for line in lines:
            self.answer_raw_input.append(line.strip())

    def extract_values_from_answer(self):
        global answer_array
        m = 58
        n = 58
        o = 58
        p = 58

        # Creates the array
        answer_array = self.create_four_d_array(0, *(m, n, o, p))

        for x in range(0, len(self.answer_raw_input)):
            # Gets w
            par = self.answer_raw_input[x].partition(' ')
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
            temp_mo_value = par[0]

            # Set the value in the current iteration to index in 4D array
            answer_array[temp_w][temp_x][temp_y][temp_z] = temp_mo_value
            answer_array[temp_x][temp_w][temp_y][temp_z] = temp_mo_value
            answer_array[temp_x][temp_w][temp_z][temp_y] = temp_mo_value
            answer_array[temp_w][temp_x][temp_z][temp_y] = temp_mo_value

            answer_array[temp_y][temp_z][temp_w][temp_x] = temp_mo_value
            answer_array[temp_y][temp_z][temp_x][temp_w] = temp_mo_value
            answer_array[temp_z][temp_y][temp_x][temp_w] = temp_mo_value
            answer_array[temp_z][temp_y][temp_w][temp_x] = temp_mo_value

    # Testing to see if print works
    def print_f_raw_input(self, indexOne):
        print(self.f_raw_input[indexOne])

    # Creates the 4d array to store the AOs in
    def create_four_d_array(self, value, *dim):
        """
        Create 4D-array
        :param dim: a tuple of dimensions - (w, x, y, z)
        :param value: value with which 4D-array is to be filled
        :return: 4D-array
        """
        return [[[[value for x in range(dim[3])] for x in range(dim[2])] for x in range(dim[1])] for x in range(dim[0])]

    def find_d(self, i, j, a, b):
        return (float(f_array[i]) + float(f_array[j]) - float(f_array[a]) - float(f_array[b]))

    def find_t(self, a, b, i, j):
        return (2 * float(answer_array[i][a][j][b]) - float(answer_array[i][b][j][a]))/(float(self.find_d(a, b, i, j)))

    def summation(self):
        # global final_answer
        final_answer = 0
        for i in range(0,5):
            for j in range(0,5):
                for a in range(5, 58):
                    for b in range (5, 58):
                        final_answer = float(final_answer) + ((float(self.find_t(i, j, a, b)) * float(answer_array[i][a][j][b])))
        print("Final Answer:", final_answer)

start = time.time()
test = ProjectThree()
test.read_f('OrbEngsV2.dat')
test.create_1d_array()
test.extract_values_from_f()
test.read_answer("answer.dat")
test.extract_values_from_answer()
test.summation()
