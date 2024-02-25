import numpy as np
from datetime import datetime

FMT = '%H:%M'


class User:
    def __init__(self, name):
        self.name = name

    def get_user(self):
        return self.name


class Booking:
    def __init__(self, user, start_time: str, end_time: str, equipment: str):
        self.user = user
        self.start_time = start_time
        self.end_time = end_time
        self.equipment = equipment

    def is_intersect(self, another_booking):
        start1 = datetime.strptime(another_booking.start_time, FMT)
        end1 = datetime.strptime(another_booking.end_time, FMT)
        start2 = datetime.strptime(self.start_time, FMT)
        end2 = datetime.strptime(self.end_time, FMT)
        if self.equipment != another_booking.equipment:
            return False
        if ((end1 <= start2) & (start1 < start2)) or ((start1 >= end2) & (end1 > end2)):
            return False
        else:
            return True


class LabEquipment:
    def __init__(self, equipment: str, booking_history: list):
        self.equipment = equipment
        self.booking_history = booking_history

    def get_user_with_intercept(self, new_booking):
        for booking in self.booking_history:
            if booking.is_intersect(new_booking):
                return booking.user

    def is_available(self, equip_type, start_time, end_time):
        potential_booking = Booking(user=None, equipment=equip_type, start_time=start_time, end_time=end_time)
        non_availability = []
        for booking in self.booking_history:
            non_availability.append(booking.is_intersect(potential_booking))
        if sum(non_availability) > 0:
            return False
        else:
            return True

    def book(self, new_booking):
        if not self.is_available(equip_type=new_booking.equipment, start_time=new_booking.start_time,
                                 end_time=new_booking.end_time):
            collegue = self.get_user_with_intercept(new_booking)
            return f'Sorry, another booking has already been made during this time by {collegue}'
        else:
            self.booking_history.append(new_booking)
            return 'Successful booking'


class GenCodeInterpreter:
    def __init__(self):
        self.memory = np.zeros(5000)
        self.cursor = 0
        self.buffer = ''

    def eval(self, code: str):
        self.cursor = 0
        self.memory = np.zeros(5000)
        self.buffer = ''
        encoded = {'A': + 1, 'T': -1, 'G': +1, 'C': -1}
        for char in code:
            if char in 'AT':
                self.cursor += encoded[char]
            if char in 'CG':
                self.memory[self.cursor] += encoded[char]
            if char == 'N':
                self.buffer += str(round(self.memory[self.cursor])) + ','
                #print(self.buffer)
            elif char not in 'ATCGN':
                print(
                    "Your genetic code can not be interpreted. Make sure it contains only values in range: A, T, C, "
                    "G, N ")
            elif self.cursor > 5000:
                print('You run out of memory. There are only 5000 memory units')

        decoded = [chr(int(n)) for n in self.buffer.split(',')[:-1]]

        return ''.join(decoded)


def meet_the_dunders():
    res = 0
    matrix = []
    
    iter_l = range(0, 100, 10).__iter__()
    while True:
        try:
            idx = (iter_l).__next__()
            matrix += [list(range(idx, (idx).__add__(10)))]
        except:
            break

    def func_1(x):
        return x in range(1, 5, 2).__iter__()

    def func_2(x):
        return [x.__getitem__(col) for col in selected_columns_indices]

    selected_columns_indices = list(filter(func_1, range((matrix).__len__())))
    selected_columns = map(func_2, matrix)

    arr = np.array(list(selected_columns))

    mask = (arr[:, 1].__divmod__(3)[1]).__eq__(0)
    new_arr = arr[mask]

    product = new_arr.__matmul__(new_arr.T)

    if (product[0].__lt__(1000)).all() and (product[2].__gt__(1000)).any():
        res = int(((np.mean.__call__(product) // 10).__divmod__(100)[1]))
    return res
